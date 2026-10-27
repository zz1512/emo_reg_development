# Age decomposition at the paper's basis size, k = 3.
#
# The cubic-regression basis s(age, k = 3) already contains a straight line.
# Putting age and s(age, k = 3) in one formula aliases the intercept, so the
# two p-values do not separate the linear and nonlinear parts.
#
# This script keeps that k = 3 basis and splits it once:
#   f(age) = beta * age + gamma * h(age)
# h(age) is the single direction in the basis orthogonal to the intercept and
# to age. Its cubic-regression penalty is passed to mgcv through paraPen, and
# the smoothing parameter is chosen by REML. Inference is the Wald p-value of
# beta and of gamma. No anova() comparison is performed.
#
# Linear covariates are gender, site_id, parent_edu, and meanFD.
# parent_edu is 父母受教育水平_父母受教育水平, aligned by 被试编号.
# meanFD is the mean of EMO and SOC mean_fd in qc_res_subs_runs_fmri.csv, aligned by sub_id.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 || !args[[3]] %in% c("terms", "matrix")) {
  stop("Usage: Rscript mgcv_age_decomposition.R <input_dir> <output_dir> <terms|matrix>")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
mode <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(mgcv))

design <- read.csv(file.path(input_dir, "design.csv"), stringsAsFactors = FALSE)
covariates <- c("gender", "site_id", "parent_edu", "meanFD")
needed <- c("sub_id", "age", covariates)
# Intercept and linear age, plus each covariate, are unpenalized.
n_unpenalized <- 2 + length(covariates)
if (any(!needed %in% names(design))) {
  stop("design.csv is missing one of: ", paste(needed, collapse = ", "))
}

age <- as.numeric(design$age)
if (any(!is.finite(age))) stop("age contains non-finite values")
age_ref <- mean(age)
age_grid <- seq(min(age), max(age), length.out = 100)

# One shared nonlinear direction, evaluated at the observed ages and on the grid.
basis <- smoothCon(s(age, k = 3, bs = "cr"), data = data.frame(age = age),
                   knots = NULL, absorb.cons = TRUE)[[1]]
raw_basis <- basis$X
residual_basis <- qr.resid(qr(cbind(1, age)), raw_basis)
h_qr <- qr(residual_basis)
if (h_qr$rank < 1) stop("k = 3 basis has no direction orthogonal to linear age")
h <- qr.Q(h_qr)[, 1]
loading <- qr.coef(qr(raw_basis), h)
loading[is.na(loading)] <- 0
reconstruction <- as.numeric(raw_basis %*% loading)
if (max(abs(reconstruction - h)) > 1e-6) {
  stop("Nonlinear direction is not inside the k = 3 basis")
}
penalty <- as.numeric(t(loading) %*% basis$S[[1]] %*% loading)
if (!is.finite(penalty) || penalty < 0) stop("Cubic-regression penalty of h(age) is not usable")
h_grid <- as.numeric(PredictMat(basis, data.frame(age = age_grid)) %*% loading)
h_ref <- as.numeric(PredictMat(basis, data.frame(age = age_ref)) %*% loading)

model_data <- data.frame(
  age = age,
  h = h,
  gender = as.numeric(design$gender),
  site_id = as.numeric(design$site_id),
  parent_edu = as.numeric(design$parent_edu),
  meanFD = as.numeric(design$meanFD)
)
full_formula <- as.formula(paste("y ~ age + h +", paste(covariates, collapse = " + ")))
reduced_formula <- as.formula(paste("y ~", paste(covariates, collapse = " + ")))
if (!all(is.finite(as.matrix(model_data)))) {
  stop("Age model data contain non-finite covariate values")
}

shape_of <- function(curve) {
  change <- diff(curve)
  net <- curve[[length(curve)]] - curve[[1]]
  tolerance <- 1e-4 * max(1, diff(range(curve)))
  if (abs(net) <= tolerance) return("flat")
  if (all(change >= -tolerance) && net > 0) return("increasing")
  if (all(change <= tolerance) && net < 0) return("decreasing")
  "nonmonotone"
}

fit_one <- function(y) {
  fit_data <- model_data
  fit_data$y <- as.numeric(y)
  if (any(!is.finite(fit_data$y))) stop("Response contains non-finite values")
  full <- gam(full_formula, data = fit_data, method = "REML",
              paraPen = list(h = list(S = matrix(penalty, 1, 1))))
  covariates_only <- lm(reduced_formula, data = fit_data)
  table <- summary(full)$p.table
  rss_full <- sum(residuals(full)^2)
  rss_covariates <- sum(residuals(covariates_only)^2)
  partial_r2 <- if (rss_covariates > 0) max(0, 1 - rss_full / rss_covariates) else NA_real_
  beta <- unname(coef(full)[["age"]])
  gamma <- unname(coef(full)[["h"]])
  centered <- beta * (age_grid - age_ref) + gamma * (h_grid - h_ref)
  net <- centered[[length(centered)]] - centered[[1]]
  sign_net <- if (abs(net) <= 1e-12) 0 else sign(net)
  data.frame(
    beta_age = beta,
    se_age = unname(table["age", "Std. Error"]),
    p_linear = unname(table["age", "Pr(>|t|)"]),
    gamma_h = gamma,
    se_h = unname(table["h", "Std. Error"]),
    p_smooth = unname(table["h", "Pr(>|t|)"]),
    edf_h = sum(full$edf) - n_unpenalized,
    edf_model = sum(full$edf),
    sp_h = unname(full$sp[[1]]),
    partial_r2 = partial_r2,
    signed_partial_r2 = sign_net * partial_r2,
    net_change = net,
    shape = shape_of(centered),
    stringsAsFactors = FALSE
  )
}

curve_of <- function(stats_row) {
  stats_row$beta_age * (age_grid - age_ref) + stats_row$gamma_h * (h_grid - h_ref)
}

add_bh <- function(frame, columns) {
  for (column in columns) {
    q_column <- sub("^p_", "q_", column)
    frame[[q_column]] <- NA_real_
    ok <- is.finite(frame[[column]])
    if (any(ok)) frame[[q_column]][ok] <- p.adjust(frame[[column]][ok], method = "BH")
  }
  frame
}

if (mode == "terms") {
parcel_y <- read.csv(file.path(input_dir, "parcel_response.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(as.character(design$sub_id), as.character(parcel_y$sub_id))) {
  stop("design.csv and parcel_response.csv must have the same subject order")
}
parcel_ids <- names(parcel_y)[names(parcel_y) != "sub_id"]
parcel_stats <- vector("list", length(parcel_ids))
parcel_curves <- vector("list", length(parcel_ids))
for (index in seq_along(parcel_ids)) {
  if (index %% 25 == 0 || index == length(parcel_ids)) {
    message("Parcel ", index, " / ", length(parcel_ids))
  }
  parcel <- parcel_ids[[index]]
  stats_row <- fit_one(parcel_y[[parcel]])
  stats_row$parcel <- as.integer(parcel)
  parcel_stats[[index]] <- stats_row
  parcel_curves[[index]] <- data.frame(
    parcel = as.integer(parcel),
    age = age_grid,
    centered_effect = curve_of(stats_row)
  )
}
parcel_stats <- add_bh(do.call(rbind, parcel_stats), c("p_linear", "p_smooth"))
parcel_curves <- do.call(rbind, parcel_curves)
write.csv(parcel_stats, file.path(output_dir, "parcel_gam_statistics.csv"), row.names = FALSE)
write.csv(parcel_curves, file.path(output_dir, "parcel_curves.csv"), row.names = FALSE)

group_y <- read.csv(file.path(input_dir, "group_response.csv"), stringsAsFactors = FALSE)
yeo_order <- read.csv(file.path(input_dir, "yeo_order.csv"), stringsAsFactors = FALSE)$group
system_order <- read.csv(file.path(input_dir, "system_order.csv"), stringsAsFactors = FALSE)$group

fit_group_means <- function(dimension, group_order) {
  subset_y <- group_y[group_y$dimension == dimension, , drop = FALSE]
  wide <- reshape(subset_y, idvar = "sub_id", timevar = "group", direction = "wide")
  names(wide) <- sub("^y\\.", "", names(wide))
  wide <- wide[match(as.character(design$sub_id), as.character(wide$sub_id)), , drop = FALSE]
  missing_groups <- setdiff(group_order, names(wide))
  if (length(missing_groups) > 0) {
    stop(dimension, " is missing groups: ", paste(missing_groups, collapse = ", "))
  }
  stats_rows <- vector("list", length(group_order))
  curve_rows <- vector("list", length(group_order))
  for (index in seq_along(group_order)) {
    group <- group_order[[index]]
    stats_row <- fit_one(wide[[group]])
    stats_row$dimension <- dimension
    stats_row$group <- group
    stats_rows[[index]] <- stats_row
    curve_rows[[index]] <- data.frame(
      dimension = dimension,
      group = group,
      age = age_grid,
      centered_effect = curve_of(stats_row),
      stringsAsFactors = FALSE
    )
  }
  list(
    stats = add_bh(do.call(rbind, stats_rows), c("p_linear", "p_smooth")),
    curves = do.call(rbind, curve_rows)
  )
}

yeo_terms <- fit_group_means("yeo", yeo_order)
system_terms <- fit_group_means("system", system_order)
write.csv(rbind(yeo_terms$stats, system_terms$stats),
          file.path(output_dir, "group_gam_statistics.csv"), row.names = FALSE)
write.csv(rbind(yeo_terms$curves, system_terms$curves),
          file.path(output_dir, "group_curves.csv"), row.names = FALSE)
write.csv(data.frame(age_reference = age_ref, basis_k = 3, basis = "cr", step = "terms"),
          file.path(output_dir, "model_reference.csv"), row.names = FALSE)
message("linear and smooth terms wrote ", output_dir)
quit(save = "no", status = 0)
}

group_y <- read.csv(file.path(input_dir, "group_response.csv"), stringsAsFactors = FALSE)
yeo_order <- read.csv(file.path(input_dir, "yeo_order.csv"), stringsAsFactors = FALSE)$group
system_order <- read.csv(file.path(input_dir, "system_order.csv"), stringsAsFactors = FALSE)$group

fit_dimension <- function(dimension, group_order) {
  subset_y <- group_y[group_y$dimension == dimension, , drop = FALSE]
  wide <- reshape(subset_y, idvar = "sub_id", timevar = "group", direction = "wide")
  names(wide) <- sub("^y\\.", "", names(wide))
  if (!identical(as.character(wide$sub_id), as.character(design$sub_id))) {
    wide <- wide[match(as.character(design$sub_id), as.character(wide$sub_id)), , drop = FALSE]
  }
  missing_groups <- setdiff(group_order, names(wide))
  if (length(missing_groups) > 0) {
    stop(dimension, " is missing groups: ", paste(missing_groups, collapse = ", "))
  }
  pair_stats <- list()
  pair_curves <- list()
  pair_index <- 0
  for (row_name in group_order) {
    for (column_name in group_order) {
      if (identical(row_name, column_name)) next
      pair_index <- pair_index + 1
      stats_row <- fit_one(wide[[row_name]] - wide[[column_name]])
      stats_row$dimension <- dimension
      stats_row$group_a <- row_name
      stats_row$group_b <- column_name
      pair_stats[[pair_index]] <- stats_row
      pair_curves[[pair_index]] <- data.frame(
        dimension = dimension,
        group_a = row_name,
        group_b = column_name,
        age = age_grid,
        centered_effect = curve_of(stats_row),
        stringsAsFactors = FALSE
      )
    }
  }
  stats_frame <- add_bh(do.call(rbind, pair_stats), c("p_linear", "p_smooth"))
  list(stats = stats_frame, curves = do.call(rbind, pair_curves))
}

yeo_fit <- fit_dimension("yeo", yeo_order)
system_fit <- fit_dimension("system", system_order)

parcel_y <- read.csv(file.path(input_dir, "parcel_response.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
parcel_ids <- names(parcel_y)[names(parcel_y) != "sub_id"]
parcel_matrix <- as.matrix(parcel_y[, parcel_ids, drop = FALSE])
storage.mode(parcel_matrix) <- "double"
n_parcels <- length(parcel_ids)
pair_n <- n_parcels * (n_parcels - 1)
message("Parcel pairs ", pair_n)
parcel_chunks <- list()
chunk <- vector("list", 400)
chunk_index <- 0
done <- 0
for (row_index in seq_len(n_parcels)) {
  if (row_index %% 10 == 0 || row_index == n_parcels) {
    message("Parcel pair row ", row_index, " / ", n_parcels)
  }
  outcome_row <- parcel_matrix[, row_index]
  for (column_index in seq_len(n_parcels)) {
    if (row_index == column_index) next
    done <- done + 1
    stats_row <- fit_one(outcome_row - parcel_matrix[, column_index])
    stats_row$dimension <- "parcel"
    stats_row$group_a <- parcel_ids[[row_index]]
    stats_row$group_b <- parcel_ids[[column_index]]
    chunk_index <- chunk_index + 1
    chunk[[chunk_index]] <- stats_row
    if (chunk_index == length(chunk)) {
      parcel_chunks[[length(parcel_chunks) + 1]] <- do.call(rbind, chunk)
      chunk_index <- 0
    }
  }
}
if (chunk_index > 0) {
  parcel_chunks[[length(parcel_chunks) + 1]] <- do.call(rbind, chunk[seq_len(chunk_index)])
}
parcel_pairs <- add_bh(do.call(rbind, parcel_chunks), c("p_linear", "p_smooth"))
if (nrow(parcel_pairs) != pair_n) {
  stop("Parcel pair count is ", nrow(parcel_pairs), ", expected ", pair_n)
}

difference_stats <- rbind(yeo_fit$stats, system_fit$stats, parcel_pairs)
difference_curves <- rbind(yeo_fit$curves, system_fit$curves)

write.csv(difference_stats, file.path(output_dir, "difference_statistics.csv"), row.names = FALSE)
write.csv(difference_curves, file.path(output_dir, "difference_curves.csv"), row.names = FALSE)
write.csv(data.frame(age_reference = age_ref, basis_k = 3, basis = "cr", step = "matrix"),
          file.path(output_dir, "model_reference.csv"), row.names = FALSE)
message("comparison matrix wrote ", output_dir)
