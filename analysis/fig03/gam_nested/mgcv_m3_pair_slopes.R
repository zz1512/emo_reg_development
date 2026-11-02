# Pairwise m3 slope comparison for Yeo networks and the four systems.
#
# For each ordered pair, y is the subject-level difference group A minus group B.
# The model is the same m3 used for the single-region fits.
# method is GCV.Cp or REML, matching the nested fit.
#   y ~ s(age, k = 3, bs = "cr") + gender + site_id + parent_edu + meanFD
# The slope is the average change per year of that m3 age curve,
# (f(age_max) - f(age_min)) / (age_max - age_min), with covariates held fixed.
# Its p-value is the Wald test of this contrast, using the penalized covariance.
# BH correction is done inside the network pairs and inside the system pairs.

resolve_method <- function(requested) {
  key <- toupper(requested)
  if (key %in% c("GCV", "GCV.CP")) return("GCV.Cp")
  if (key == "REML") return("REML")
  stop("method must be GCV or REML, got: ", requested)
}

resolve_basis <- function(requested) {
  key <- tolower(requested)
  if (key %in% c("cr", "ps")) return(key)
  stop("basis must be cr or ps, got: ", requested)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2, 3, 4, 5)) {
  stop("Usage: Rscript mgcv_m3_pair_slopes.R <input_dir> <output_dir> [GCV|REML] [k] [cr|ps]")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
smooth_method <- if (length(args) >= 3) resolve_method(args[[3]]) else "REML"
basis_k <- 3L
if (length(args) >= 4) {
  basis_k <- suppressWarnings(as.integer(args[[4]]))
  if (is.na(basis_k) || basis_k < 3) stop("k must be an integer >= 3, got: ", args[[4]])
}
smooth_basis <- if (length(args) >= 5) resolve_basis(args[[5]]) else "cr"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(mgcv))

design <- read.csv(file.path(input_dir, "design.csv"), stringsAsFactors = FALSE)
if ("in_analysis" %in% names(design)) {
  design <- design[!is.na(design$in_analysis) & as.integer(design$in_analysis) == 1, , drop = FALSE]
  rownames(design) <- NULL
}
covariates <- c("gender", "site_id", "parent_edu", "meanFD")
needed <- c("sub_id", "age", covariates)
if (any(!needed %in% names(design))) {
  stop("design.csv is missing one of: ", paste(needed, collapse = ", "))
}

model_data <- data.frame(
  age = as.numeric(design$age),
  gender = as.numeric(design$gender),
  site_id = as.numeric(design$site_id),
  parent_edu = as.numeric(design$parent_edu),
  meanFD = as.numeric(design$meanFD)
)
if (!all(is.finite(as.matrix(model_data)))) {
  stop("Age model data contain non-finite covariate values")
}

age <- model_data$age
age_min <- min(age)
age_max <- max(age)
age_span <- age_max - age_min
if (!is.finite(age_span) || age_span <= 0) stop("Age has no range")
age_ref <- mean(age)
age_grid <- seq(age_min, age_max, length.out = 100)
reference_row <- data.frame(
  age = age_ref,
  gender = mean(model_data$gender),
  site_id = mean(model_data$site_id),
  parent_edu = mean(model_data$parent_edu),
  meanFD = mean(model_data$meanFD)
)
prediction_grid <- reference_row[rep(1, length(age_grid)), , drop = FALSE]
prediction_grid$age <- age_grid
rownames(prediction_grid) <- NULL
endpoints <- reference_row[c(1, 1), , drop = FALSE]
endpoints$age <- c(age_min, age_max)
rownames(endpoints) <- NULL

formula_m3 <- as.formula(sprintf(
  "y ~ s(age, k = %d, bs = \"%s\") + gender + site_id + parent_edu + meanFD",
  basis_k, smooth_basis
))

fit_difference <- function(y) {
  fit_data <- model_data
  fit_data$y <- as.numeric(y)
  if (any(!is.finite(fit_data$y))) stop("A pairwise difference contains non-finite values")
  m3 <- gam(formula_m3, data = fit_data, method = smooth_method)
  contrast <- predict(m3, newdata = endpoints, type = "lpmatrix")
  linear_combination <- (contrast[2, , drop = TRUE] - contrast[1, , drop = TRUE]) / age_span
  coefficients <- coef(m3)
  covariance <- vcov(m3)
  slope <- as.numeric(sum(linear_combination * coefficients))
  se <- sqrt(as.numeric(t(linear_combination) %*% covariance %*% linear_combination))
  if (!is.finite(slope) || !is.finite(se) || se <= 0) {
    stop("The m3 slope contrast did not return a finite standard error")
  }
  z_value <- slope / se
  p_value <- 2 * pnorm(-abs(z_value))
  centered <- as.numeric(predict(m3, newdata = prediction_grid)) -
    as.numeric(predict(m3, newdata = reference_row))
  list(slope = slope, se = se, z = z_value, p = p_value, curve = centered)
}

wide_groups <- function(dimension, order_file) {
  group_y <- read.csv(file.path(input_dir, "group_response.csv"), stringsAsFactors = FALSE)
  order <- read.csv(file.path(input_dir, order_file), stringsAsFactors = FALSE)$group
  subset_y <- group_y[group_y$dimension == dimension, c("sub_id", "group", "y"), drop = FALSE]
  wide <- reshape(subset_y, idvar = "sub_id", timevar = "group", direction = "wide")
  names(wide) <- sub("^y\\.", "", names(wide))
  wide <- wide[match(as.character(design$sub_id), as.character(wide$sub_id)), , drop = FALSE]
  if (any(is.na(wide$sub_id)) || !identical(as.character(wide$sub_id), as.character(design$sub_id))) {
    stop(dimension, " subjects do not match design.csv")
  }
  missing <- setdiff(order, names(wide))
  if (length(missing)) stop(dimension, " responses are missing: ", paste(missing, collapse = ", "))
  values <- as.matrix(wide[, order, drop = FALSE])
  storage.mode(values) <- "double"
  list(order = order, values = values)
}

compare_dimension <- function(dimension, order_file) {
  prepared <- wide_groups(dimension, order_file)
  order <- prepared$order
  values <- prepared$values
  stats_rows <- list()
  curve_rows <- list()
  row_index <- 0
  for (left in seq_along(order)) {
    for (right in seq_along(order)) {
      if (left == right) next
      row_index <- row_index + 1
      message(dimension, " pair ", row_index, ": ", order[[left]], " - ", order[[right]])
      fitted <- fit_difference(values[, left] - values[, right])
      stats_rows[[row_index]] <- data.frame(
        dimension = dimension,
        group_a = order[[left]],
        group_b = order[[right]],
        slope = fitted$slope,
        se_slope = fitted$se,
        z_slope = fitted$z,
        p_slope = fitted$p,
        stringsAsFactors = FALSE
      )
      curve_rows[[row_index]] <- data.frame(
        dimension = dimension,
        group_a = order[[left]],
        group_b = order[[right]],
        age = age_grid,
        centered_effect = fitted$curve,
        stringsAsFactors = FALSE
      )
    }
  }
  stats <- do.call(rbind, stats_rows)
  stats$q_slope <- p.adjust(stats$p_slope, method = "BH")
  list(stats = stats, curves = do.call(rbind, curve_rows))
}

yeo <- compare_dimension("yeo", "yeo_order.csv")
system <- compare_dimension("system", "system_order.csv")
write.csv(rbind(yeo$stats, system$stats),
          file.path(output_dir, "m3_pair_slopes.csv"), row.names = FALSE)
write.csv(rbind(yeo$curves, system$curves),
          file.path(output_dir, "m3_pair_curves.csv"), row.names = FALSE)
message("m3 pairwise slopes wrote ", output_dir)
