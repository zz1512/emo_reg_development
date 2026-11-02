# Nested age models.
#
# m1: covariates only
# m2: covariates + linear age
# m3: covariates + s(age, k = basis_k, bs = "cr"); basis_k defaults to 3
#
# All three are mgcv::gam fits. method is GCV.Cp or REML. The reported tests are
# nested F comparisons from anova():
#   m2 vs m1  linear age
#   m3 vs m1  any age effect
#   m3 vs m2  smooth beyond a straight line
# BH correction is applied inside each dimension, and separately for each test.
# The summary p-value of s(age) is not used. It includes the linear part.

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
  stop("Usage: Rscript mgcv_nested_models.R <input_dir> <output_dir> [GCV|REML] [k] [cr|ps]")
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
age_ref <- mean(age)
age_grid <- seq(min(age), max(age), length.out = 100)
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

formula_m1 <- y ~ gender + site_id + parent_edu + meanFD
formula_m2 <- y ~ age + gender + site_id + parent_edu + meanFD
formula_m3 <- as.formula(sprintf(
  "y ~ s(age, k = %d, bs = \"%s\") + gender + site_id + parent_edu + meanFD",
  basis_k, smooth_basis
))

nested_test <- function(smaller, larger) {
  tab <- anova(smaller, larger, test = "F")
  p_col <- if ("Pr(>F)" %in% names(tab)) "Pr(>F)" else if ("Pr(>Chi)" %in% names(tab)) "Pr(>Chi)" else ""
  stat_col <- if ("F" %in% names(tab)) "F" else if ("Chisq" %in% names(tab)) "Chisq" else ""
  if (!nzchar(p_col) || !nzchar(stat_col) || !"Df" %in% names(tab)) {
    stop("anova.gam returned columns: ", paste(names(tab), collapse = ", "))
  }
  row <- nrow(tab)
  statistic <- unname(tab[[stat_col]][row])
  df <- unname(tab[["Df"]][row])
  p_value <- unname(tab[[p_col]][row])
  if (!is.finite(p_value) || !is.finite(statistic)) {
    stop("Nested comparison did not return a finite statistic")
  }
  c(statistic = statistic, df = df, p = p_value)
}

fit_three <- function(y) {
  fit_data <- model_data
  fit_data$y <- as.numeric(y)
  if (any(!is.finite(fit_data$y))) stop("Response contains non-finite values")
  m1 <- gam(formula_m1, data = fit_data, method = smooth_method)
  m2 <- gam(formula_m2, data = fit_data, method = smooth_method)
  m3 <- gam(formula_m3, data = fit_data, method = smooth_method)
  linear <- nested_test(m1, m2)
  any_age <- nested_test(m1, m3)
  beyond <- nested_test(m2, m3)
  age_row <- summary(m2)$p.table["age", , drop = TRUE]
  m3_at_grid <- as.numeric(predict(m3, newdata = prediction_grid))
  m3_at_ref <- as.numeric(predict(m3, newdata = reference_row))
  stats <- data.frame(
    n = length(fit_data$y),
    beta_age = unname(age_row[["Estimate"]]),
    se_age = unname(age_row[["Std. Error"]]),
    edf_m1 = sum(m1$edf),
    edf_m2 = sum(m2$edf),
    edf_m3 = sum(m3$edf),
    deviance_m1 = deviance(m1),
    deviance_m2 = deviance(m2),
    deviance_m3 = deviance(m3),
    F_m2_vs_m1 = unname(linear[["statistic"]]),
    df_m2_vs_m1 = unname(linear[["df"]]),
    p_m2_vs_m1 = unname(linear[["p"]]),
    F_m3_vs_m1 = unname(any_age[["statistic"]]),
    df_m3_vs_m1 = unname(any_age[["df"]]),
    p_m3_vs_m1 = unname(any_age[["p"]]),
    F_m3_vs_m2 = unname(beyond[["statistic"]]),
    df_m3_vs_m2 = unname(beyond[["df"]]),
    p_m3_vs_m2 = unname(beyond[["p"]]),
    stringsAsFactors = FALSE
  )
  curves <- rbind(
    data.frame(model = "m2", age = age_grid,
               centered_effect = unname(stats$beta_age) * (age_grid - age_ref),
               stringsAsFactors = FALSE),
    data.frame(model = "m3", age = age_grid,
               centered_effect = m3_at_grid - m3_at_ref,
               stringsAsFactors = FALSE)
  )
  list(stats = stats, curves = curves)
}

add_bh <- function(frame) {
  for (column in c("p_m2_vs_m1", "p_m3_vs_m1", "p_m3_vs_m2")) {
    q_column <- sub("^p_", "q_", column)
    frame[[q_column]] <- NA_real_
    ok <- is.finite(frame[[column]])
    if (!all(ok)) stop(column, " contains non-finite p-values")
    frame[[q_column]] <- p.adjust(frame[[column]], method = "BH")
  }
  frame
}

fit_columns <- function(dimension, item_names, y_matrix) {
  stats_rows <- vector("list", length(item_names))
  curve_rows <- vector("list", length(item_names))
  for (index in seq_along(item_names)) {
    item <- item_names[[index]]
    if (index %% 25 == 0 || index == length(item_names) || length(item_names) < 25) {
      message(dimension, " ", index, " / ", length(item_names), " ", item)
    }
    fitted <- fit_three(y_matrix[, index])
    stats_row <- fitted$stats
    stats_row$dimension <- dimension
    stats_row$name <- item
    stats_rows[[index]] <- stats_row
    curve_rows[[index]] <- data.frame(
      dimension = dimension,
      name = item,
      fitted$curves,
      stringsAsFactors = FALSE
    )
  }
  list(stats = add_bh(do.call(rbind, stats_rows)),
       curves = do.call(rbind, curve_rows))
}

same_subjects <- function(response_ids, label) {
  if (!identical(as.character(design$sub_id), as.character(response_ids))) {
    stop(label, " subject order does not match design.csv")
  }
}

parcel_y <- read.csv(file.path(input_dir, "parcel_response.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
same_subjects(parcel_y$sub_id, "parcel_response.csv")
parcel_names <- names(parcel_y)[names(parcel_y) != "sub_id"]
parcel_matrix <- as.matrix(parcel_y[, parcel_names, drop = FALSE])
storage.mode(parcel_matrix) <- "double"
parcel_fit <- fit_columns("parcel", parcel_names, parcel_matrix)

read_group <- function(dimension, order_file) {
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
  fit_columns(dimension, order, values)
}

stats <- parcel_fit$stats
curves <- parcel_fit$curves
# Parcel-only inputs omit these files. The 7-network run still writes all three.
if (file.exists(file.path(input_dir, "yeo_order.csv"))) {
  yeo_fit <- read_group("yeo", "yeo_order.csv")
  stats <- rbind(stats, yeo_fit$stats)
  curves <- rbind(curves, yeo_fit$curves)
}
if (file.exists(file.path(input_dir, "system_order.csv"))) {
  system_fit <- read_group("system", "system_order.csv")
  stats <- rbind(stats, system_fit$stats)
  curves <- rbind(curves, system_fit$curves)
}
if (file.exists(file.path(input_dir, "amygdala_response.csv"))) {
  amygdala_y <- read.csv(file.path(input_dir, "amygdala_response.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  same_subjects(amygdala_y$sub_id, "amygdala_response.csv")
  amygdala_names <- read.csv(file.path(input_dir, "amygdala_order.csv"), stringsAsFactors = FALSE)$roi
  missing_amygdala <- setdiff(amygdala_names, names(amygdala_y))
  if (length(missing_amygdala)) {
    stop("amygdala_response.csv is missing: ", paste(missing_amygdala, collapse = ", "))
  }
  if (length(amygdala_names) != 6) {
    stop("Expected 6 amygdala ROIs, found ", length(amygdala_names))
  }
  amygdala_matrix <- as.matrix(amygdala_y[, amygdala_names, drop = FALSE])
  storage.mode(amygdala_matrix) <- "double"
  amygdala_fit <- fit_columns("amygdala", amygdala_names, amygdala_matrix)
  stats <- rbind(stats, amygdala_fit$stats)
  curves <- rbind(curves, amygdala_fit$curves)
}
write.csv(stats, file.path(output_dir, "nested_model_comparison.csv"), row.names = FALSE)
write.csv(curves, file.path(output_dir, "nested_curves.csv"), row.names = FALSE)
write.csv(data.frame(
  age_reference = age_ref,
  basis_k = basis_k,
  basis = smooth_basis,
  method = smooth_method,
  comparison = "anova F: m2 vs m1, m3 vs m1, m3 vs m2",
  bh = "within dimension and within comparison"
), file.path(output_dir, "model_reference.csv"), row.names = FALSE)
message("nested models wrote ", output_dir)
