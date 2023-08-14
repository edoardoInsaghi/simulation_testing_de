
# compute_general_metrics = function(results, de.genes, beta = 0.5, lfc_cut = 0, p_value_cut = 0.05) {
#
#   res <- results %>%
#     dplyr::filter(p_val_adj <= p_value_cut) %>%
#     dplyr::filter(abs(avg_log2FC) >= lfc_cut)
#
#   precision <- compute_precision(rownames(results), de.genes)
#   recall <- compute_recall(rownames(results), de.genes)
#   f_beta_score <- compute_fbeta_score(rownames(results), de.genes, beta = beta)
#
#   return(list(
#     precision = precision,
#     recall = recall,
#     f_beta_score = f_beta_score
#   ))
# }

get_roc_curve <- function(d, de.genes) {
  assertthat::assert_that(all(colnames(d) == c("gene", "f")))

  d <- d %>%
    arrange(f) %>%
    mutate(is_positive = gene %in% de.genes)

  P <- length(de.genes)
  N <- nrow(d) - P

  L_sorted <- d$is_positive
  f <- d$f
  FP <- 0
  TP <- 0
  FPR <- c()
  TPR <- c()
  f_prev <- -Inf
  i <- 1

  while (i <= length(L_sorted)) {
    if (f[i] != f_prev) {
      FPR <- c(FPR, FP / N)
      TPR <- c(TPR, TP / P)
      f_prev <- f[i]
    }

    if (L_sorted[i]) {
      TP <- TP + 1
    } else {
      FP <- FP + 1
    }

    i <- i + 1
  }

  FPR <- c(FPR, FP / N)
  TPR <- c(TPR, TP / P)

  dplyr::tibble(x=FPR, y=TPR) %>% stats::na.omit()
}


get_precision_recall_curve <- function(d, de.genes) {
  assertthat::are_equal(colnames(d), c("gene", "f"))

  d <- d %>%
    arrange(f) %>%
    mutate(is_positive = gene %in% de.genes)

  P <- length(de.genes)
  N <- nrow(d) - P

  L_sorted <- d$is_positive
  f <- d$f
  PP <- 0
  TP <- 0
  PPV <- c(1)
  TPR <- c(0)
  f_prev <- -Inf
  i <- 1

  while (i <= length(L_sorted)) {
    if (f[i] != f_prev) {
      PPV <- c(PPV, TP / PP)
      TPR <- c(TPR, TP / P)
      f_prev <- f[i]
    }

    if (L_sorted[i]) TP <- TP + 1

    PP <- PP + 1

    i <- i + 1
  }

  # FPR <- c(FPR, FP / N)
  PPV <- c(PPV, TP / PP)
  TPR <- c(TPR, TP / P)

  dplyr::tibble(x=TPR, y=PPV) %>% stats::na.omit()
}

approximate_AUC <- function(curve_df) {
  assertthat::are_equal(c("x","y"), colnames(curve_df))
  x <- curve_df$x
  y <- curve_df$y

  # Calculate the width of each trapezoid
  width <- diff(x)

  # Calculate the average height of adjacent Y values
  avg_height <- (y[-1] + y[-length(y)]) / 2

  # Calculate the area of each trapezoid
  trapezoid_area <- width * avg_height

  # Sum the areas of all trapezoids to get the approximate AUC
  auc <- sum(trapezoid_area)

  return(auc)
}

# compute_precision = function(pred, true) {
#   tp <- length( intersect(pred, true) )
#   fp <- length( setdiff(pred, true) )
#   tp / (tp + fp)
# }
#
# compute_recall = function(pred, true) {
#   tp <- length( intersect(pred, true) )
#   fn <- length( setdiff(true, pred) )
#   tp / (tp + fn)
# }
#
# compute_fbeta_score = function(pred, true, beta) {
#   precision = compute_precision(pred, true)
#   recall = compute_recall(pred, true)
#   (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall)
# }

