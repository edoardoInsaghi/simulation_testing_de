library(rdevil)
require(magrittr)

devil_fit <- function(counts, model_matrix, group_matrix = NULL){

  UMI <- colSums(sim$counts) / mean(colSums(sim$counts))

  devil.fit <- rdevil::fit_linear_model(
    input_matrix = counts,
    model_matrix = model_matrix,
    group_matrix = group_matrix,
    ncounts = UMI,
    gene_names = rownames(counts),
    method_specific_args = list(
      optimizer_name = "ClippedAdam",
      steps = as.integer(500),
      lr = 0.5,
      gamma_lr = 0.0001,
      cuda = TRUE,
      jit_compile = FALSE,
      batch_size = 10000L,
      full_cov = TRUE,
      prior_loc = 10,
      theta_bounds = c(0., 1e16),
      init_loc = 10,
      threshold = 0
    )
  )
  return(devil.fit)
}

devil_test <- function(fit, contrast = as.array(c(1, -1)), ROPE = F){

  if(ROPE){
    results <- rdevil::test_posterior_ROPE(devil.fit, contrast = contrast, LFC = 1e-10) %>%
      dplyr::rename(avg_log2FC = log_FC, p_val_adj = ROPE)

    rownames(results) <- results$gene

    return(results)
  }

  results <- rdevil::test_posterior_null(fit, contrast = contrast) %>%
    dplyr::rename(avg_log2FC = log_FC, p_val_adj = p_value_adj)

  rownames(results) <- results$gene

  return(results %>% dplyr::select("gene", "p_val_adj", "avg_log2FC"))
}

devil_ft <- function(counts, model_matrix, contrast = as.array(c(1, -1)), ROPE = F){

  fit <- devil_fit(counts, model_matrix)
  results <- devil_test(fit, contrast)

  return(results)
}
