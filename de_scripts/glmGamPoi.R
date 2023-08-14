require(glmGamPoi)
require(magrittr)

glm_fit <- function(counts, model_matrix){
  glm.fit <- glmGamPoi::glm_gp(counts %>% as.matrix(),
                               design = model_matrix,
                               on_disk = F,
                               overdispersion = T)
  return(glm.fit)
}

glm_test <- function(fit, contrast = c(1, -1)){

  glm.results <- glmGamPoi::test_de(fit, contrast = contrast)
  colnames(glm.results) <- c("gene", "p_val", "p_val_adj", "f_statistics", "df1", "df2", "avg_log2FC")
  rownames(glm.results) <- glm.results$gene

  return(glm.results %>% dplyr::select("gene", "p_val_adj", "avg_log2FC"))
}

glm_ft <- function(counts, model_matrix, contrast = c(1, -1)){

  fit <- glm_fit(counts, model_matrix)
  results <- glm_test(fit, contrast)

  return(results)
}
