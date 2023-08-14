rm(list = ls())
require(dplyr)

source("de_scripts/devil.R")
source("de_scripts/glmGamPoi.R")
source("de_scripts/MAST.R")
source("de_scripts/wilcoxon.R")
source("utils/metrics.R")

# Path to de_results
path_de_res <- "de_results/"

# Path to general results
path_general_res <- "general_results/"

# Leggi lista dei files
path_to_sim_folder <- "simulations/"
sim_paths <- list.files(path_to_sim_folder, full.names = T)

for (p in sim_paths) {
  sim <- readRDS(p)

  # Prepare counts
  X <- sim$counts.data %>% as.matrix()
  model_matrix <- model.matrix(~ Group - 1, sim$cell.data)
  group_matrix <- model.matrix(~ Batch - 1, sim$cell.data)
  UMI <- colSums(X) / mean(colSums(X))
  contrast = as.array(c(1,-1))

  # Fit data
  # Devil
  s <- Sys.time()
  devil.fit <- devil_fit(X, model_matrix)
  e <- Sys.time()
  time.devil.fit <- e - s

  s <- Sys.time()
  devil.res <- devil_test(devil.fit, contrast) %>% dplyr::mutate(name = "devil")
  e <- Sys.time()
  time.devil.test <- e - s

  s <- Sys.time()
  devil.rope.res <- devil_test(devil.fit, contrast, ROPE = T) %>% dplyr::mutate(name = "devil rope")
  e <- Sys.time()
  time.devil.test.rope <- e - s

  # Glm
  s <- Sys.time()
  glm.fit <- glm_fit(X, model_matrix)
  e <- Sys.time()
  time.glm.fit <- e - s

  s <- Sys.time()
  glm.res <- glm_test(glm.fit, contrast) %>% dplyr::mutate(name = "glm")
  e <- Sys.time()
  time.glm.test <- e - s

  # Mast
  s <- Sys.time()
  mast.fit <- mast_fit(sim, stats::formula(~ Group))
  e <- Sys.time()
  time.mast.fit <- e - s

  s <- Sys.time()
  mast.res <- mast_test(mast.fit, "GroupGroup2") %>% dplyr::mutate(name = 'mast')
  e <- Sys.time()
  time.mast.test <- e - s

  # Results
  results <- list(
    "glm" = glm.res,
    "devil" = devil.res,
    "devil rope" = devil.rope.res,
    "mast" = mast.res
  )

  # Prepare time dataframe
  timings <- dplyr::tibble(name = names(results),
                           time_fit = c(time.glm.fit, time.devil.fit, time.devil.fit, time.mast.fit),
                           time_test = c(time.glm.test, time.devil.test, time.devil.test.rope, time.mast.test)
                           ) %>%
    dplyr::mutate(time_fit = as.numeric(time_fit), time_test = as.numeric(time_test)) %>%
    dplyr::mutate(time_total = time_fit + time_test)

  # Compute AUPR and AUROC for all models

  pr_curves_df <- lapply(names(results), function(n) {
    res <- results[[n]] %>% dplyr::select(gene, p_val_adj) %>% `colnames<-`(c("gene", "f")) %>% na.omit()
    get_precision_recall_curve(res, de.genes = sim$de.genes) %>% dplyr::mutate(name = n)
  }) %>% do.call("bind_rows", .)

  roc_curves_df <- lapply(names(results), function(n) {
    res <- results[[n]] %>% dplyr::select(gene, p_val_adj) %>% `colnames<-`(c("gene", "f")) %>% na.omit()
    get_roc_curve(res, de.genes = sim$de.genes) %>%
      dplyr::mutate(name = n)
  }) %>% do.call("bind_rows", .)

  AUPRs <- lapply(pr_curves_df$name %>% unique(), function(n) {
    dplyr::tibble(
      aupr = approximate_AUC(pr_curves_df %>% dplyr::filter(name == n)),
      name = n
    )
  }) %>% do.call("rbind", .) %>% arrange(- aupr)
  AUPRs <- AUPRs %>% dplyr::mutate(aupr = round(aupr, 3))

  AUROCs <- lapply(roc_curves_df$name %>% unique(), function(n) {
    dplyr::tibble(
      auroc = approximate_AUC(roc_curves_df %>% dplyr::filter(name == n)),
      name = n
    )
  }) %>% do.call("bind_rows", .) %>% arrange(- auroc)
  AUROCs <- AUROCs %>% dplyr::mutate(auroc = round(auroc, 3))

  # Merge general results
  full_table <- dplyr::full_join(AUPRs, AUROCs, by='name') %>% dplyr::select(name, auroc, aupr) %>% `rownames<-`(NULL)
  full_table <- dplyr::full_join(full_table, timings, by='name') %>% `rownames<-`(NULL)

  # Merge de Results
  de_results <- do.call("bind_rows", results)

  # Save
  nn <- strsplit(strsplit(p, "//")[[1]][2], ".rds")[[1]][1]
  de_name <- paste0(path_de_res, "de_results_", nn, ".rds")
  res_name <- paste0(path_general_res, "general_results_", nn, ".rds")

  saveRDS(full_table, res_name)
  saveRDS(de_results, de_name)
}

