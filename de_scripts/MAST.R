require(MAST)
require(magrittr)

mast_fit <- function(sim, formula){
  sca <- MAST::FromMatrix(sim$counts.data %>% as.matrix(), check_sanity = F)
  for (n in names(sim$cell.data)) {
    SummarizedExperiment::colData(sca)[n] <- sim$cell.data[[n]]
  }

  mast.fit <- MAST::zlm(formula, sca)
  return(mast.fit)
}

mast_test <- function(fit, contrast = "groupsGroup2"){

  summaryCond <- MAST::summary(fit, doLRT=contrast)
  summaryCond <- summaryCond$datatable

  fcHurdle <- merge(summaryCond[contrast==contrast & component=='H',.(primerid, `Pr(>Chisq)`)],
                    summaryCond[contrast==contrast & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')

  fcHurdle[,fdr := stats::p.adjust(`Pr(>Chisq)`, 'fdr')]

  results <- fcHurdle %>%
    dplyr::select(primerid, `Pr(>Chisq)`, coef)

  colnames(results) <- c("gene", "p_val_adj", "avg_log2FC")
  rownames(results) <- results$gene

  return(results %>% dplyr::select("gene", "p_val_adj", "avg_log2FC"))
}

mast_ft <- function(sim, formula, contrast = "groupsGroup2"){

  fit <- mast_fit(sim, formula)
  results <- mast_test(fit, contrast)

  return(results)
}
