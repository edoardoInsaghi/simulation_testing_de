rm(list = ls())
require(splatter)

simulate <- function(
    id = 1,
    sparsity.cut = 0.95,
    nGenes = 5000,
    batchCells = c(500,500),
    group.prob = c(0.5, 0.5),
    de.prob = c(0.1, 0.1),
    de.downProb = c(0.4, 0.4),
    de.facLoc = 0.1,
    de.facScale = 0.1,
    batch.facLoc = 0.4,
    batch.facScale = 0.4,
    filter.data = TRUE){

  set.seed(id)

  sim <- splatter::splatSimulate(method="groups",
                       nGenes = nGenes,
                       batchCells = batchCells,
                       group.prob = group.prob,
                       de.prob = de.prob,
                       de.downProb = de.downProb,
                       de.facLoc = de.facLoc,
                       de.facScale = de.facScale,
                       batch.facLoc = batch.facLoc,
                       batch.facScale = batch.facScale,
                       lib.loc = 10.5,
                       lib.scale = 0.5)

  counts.data <- as.data.frame(counts(sim))
  gene.data <- as.data.frame(rowData(sim))
  cell.data <- as.data.frame(colData(sim))

  # Removing Sparse Genes
  if (filter.data){

    # counting the number of zero occurrences of the genes in each cell, then
    # dividing by the total number of cells.
    row.mean <- rowSums(counts(sim) == 0) / ncol(counts(sim))
    non_sparse_genes <- rownames(counts(sim))[which(row.mean < sparsity.cut)]


    counts.data <- counts.data[non_sparse_genes, ]
    gene.data <- gene.data[non_sparse_genes, ]
  }

  # List of Differentially Expressed Genes
  de.genes <- rownames(counts.data)[which(gene.data$DEFacGroup1 != gene.data$DEFacGroup2)]


  return( list(simulation = sim,
               counts.data = counts.data,
               gene.data = gene.data,
               cell.data = cell.data,
               de.genes = de.genes) )
}


n_genes <- list(5000, 10000, 20000)
batch_cells <- list(rep(500, 2), rep(5000, 2), rep(10000, 2))
n_repetitions <- 50
path_sim_folder <- "simulations/"

id <- 1
for (i in 1:length(n_genes)) {
  for (j in 1:length(batch_cells)) {
    for (k in 1:n_repetitions) {
      sim <- simulate(id=id, nGenes = n_genes[[i]], batchCells = batch_cells[[j]])
      name <- paste0(path_sim_folder, "sim_", "ng_", n_genes[[i]], "_nc_", sum(batch_cells[[j]]), "_", id, ".rds")
      saveRDS(sim, name)
      id <- id + 1
    }
  }
}

