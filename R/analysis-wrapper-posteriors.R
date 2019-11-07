source("R/PDR-fitting-posteriors.R")

## For the following clades, run analyses across a posterior set of trees
post_trees_list <- c("Mammalia", "Squamates", "Amphibians", "Aves", "Fish", "Sharks")

## read in data description
tree_info <- read.csv("data/tree_descriptions.csv")


for (j in 1:length(post_trees_list)){
  clade <- post_trees_list[j]
  fit <- readRDS(paste0("output/pdr_boot/", clade, ".rds"))
  rho <- as.numeric(tree_info[which(tree_info$taxon == clade), "rho"])
  dist <- readRDS(paste0("data/tree_dist/", clade, ".tre"))
  pdist <- fit_pdr_dist_trees(dist, fit, rho, ntry_fit=4, nthreads=1)
  saveRDS(res, paste0("output/pdr_boot_dist/", clade, ".rds"))
}

