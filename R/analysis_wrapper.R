source("R/PDR-fitting.R")

## read in data description
tree_info <- read.csv("data/tree_descriptions.csv")

for (i in 1:nrow(tree_info)){
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  res  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=4, 
                                ntry_fit=5, ntry_search=5,
                                nboot=100, ntry_boot=2, nthreads=1)
  saveRDS(res, paste0("output/pdr_boot/", clade, ".rds"))
}

