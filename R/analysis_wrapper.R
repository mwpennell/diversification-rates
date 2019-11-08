source("R/PDR-fitting.R")

## read in data description
tree_info <- read.csv("data/tree_descriptions.csv")

## define ages to trim tree
trim <- c(0, 1, 2, 5)

## definte variables for fitting
ntry_fit <- 5
nthreads <- 1

 
for (i in 1:nrow(tree_info)){
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  
  ## for each age0, estimate a time-variable birth-death model and a constant model
  for (j in 1:length(trim)){
    ## fit time variable rates - calling funtion from PDR-fitting
    ## Compute bootstraps
    res  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=4, age0=trim[j],
                                  ntry_fit=ntry_fit, ntry_search=5,
                                  nboot=100, ntry_boot=1, nthreads=nthreads)
    saveRDS(res, paste0("output/pdr_boot/", trim[j], "_", clade, ".rds"))
    
    ## fit single rate models
    srm  <- fit_hbd_pdr_on_grid(tree, Ntrials = ntry_fit, age_grid=NULL, age0=trim[j],
                                min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                                max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
    saveRDS(srm, paste0("output/srm/", trim[j], "_", clade, ".rds"))
  }

}

