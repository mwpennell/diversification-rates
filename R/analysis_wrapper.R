source("R/PDR-fitting.R")

## set up directories for putting results
## for "pdr_boot"
ifelse(!dir.exists(file.path("output", "pdr_boot")), dir.create(file.path("output", "pdr_boot")), FALSE)
## for "srm"
ifelse(!dir.exists(file.path("output", "srm")), dir.create(file.path("output", "srm")), FALSE)


## read in data description
tree_info <- read.csv("data/tree_descriptions.csv")

## define ages to trim tree
trim <- c(0, 2)

## definte variables for fitting
ntry_fit <- 5
nthreads <- 1
nboot <- 100
ntry_boot <- 1
ntry_search <- 5
sgs <- 4 ## starting grid size
 
for (i in 1:nrow(tree_info)){
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  
  ## for each age0, estimate a time-variable birth-death model and a constant model
  for (j in 1:length(trim)){
    ## fit time variable rates - calling funtion from PDR-fitting
    ## Compute bootstraps
    vrm  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=sgs, age0=trim[j],
                                  ntry_fit=ntry_fit, ntry_search=ntry_search,
                                  nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads)
    saveRDS(vrm, paste0("output/pdr_boot/", trim[j], "_", clade, ".rds"))
    
    ## fit single rate models
    ## estimate lambda0
    lambda0 <- vrm$fit$fitted_rholambda0/rho
    max_age <- get_tree_span(tree)$max_distance
    srm <- fit_hbd_model_on_grid(tree,
                                  oldest_age = max_age,
                                  age0=trim[j],
                                  age_grid = 1,
                                  min_mu = 0,
                                  max_mu = 5,
                                  guess_mu= lambda0/2,
                                  fixed_lambda = lambda0,
                                  fixed_rho = rho,
                                  condition="crown",
                                  Ntrials=ntry_fit, Nthreads=nthreads,
                                  max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4),
                                  control=list(eval.max=500, iter.max=200, rel.tol=1e-6))
    
    saveRDS(srm, paste0("output/srm/", trim[j], "_", clade, ".rds"))
  }

}

