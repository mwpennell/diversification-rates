source("R/PDR-fitting.R")

## set up directories for putting results
## for "vrm_pdr_boot"
ifelse(!dir.exists(file.path("output", "vrm_pdr_boot")), dir.create(file.path("output", "vrm_pdr_boot")), FALSE)
## for "vrm_pdr_fit"
ifelse(!dir.exists(file.path("output", "vrm_pdr_fit")), dir.create(file.path("output", "vrm_pdr_fit")), FALSE)
## for "vrm_fixed_lambda"
ifelse(!dir.exists(file.path("output", "vrm_fixed_lambda")), dir.create(file.path("output", "vrm_fixed_lambda")), FALSE)
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
  
  ## for each age0, estimate a time-variable birth-death models and a constant model
  for (j in 1:length(trim)){
    ## fit time variable rates - calling funtion from PDR-fitting
    ## Compute bootstraps if age0 = 0
    if (trim[j] == 0){
      vrm  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=sgs, age0=trim[j],
                                    ntry_fit=ntry_fit, ntry_search=ntry_search,
                                    nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads)
      saveRDS(vrm, paste0("output/vrm_pdr_boot/", clade, ".rds"))
    } else { ## if age0 != 0 only fit the model
      ## don't impute the rest of the time period (to save space)
      vrm  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=sgs, age0=trim[j],
                                    ntry_fit=ntry_fit, ntry_search=ntry_search,
                                    nboot=NULL, ntry_boot=ntry_boot, nthreads=nthreads,
                                    to_return="fit")  
      saveRDS(vrm, paste0("output/vrm_pdr_fit/", trim[j], "_", clade, ".rds"))
    }
    
    ## fit variable extinciton model with fixed lambda value
    ## use same grid as PDR model
    ## estimate lambda0
    lambda0 <- vrm$fit$fitted_rholambda0/rho
    ## get max age
    max_age <- get_tree_span(tree)$max_distance
    ## age_grid
    age_grid <- vrm$fit$age_grid
    
    vmu <- fit_hbd_model_on_grid(tree,
                                 oldest_age = max_age,
                                 age0=trim[j],
                                 age_grid = age_grid,
                                 min_mu = 0,
                                 max_mu = 5,
                                 guess_mu= lambda0/2,
                                 fixed_lambda = lambda0,
                                 fixed_rho = rho,
                                 condition="crown",
                                 Ntrials=ntry_fit, Nthreads=nthreads,
                                 max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
    saveRDS(srm, paste0("output/vrm_fixed_lambda/", trim[j], "_", clade, ".rds"))
    
    ## fit single rate models
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
                                  max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
    
    saveRDS(srm, paste0("output/srm/", trim[j], "_", clade, ".rds"))
  }

}

