source("R/PDR-fitting.R")

## set up directories for putting results
## for "vrm_pdr_boot"
ifelse(!dir.exists(file.path("output", "vrm_pdr_boot")), dir.create(file.path("output", "vrm_pdr_boot")), FALSE)
## for "vrm_pdr_fit"
ifelse(!dir.exists(file.path("output", "vrm_pdr_fit")), dir.create(file.path("output", "vrm_pdr_fit")), FALSE)
## for "srm_fixed_lambda"
ifelse(!dir.exists(file.path("output", "srm_fixed_lambda")), dir.create(file.path("output", "srm_fixed_lambda")), FALSE)
## for "srm_est_lambda"
ifelse(!dir.exists(file.path("output", "srm_est_lambda")), dir.create(file.path("output", "srm_est_lambda")), FALSE)


## read in data description
tree_info <- read.csv("data/tree_descriptions.csv")

## define ages to trim tree
trim <- 2

## definte variables for fitting
ntry_fit <- 10
nthreads <- 10
nboot <- 100 
ntry_boot <- 2 
ntry_search <- 10


## starting grid size function
sgs_f <- function(tree)
  1 + round(log(length(tree$tip.label)))

## function for setting the maximum time
maxt_f <- function(tree)
  max(5,length(tree$tip.label)/ 5e4)
 
for (i in 1:nrow(tree_info)){
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  
  ## fit time variable rates - calling funtion from PDR-fitting
  sgs <- sgs_f(tree)
  maxt <- maxt_f(tree)
  vrm  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=sgs, age0=0,
                                ntry_fit=ntry_fit, ntry_search=ntry_search,
                                nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                                max_time = maxt)
  saveRDS(vrm, paste0("output/vrm_pdr_boot/", clade, ".rds"))
  
  
  ## Estimate parameters from the trimmed trees 
  ## use grid size of vrm - 1 as starting grid size
  gs <- length(vrm$fit$age_grid) - 1
    
  if (gs < 4){
    gs <- 4 ## maintain minimum grid size
  }
  vrm_trim  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=gs, age0=trim,
                                    ntry_fit=ntry_fit, ntry_search=ntry_search,
                                    nboot=NULL, ntry_boot=ntry_boot, nthreads=nthreads,
                                    max_time = maxt, to_return="fit")  
  saveRDS(vrm_trim, paste0("output/vrm_pdr_fit/", trim, "_", clade, ".rds"))
    
  ## estimate lambda0
  lambda0 <- vrm$fit$fitted_rholambda0/rho
  ## get max age
  max_age <- get_tree_span(tree)$max_distance
 
  ## fit single rate models
  ## Using estimate lambda0 from VRM
  srm <- fit_hbd_model_on_grid(tree,
                                oldest_age = max_age,
                                age0=0,
                                age_grid = 1,
                                min_mu = 0,
                                max_mu = 5,
                                guess_mu= lambda0/2,
                                fixed_lambda = lambda0,
                                fixed_rho = rho,
                                condition="crown",
                                Ntrials=ntry_fit, Nthreads=nthreads,
                                max_model_runtime = maxt)
    
  saveRDS(srm, paste0("output/srm_fixed_lambda/", clade, ".rds"))
  
  ## Estimating lambda0 and mu0 together 
  ## Using this for the final paper
  srm <- fit_hbd_model_on_grid(tree,
                               oldest_age = max_age,
                               age0=0,
                               age_grid = 1,
                               min_mu = 0,
                               max_mu = 5,
                               guess_mu= 0.5,
                               min_lambda = 0,
                               max_lambda = 5,
                               guess_lambda = 1,
                               fixed_rho = rho,
                               condition="crown",
                               Ntrials=ntry_fit, Nthreads=nthreads,
                               max_model_runtime = maxt)
  saveRDS(srm, paste0("output/srm_est_lambda/", clade, ".rds"))
  
}

## Re-estimate a CRBD model using estimated lambda rather 
## than fixing it from the VRM

for (i in 1:nrow(tree_info)){
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  max_age <- get_tree_span(tree)$max_distance

}
  
