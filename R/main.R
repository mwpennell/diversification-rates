## Main analysis script

## Load in castor
library(castor)

## load in functions for fitting PDR on grid
source("R/PDR-fitting.R")

## Make directories for output
## for "vrm" - Variable rates model
ifelse(!dir.exists(file.path("output", "vrm")), dir.create(file.path("output", "vrm")), FALSE)
## for "srm" - Single rates model
ifelse(!dir.exists(file.path("output", "srm")), dir.create(file.path("output", "srm")), FALSE)

## Read in data description
tree_info <- read.csv("data/tree_descriptions.csv")

## Define variables for PDR fitting
ntry_fit <- 10
nthreads <- 10
nboot <- 100 
ntry_boot <- 2 
ntry_search <- 10

## Run through all trees in the dataset
for (i in 1:nrow(tree_info)){
  ## Read in trees; get metadata info from tree descriptions
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  
  ## Starting parameters
  sgs <- sgs_f(tree)
  maxt <- maxt_f(tree)
  
  ##  Fit variable rate models (with bootstraps) 
  vrm  <- fit_pdr_variable_grid(tree, rho=rho, starting_grid_size=sgs, age0=0,
                                ntry_fit=ntry_fit, ntry_search=ntry_search,
                                nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                                max_time = maxt)
  saveRDS(vrm, paste0("output/vrm/", clade, ".rds"))
  
  ## Fit constant rate birth death model
  max_age <- get_tree_span(tree)$max_distance
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
  saveRDS(srm, paste0("output/srm/", clade, ".rds"))
}
  
