## modified version of main script
## for getting confidence intervals after the fact

## Main analysis script

## Load in castor
library(castor)

## load in functions for fitting PDR on grid
source("R/PDR-fitting.R")

## Make directories for output
## for "vrm_ci" - Variable rates model
ifelse(!dir.exists(file.path("output", "vrm_ci")), dir.create(file.path("output", "vrm_ci")), FALSE)

## Read in data description
tree_info <- read.csv("data/tree_descriptions.csv")


## Collect data on grid size
## all_clades <- vector(length=nrow(tree_info))
## grid_size  <- vector(length=nrow(tree_info))
## for (i in 1:nrow(tree_info)){
##  cl <- as.character(tree_info[i, "taxon"])
##  all_clades[i] <- cl
##  tmp <- readRDS(paste0("output/vrm/", cl, ".rds"))
##  grid_size[i] <- length(tmp$fit$age_grid)
## }
## gs <- data.frame(clade = all_clades, grid_points = grid_size)
## write.csv(gs, "output/grid_sizes.csv")

## Read in data on grid size used for fitting
grid_size <- read.csv("output/grid_sizes.csv")

## Define variables for PDR fitting
ntry_fit <- 10
nthreads <- parallel::detectCores() - 2
nboot <- 100 
ntry_boot <- 2 
## ntry_search <- 10

## Run through all trees in the dataset
for (i in 1:nrow(tree_info)){
  ## Read in trees; get metadata info from tree descriptions
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  
  ## Starting parameters
  sgs <- as.numeric(grid_size[which(grid_size$clade == clade), "grid_points"])
  maxt <- maxt_f(tree)
  max_age <- get_tree_span(tree)$max_distance
  LTT <- count_lineages_through_time(tree, Ntimes=1000)
  age_grid <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                 Ngrid=sgs, densityX=rev(max_age-LTT$times),
                                                 densityY=sqrt(rev(LTT$lineages)))
  age_grid[1] <- 0
  vrm <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid, age0=0, Ntrials = ntry_fit,
                                min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                                Nbootstraps = nboot, Ntrials_per_bootstrap = ntry_boot,
                                max_model_runtime = maxt)
  saveRDS(vrm, paste0("output/vrm_ci/", clade, ".rds")) }

