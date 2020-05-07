## Checking robustness of conclusions to tree uncertainty

## Load in castor
library(castor)
## check version number
pv <- as.numeric(gsub(".", "", packageVersion("castor"), fixed=TRUE))
if (pv < 151)
  stop("Update package 'castor' to at least '1.5.1'")

## load in functions for fitting PDR on grid
source("R/PDR-fitting.R")

## Additional functions for computing stats
get_pulled_ext <- function(f, rho) {
  lambda_t <- f$fit$fitted_rholambda0 / rho
  rp_t <- f$fit$fitted_PDR[1]
  lambda_t - rp_t
}

## function for estimating epislon*
get_epsilon <- function(f){
  f$fitted_mu[1] / f$fitted_lambda[1]
}


## Make directory for output
ifelse(!dir.exists(file.path("output", "post")), dir.create(file.path("output", "post")), FALSE)

## Define variables for PDR fitting
ntry_fit <- 10
nthreads <- 2
sample_size <- 10
ntry_search <- 10




## Clade by clade analysis starts here

## Mammals

## Load in all trees
trees <- ape::read.nexus("data/trees-distributions/mammals/output.nex")
## Sample n trees
tmp <- sample(c(1:length(trees)), sample_size)
trees <- trees[tmp]
## Get parameters from tree files
rho <- 0.693

## Starting parameters
sgs <- sgs_f(trees[[1]])
maxt <- maxt_f(trees[[1]])

mup_0_mamm <- vector(length=sample_size)
eps_0_mamm <- vector(length=sample_size)

for (i in 1:sample_size){
  ##  Fit variable rate models (without bootstraps) 
  vrm  <- fit_pdr_variable_grid(trees[[i]], rho=rho, starting_grid_size=sgs, age0=0,
                              ntry_fit=ntry_fit, ntry_search=ntry_search,
                              nboot=NULL, nthreads=nthreads,
                              max_time = maxt)
  
  ## Fit single rate model (without bootstraps)
  max_age <- get_tree_span(trees[[i]])$max_distance
  srm <- fit_hbd_model_on_grid(trees[[i]],
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
  ## Compute summary stats
  mup_0_mamm[i] <- get_pulled_ext(vrm, rho)
  eps_0_mamm[i] <- get_epsilon(srm)
}

mamm <- data.frame(clade = rep("Mammals", sample_size), 
                   mup_0 = mup_0_mamm,
                   eps_0 = eps_0_mamm)
write.csv(mamm, "output/post/res-mammals.csv")


## Birds

## Load in all trees
trees <- ape::read.nexus("data/trees-distributions/birds/output.nex")
## Sample n trees
tmp <- sample(c(1:length(trees)), sample_size)
trees <- trees[tmp]
## Get parameters from tree files
rho <- 0.667

## Starting parameters
sgs <- sgs_f(trees[[1]])
maxt <- maxt_f(trees[[1]])

mup_0_bird <- vector(length=sample_size)
eps_0_bird <- vector(length=sample_size)

for (i in 1:sample_size){
  ##  Fit variable rate models (without bootstraps) 
  vrm  <- fit_pdr_variable_grid(trees[[i]], rho=rho, starting_grid_size=sgs, age0=0,
                                ntry_fit=ntry_fit, ntry_search=ntry_search,
                                nboot=NULL, nthreads=nthreads,
                                max_time = maxt)
  
  ## Fit single rate model (without bootstraps)
  max_age <- get_tree_span(trees[[i]])$max_distance
  srm <- fit_hbd_model_on_grid(trees[[i]],
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
  ## Compute summary stats
  mup_0_bird[i] <- get_pulled_ext(vrm, rho)
  eps_0_bird[i] <- get_epsilon(srm)
}

bird <- data.frame(clade = rep("Birds", sample_size), 
                   mup_0 = mup_0_bird,
                   eps_0 = eps_0_bird)
write.csv(bird, "output/post/res-birds.csv")


## Amphibians

## Load in all trees
trees <- ape::read.nexus("data/trees-distributions/amphibians/output.nex")
## Sample n trees
tmp <- sample(c(1:length(trees)), sample_size)
trees <- trees[tmp]
## Get parameters from tree files
rho <- 0.561

## Starting parameters
sgs <- sgs_f(trees[[1]])
maxt <- maxt_f(trees[[1]])

mup_0_amph <- vector(length=sample_size)
eps_0_amph <- vector(length=sample_size)

for (i in 1:sample_size){
  ##  Fit variable rate models (without bootstraps) 
  vrm  <- fit_pdr_variable_grid(trees[[i]], rho=rho, starting_grid_size=sgs, age0=0,
                                ntry_fit=ntry_fit, ntry_search=ntry_search,
                                nboot=NULL, nthreads=nthreads,
                                max_time = maxt)
  
  ## Fit single rate model (without bootstraps)
  max_age <- get_tree_span(trees[[i]])$max_distance
  srm <- fit_hbd_model_on_grid(trees[[i]],
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
  ## Compute summary stats
  mup_0_amph[i] <- get_pulled_ext(vrm, rho)
  eps_0_amph[i] <- get_epsilon(srm)
}

amph <- data.frame(clade = rep("Amphibians", sample_size), 
                   mup_0 = mup_0_amph,
                   eps_0 = eps_0_amph)
write.csv(amph, "output/post/res-amphibians.csv")


## Squamates

## Load in all trees
trees <- ape::read.nexus("data/trees-distributions/squamates/output.nex")
## Sample n trees
tmp <- sample(c(1:length(trees)), sample_size)
trees <- trees[tmp]
## Get parameters from tree files
rho <- 0.551

## Starting parameters
sgs <- sgs_f(trees[[1]])
maxt <- maxt_f(trees[[1]])

mup_0_squam <- vector(length=sample_size)
eps_0_squam <- vector(length=sample_size)

for (i in 1:sample_size){
  ##  Fit variable rate models (without bootstraps) 
  vrm  <- fit_pdr_variable_grid(trees[[i]], rho=rho, starting_grid_size=sgs, age0=0,
                                ntry_fit=ntry_fit, ntry_search=ntry_search,
                                nboot=NULL, nthreads=nthreads,
                                max_time = maxt)
  
  ## Fit single rate model (without bootstraps)
  max_age <- get_tree_span(trees[[i]])$max_distance
  srm <- fit_hbd_model_on_grid(trees[[i]],
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
  ## Compute summary stats
  mup_0_squam[i] <- get_pulled_ext(vrm, rho)
  eps_0_squam[i] <- get_epsilon(srm)
}

squam <- data.frame(clade = rep("Squamates", sample_size), 
                   mup_0 = mup_0_squam,
                   eps_0 = eps_0_squam)
write.csv(squam, "output/post/res-squamates.csv")



## Sharks

## Load in all trees
trees <- ape::read.nexus("data/trees-distributions/sharks/output.nex")
## Sample n trees
tmp <- sample(c(1:length(trees)), sample_size)
trees <- trees[tmp]
## Get parameters from tree files
rho <- 0.511

## Starting parameters
sgs <- sgs_f(trees[[1]])
maxt <- maxt_f(trees[[1]])

mup_0_shark <- vector(length=sample_size)
eps_0_shark <- vector(length=sample_size)

for (i in 1:sample_size){
  ##  Fit variable rate models (without bootstraps) 
  vrm  <- fit_pdr_variable_grid(trees[[i]], rho=rho, starting_grid_size=sgs, age0=0,
                                ntry_fit=ntry_fit, ntry_search=ntry_search,
                                nboot=NULL, nthreads=nthreads,
                                max_time = maxt)
  
  ## Fit single rate model (without bootstraps)
  max_age <- get_tree_span(trees[[i]])$max_distance
  srm <- fit_hbd_model_on_grid(trees[[i]],
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
  ## Compute summary stats
  mup_0_shark[i] <- get_pulled_ext(vrm, rho)
  eps_0_shark[i] <- get_epsilon(srm)
}

shark <- data.frame(clade = rep("Sharks", sample_size), 
                    mup_0 = mup_0_shark,
                    eps_0 = eps_0_shark)
write.csv(shark, "output/post/res-sharks.csv")

