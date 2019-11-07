
## Modification of function for running analyses across a posterior of trees
## Reads in a ML fit
fit_pdr_dist_trees <- function(tree_dist, mle_fit, rho, ntry_fit=4, nthreads=1){
  lapply(tree_dist,function(j){
    get_pdr_single_fxn(j, mle_fit, rho, ntry_fit, nthreads)
  })
}


## little function which runs the analysis on a single tree
get_pdr_single_fxn <- function(x, mle_fit, rho, ntry_fit, nthreads){
  ## Pick the number of data points from the age grid
  n_grid <- length(mle_fit$fit$age_grid)
  ## get age grid and ltt from empirical tree
  LTT <- count_lineages_through_time(tree, Ntimes=1000)
  max_age <- get_tree_span(tree)$max_distance
  age_grid <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                 Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                 densityY=sqrt(rev(LTT$lineages)))
  age_grid[1] <- 0
  f <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid, Ntrials = ntry_fit,
                           min_PDR = -5, max_PDR=5, Nthreads = nthreads,
                           max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
  f
}

