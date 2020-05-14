## PDR fitting fxns
## using AIC to minimize grid size
## using non-homogeneous grid size so density is higher closer to the present

fit_pdr_variable_grid <- function(tree, rho=1, starting_grid_size, age0=0,
                                  ntry_fit, ntry_search,
                                  nboot, ntry_boot, nthreads, 
                                  max_time, to_return="sim"){
  max_age <- get_tree_span(tree)$max_distance
  n_grid <- starting_grid_size
  LTT <- count_lineages_through_time(tree, Ntimes=1000)
  age_grid_i <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                   Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                   densityY=sqrt(rev(LTT$lineages)))
  age_grid_i[1] <- 0
  f_i <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_i, age0=age0, Ntrials = ntry_search,
                             min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                             max_model_runtime = max_time)
  age_grid_j <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                   Ngrid=n_grid+1, densityX=rev(max_age-LTT$times),
                                                   densityY=sqrt(rev(LTT$lineages)))
  age_grid_j[1] <- 0
  f_j <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_j, age0=age0, Ntrials = ntry_search,
                             min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                             max_model_runtime = max_time)
  aic_i <- f_i$AIC
  aic_j <- f_j$AIC
  
  ## set threshold of AIC 
  while (aic_j <= aic_i){
    aic_i <- aic_j
    n_grid <- n_grid + 1
    age_grid_j <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                     Ngrid=n_grid+1, densityX=rev(max_age-LTT$times),
                                                     densityY=sqrt(rev(LTT$lineages)))
    age_grid_j[1] <- 0
    f_j <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_j, age0=age0, Ntrials = ntry_search,
                               min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                               max_model_runtime = max_time)
    aic_j <- f_j$AIC
  }
  ## use the estimated age grid 
  age_grid <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                 Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                 densityY=sqrt(rev(LTT$lineages)))
  age_grid[1] <- 0
  
  fit_hbd_pdr_on_grid(tree, age_grid = age_grid, age0=age0, Ntrials = ntry_fit,
                             min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                             Nbootstraps = nboot, Ntrials_per_bootstrap = ntry_boot,
                             max_model_runtime = max_time)

}
  

## Misc functions for running analyses

## starting grid size function
sgs_f <- function(tree)
  1 + round(log(length(tree$tip.label)))

## function for setting the maximum time
maxt_f <- function(tree)
  max(5,length(tree$tip.label)/ 5e4)

