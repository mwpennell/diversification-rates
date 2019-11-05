## PDR fitting fxns
## using AIC to minimize grid size
## using non-homogeneous grid size so density is higher closer to the present
library(castor)

fit_pdr_variable_grid <- function(tree, rho=1, starting_grid_size=3){
  max_age <- get_tree_span(tree)$max_distance
  n_grid <- starting_grid_size
  LTT <- count_lineages_through_time(tree, Ntimes=1000)
  age_grid_i <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                   Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                   densityY=sqrt(rev(LTT$lineages)))
  age_grid_i[1] <- 0
  f_i <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_i, Ntrials = 1,
                             min_PDR = -5, max_PDR=5, Nthreads = 1,
                             max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
  age_grid_j <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                   Ngrid=n_grid+1, densityX=rev(max_age-LTT$times),
                                                   densityY=sqrt(rev(LTT$lineages)))
  age_grid_j[1] <- 0
  f_j <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_j, Ntrials = 1,
                             min_PDR = -5, max_PDR=5, Nthreads = 1,
                             max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
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
    f_j <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_j, Ntrials = 1,
                               min_PDR = -5, max_PDR=5, Nthreads = 1,
                               max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
    aic_j <- f_j$AIC
  }
  ## use the estimated age grid 
  age_grid <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                 Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                 densityY=sqrt(rev(LTT$lineages)))
  age_grid[1] <- 0
  
  ## bootstrap the data
  f_boot <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid, Ntrials = 4,
                                min_PDR = -5, max_PDR=5, Nthreads = 2,
                                Nbootstraps = 100, Ntrials_per_bootstrap = 1,
                                max_model_runtime = max(0.5,length(tree$tip.label)/ 5e4))
  ## fill out the curve
  ntip <- length(tree$tip.label)
  pdr_mle <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age,
                                        age_grid = f_boot$age_grid, lambda0=f_boot$fitted_rholambda0/rho,
                                        mu=0, PDR=f_boot$fitted_PDR, splines_degree = 1)
  pdr_lower <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age,
                                          age_grid = f_boot$age_grid, lambda0=f_boot$CI95lower$rholambda0/rho,
                                          mu=0, PDR=f_boot$CI95lower$PDR, splines_degree = 1)
  pdr_upper <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age,
                                          age_grid = f_boot$age_grid, lambda0=f_boot$CI95upper$rholambda0/rho,
                                          mu=0, PDR=f_boot$CI95upper$PDR, splines_degree = 1)
  list(fit = f_boot, pdr_mle = pdr_mle, pdr_lower = pdr_lower, pdr_upper = pdr_upper)
}
