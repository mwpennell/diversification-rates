## Plotting functions
## Called by plots.R


## Function for thinning data down to alleviate strain on ggplot

trim_sims <- function(x, out=10000){
  samp_mle <- round(seq(1, length(x$pdr_mle$PDR), length.out = out))
  x$pdr_mle$PDR <- x$pdr_mle$PDR[samp_mle]
  x$pdr_mle$LTT <- x$pdr_mle$LTT[samp_mle]
  x$pdr_mle$ages <- x$pdr_mle$ages[samp_mle]
  samp_up <- seq(1, length(x$pdr_upper$PDR), length.out = out)
  x$pdr_upper$PDR <- x$pdr_upper$PDR[samp_up]
  x$pdr_upper$LTT <- x$pdr_upper$LTT[samp_up]
  x$pdr_upper$ages <- x$pdr_upper$ages[samp_up]
  samp_lo <- seq(1, length(x$pdr_lower$PDR), length.out = out)
  x$pdr_lower$PDR <- x$pdr_lower$PDR[samp_lo]
  x$pdr_lower$LTT <- x$pdr_lower$LTT[samp_lo]
  x$pdr_lower$ages <- x$pdr_lower$ages[samp_lo]
  x
}




## Function for dltt and ltt comparison

fig_dltt_ltt <- function(tree, dltt, clade_name, cols){
  tree_age <- get_tree_span(tree)$max_distance
  ## get the ltt at every point
  
  ltt <- castor::count_lineages_through_time(tree, Ntimes=1000)
  lineages <- ltt$lineages[1:(length(ltt$lineages) - 1)]
  lineage_times <- ltt$times[1:(length(ltt$times) - 1)]
  
  ltt_comp <- data.frame(lineages=c(lineages, rev(dltt$LTT)),
                         ages=c(lineage_times, dltt$ages),
                         Type=c(rep("LTT", length(lineages)),
                                rep("dLTT", length(dltt$LTT))))
  ltt_comp$ages <- max(ltt_comp$ages) - ltt_comp$ages
  
  
  ggplot(ltt_comp, aes(x=ages, y=lineages, colour=Type)) + 
    geom_line() + scale_y_log10() + theme_cowplot() +
    scale_colour_manual(values = cols) + labs(title = clade_name) +
    xlab("Time before present (My)") + ylab("Lineages") +
    theme(axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
          axis.text = element_text(size=7), plot.title = element_text(size=10),
          legend.position = "none")
}


## Slightly different text sizes
fig_dltt_ltt_small <- function(tree, dltt, clade_name, cols){
  tree_age <- get_tree_span(tree)$max_distance
  ## get the ltt at every point
  
  ltt <- castor::count_lineages_through_time(tree, Ntimes=1000)
  lineages <- ltt$lineages[1:(length(ltt$lineages) - 1)]
  lineage_times <- ltt$times[1:(length(ltt$times) - 1)]
  
  ltt_comp <- data.frame(lineages=c(lineages, rev(dltt$LTT)),
                         ages=c(lineage_times, dltt$ages),
                         Type=c(rep("LTT", length(lineages)),
                                rep("dLTT", length(dltt$LTT))))
  ltt_comp$ages <- max(ltt_comp$ages) - ltt_comp$ages
  
  
  ggplot(ltt_comp, aes(x=ages, y=lineages, colour=Type)) + 
    geom_line() + scale_y_log10() + theme_cowplot() +
    scale_colour_manual(values = cols) + labs(title = clade_name) +
    xlab("Time before present (My)") + ylab("Lineages") +
    theme(axis.title.x = element_text(size=5), axis.title.y = element_text(size=5),
          axis.text = element_text(size=5), plot.title = element_text(size=8),
          legend.position = "none")
}



## Wrapper function for iterating the above
dltt_wrapper <- function(i, tree_list, dltt_out){
  t <- read_tree(file=paste0("data/trees/", 
                             as.character(tree_list[i]), ".tre"))
  clade <- as.character(names(tree_list)[i])
  f <- dltt_out[[clade]]$pdr_mle
  list(tree=t, dltt=f, clade=clade)
}



## function for extracting pulled extinction rates at time t 
get_pulled_ext <- function(f, rho, type="mle") {
  if (type == "mle"){
    lambda_t <- f$fit$fitted_rholambda0 / rho
    rp_t <- f$fit$fitted_PDR[1]
    mu_p <- lambda_t - rp_t
  } else {
    stop("only 'mle' is implemented at this point")
  }
  mu_p
}

## function for estimating epislon*
get_epsilon <- function(f){
  f$fitted_mu[1] / f$fitted_lambda[1]
}

## Function for making plot comparing mup to epsilon0
mup_eps0 <- function(mu_df, cols){
  ggplot(mu_df, aes(x=mup, y=eps)) + 
    geom_point(colour=cols[1], size=3, alpha=0.9) +
    theme_cowplot() + xlab(expression(mu[p](0))) + 
    ylab(expression(epsilon[o]~"*")) + 
    geom_vline(xintercept=0, colour=cols[2],linetype='dashed') + 
    geom_hline(yintercept=0, colour=cols[2],linetype='dashed') +
    geom_hline(yintercept=1, colour=cols[2],linetype='dashed') +
    theme(axis.text = element_text(size=12))
}




## Function for visualizing PDR through time
pdr_time_boot <- function(res, clade_name, cols){
  df <- data.frame(PDR=c(res$pdr_lower$PDR, res$pdr_upper$PDR, res$pdr_mle$PDR), 
                   type=c(rep("lower", length(res$pdr_lower$PDR)), rep("upper", length(res$pdr_upper$PDR)), rep("mle", length(res$pdr_mle$PDR))),
                   ages=c(res$pdr_lower$ages, res$pdr_upper$ages, res$pdr_mle$ages))
  ggplot(df, aes(x=ages, y=PDR, color=type)) + geom_line() + #stat_smooth(se=FALSE) + 
    scale_colour_manual(values = c(cols[1], cols[2], cols[1])) + 
    theme_cowplot() + theme(legend.position = "none") + 
    xlab("Time before present (My)") + ylab(expression(r[p])) +
    labs(title = clade_name) +
    theme(axis.title.x = element_text(size=8), axis.title.y = element_text(size=8),
          axis.text = element_text(size=6), plot.title = element_text(size=10))
}



## Function for examining variation in mu_p0 across posterior
mu0p_post <- function(df){
  ggplot(df, aes(x=clade, y=mup_0)) + geom_point()
}
