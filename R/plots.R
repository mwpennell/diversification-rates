## Figures for MS and suppmat
library(ggplot2)
library(cowplot)
library(castor)
library(ggpubr)

## Functions for making plots

## Define universal variables
cols <- c("#5475BA", "#2E3C61")



## Figure 1: PDR through time -- megaphylogenies only

mega_phylo <- c("Mammals", "Birds", "Squamates", "Amphibians",
                "Fish", "Chondrichthyans", "Agaricomycetes",
                "Vascular Plants", "Bacteria")

## read in vrm boot trees
d_vrm_boot <- "output/vrm_pdr_boot/"
mega_id <- dir(d_vrm_boot)[sapply(mega_phylo, function(x) grep(x, dir(d_vrm_boot)))]
vrm_tt <- vector(mode = "list", length = length(mega_id))

for (i in 1:length(mega_id)){
  vrm_tt[[i]] <- readRDS(paste0(d_vrm_boot, mega_id[i]))
}

names(vrm_tt) <- mega_phylo

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

vrm_tt <- lapply(vrm_tt, function(x) trim_sims(x))


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
          axis.text = element_text(size=6), plot.title = element_text(size=8))
}


figs_1 <- lapply(1:length(vrm_tt), function(x) 
  pdr_time_boot(vrm_tt[[x]], clade_name=names(vrm_tt)[x], cols=cols))

ggarrange(figs_1[[1]], figs_1[[2]], figs_1[[3]],
          figs_1[[4]], figs_1[[5]], figs_1[[6]],
          figs_1[[7]], figs_1[[8]], figs_1[[9]], 
          ncol=3, nrow=3)

ggsave("figs/ms/pdr_time_megaphylo.pdf", width=7, height=6) ## need to adjust proportions








## Figure 2: PDR through time -- everything else (Henao Diaz trees)
## cut out Seed Plants for the time being
mega_id_sp <- c(mega_id, "Seed Plants.rds")
hd_id <- dir(d_vrm_boot)[-which(dir(d_vrm_boot) %in% mega_id_sp)]
vrm_tt_supp <- vector(mode = "list", length = length(hd_id))
clade_names_supp <- vector(length=length(hd_id))

for (i in 1:length(hd_id)){
  vrm_tt_supp[[i]] <- readRDS(paste0(d_vrm_boot, hd_id[i]))
  clade_names_supp[i] <- as.character(strsplit(hd_id[i], ".", fixed=TRUE)[[1]][1])
}

names(vrm_tt_supp) <- clade_names_supp

vrm_tt_supp <- lapply(vrm_tt_supp, function(x) trim_sims(x))

figs_s1 <- lapply(1:20, function(x) 
  pdr_time_boot(vrm_tt_supp[[x]], clade_name=names(vrm_tt_supp)[x], cols=cols))

## Break figure up into two
ggarrange(figs_s1[[1]], figs_s1[[2]], figs_s1[[3]],
          figs_s1[[4]], figs_s1[[5]], figs_s1[[6]],
          figs_s1[[7]], figs_s1[[8]], figs_s1[[9]],
          figs_s1[[10]],
          figs_s1[[11]], figs_s1[[12]], figs_s1[[13]],
          figs_s1[[14]], figs_s1[[15]], figs_s1[[16]],
          figs_s1[[17]], figs_s1[[18]], figs_s1[[19]],
          figs_s1[[20]], 
          ncol=5, nrow=4)

ggsave("figs/ms/pdr_time_henaodiaz1.pdf", width=11, height=8.5) ## need to adjust proportions

figs_s1b <- lapply(21:39, function(x) 
  pdr_time_boot(vrm_tt_supp[[x]], clade_name=names(vrm_tt_supp)[x], cols=cols))

ggarrange(figs_s1b[[1]], figs_s1b[[2]], figs_s1b[[3]],
         figs_s1b[[4]], figs_s1b[[5]], figs_s1b[[6]],
         figs_s1b[[7]], figs_s1b[[8]], figs_s1[[9]],
         figs_s1b[[10]],
         figs_s1b[[11]], figs_s1b[[12]], figs_s1b[[13]],
         figs_s1b[[14]], figs_s1b[[15]], figs_s1b[[16]],
         figs_s1b[[17]], figs_s1b[[18]], figs_s1b[[19]],
         ncol=5, nrow=4)

ggsave("figs/ms/pdr_time_henaodiaz2.pdf", width=11, height=8.5) ## need to adjust proportions






## Figure 3 and 4: Comparison of LTT and dLTT plots
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
    theme(axis.title.x = element_text(size=7), axis.title.y = element_text(size=7),
          axis.text = element_text(size=6), plot.title = element_text(size=8),
          legend.position = "none")
}

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
    xlab("Mya") + ylab("Lineages") +
    theme(axis.title.x = element_text(size=5), axis.title.y = element_text(size=5),
          axis.text = element_text(size=5), plot.title = element_text(size=6),
          legend.position = "none")
}


tree_info <- read.csv("data/tree_descriptions.csv")



## 3. Megaphylogenies only

tree_info_mega <- tree_info[which(tree_info$taxon %in% mega_phylo), ]
tree_list_mega <- as.character(tree_info_mega$tree_name)
names(tree_list_mega) <- as.character(tree_info_mega$taxon)
## reorder
tree_list_mega <- tree_list_mega[mega_phylo]


dltt_wrapper <- function(i, tree_list, dltt_out){
  t <- read_tree(file=paste0("data/trees/", 
                             as.character(tree_list[i]), ".tre"))
  clade <- as.character(names(tree_list)[i])
  f <- dltt_out[[clade]]$pdr_mle
  list(tree=t, dltt=f, clade=clade)
}


figs_3 <- lapply(c(1:length(tree_list_mega)), function(x) {
  tmp <- dltt_wrapper(x, tree_list_mega, vrm_tt)
  fig_dltt_ltt(tmp$tree, dltt=tmp$dltt, clade_name=tmp$clade, cols=cols)
})

ggarrange(figs_3[[1]], figs_3[[2]], figs_3[[3]],
          figs_3[[4]], figs_3[[5]], figs_3[[6]],
          figs_3[[7]], figs_3[[8]], figs_3[[9]],
          ncol=3, nrow=3)

ggsave("figs/ms/ltt_megaphylo.pdf", width=8, height=6) ## need to adjust proportions



## Figure 4: Comparison of LTT and dLTT plots
## Henao Diaz trees only
## remove: Seedplants 
mega_phylo_sp <- c(mega_phylo, "Seed Plants")
tree_info_hd <- tree_info[-which(tree_info$taxon %in% mega_phylo_sp), ]
tree_list_hd <- as.character(tree_info_hd$tree_name)
names(tree_list_hd) <- as.character(tree_info_hd$taxon)
## reorder
tree_list_hd <- tree_list_hd[clade_names_supp]

figs_s2 <- lapply(c(1:20), function(x) {
  tmp <- dltt_wrapper(x, tree_list_hd, vrm_tt_supp)
  fig_dltt_ltt_small(tmp$tree, dltt=tmp$dltt, clade_name=tmp$clade, cols=cols)
})

## Break figure up into two
ggarrange(figs_s2[[1]], figs_s2[[2]], figs_s2[[3]],
          figs_s2[[4]], figs_s2[[5]], figs_s2[[6]],
          figs_s2[[7]], figs_s2[[8]], figs_s2[[9]],
          figs_s2[[10]],
          figs_s2[[11]], figs_s2[[12]], figs_s2[[13]],
          figs_s2[[14]], figs_s2[[15]], figs_s2[[16]],
          figs_s2[[17]], figs_s2[[18]], figs_s2[[19]],
          figs_s2[[20]], 
          ncol=5, nrow=4)

ggsave("figs/ms/ltt_henaodiaz1.pdf", width=11, height=8.5) ## need to adjust proportions

figs_s2b <- lapply(c(21:39), function(x) {
  tmp <- dltt_wrapper(x, tree_list_hd, vrm_tt_supp)
  fig_dltt_ltt_small(tmp$tree, dltt=tmp$dltt, clade_name=tmp$clade, cols=cols)
})

ggarrange(figs_s2b[[1]], figs_s2b[[2]], figs_s2b[[3]],
         figs_s2b[[4]], figs_s2b[[5]], figs_s2b[[6]],
         figs_s2b[[7]], figs_s2b[[8]], figs_s2b[[9]],
         figs_s2b[[10]],
         figs_s2b[[11]], figs_s2b[[12]], figs_s2b[[13]],
         figs_s2b[[14]], figs_s2b[[15]], figs_s2b[[16]],
         figs_s2b[[17]], figs_s2b[[18]], figs_s2b[[19]],
         ncol=5, nrow=4)

ggsave("figs/ms/ltt_henaodiaz2.pdf", width=11, height=8.5) ## need to adjust proportions







## Figure 5: mu_p(x) vs. epsilon* 


## function for extracting pulled extinction rates at time t 
get_pulled_ext <- function(f, rho) {
  lambda_t <- f$fit$fitted_rholambda0 / rho
  rp_t <- f$fit$fitted_PDR[1]
  lambda_t - rp_t
}

## function for estimating epislon*
get_epsilon <- function(f){
  f$fitted_mu[1] / f$fitted_lambda[1]
}

taxa <- as.character(tree_info$taxon)
## tmp: remove Seed Plants
sp <- grep("Seed Plants", taxa)
taxa <- taxa[-sp]
rhos <- tree_info$rho[-sp]

mup_0 <- rep(NA, length=length(taxa))
eps_0 <- rep(NA, length=length(taxa))
## Estimate rp[0] just for additional info
rp_0 <- rep(NA, length=length(taxa))

for (i in 1:length(taxa)){
  file_id <- dir("output/vrm_pdr_boot/")[grep(taxa[i], dir("output/vrm_pdr_boot/"))]
  tmp <- readRDS(paste0("output/vrm_pdr_boot/", file_id))
  mup_0[i] <- get_pulled_ext(tmp, rhos[i])
  rp_0[i] <- tmp$fit$fitted_PDR[1]
  file_id <- dir("output/srm_est_lambda/")[grep(taxa[i], dir("output/srm_fixed_lambda/"))]
  tmp <- readRDS(paste0("output/srm_est_lambda/", file_id))
  if (tmp$success == TRUE){
    eps_0[i] <- get_epsilon(tmp)
  }
}

mu_df <- data.frame(taxa= taxa, 
                    mup=mup_0,
                    eps=eps_0,
                    rp_0=rp_0)

## some summary stats for the paper
## how many mu_p are negative

## for age0
## absolute numbers
length(which(mu_df$mup < 0)) 
## as a percentage
length(which(mu_df$mup < 0)) / nrow(mu_df)

## summary stats of rp[0]
length(which(mu_df$rp_0 > 0))
median(mu_df$rp_0)
sd(mu_df$rp_0)

ggplot(mu_df, aes(x=mup, y=eps)) + geom_point(colour=cols[1], size=4) + 
  theme_cowplot() + xlab(expression(mu[p](0))) + 
  ylab(expression(epsilon)) + 
  geom_vline(xintercept=0, colour=cols[2],linetype='dashed') + 
  geom_hline(yintercept=0, colour=cols[2],linetype='dashed') +
  geom_hline(yintercept=1, colour=cols[2],linetype='dashed') +
  theme(axis.text = element_text(size=12))

ggsave("figs/ms/epsilon_est_mup_age0.pdf") ## need to fix dimensions










