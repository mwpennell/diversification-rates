## Figures for MS and suppmat
library(ggplot2)
library(cowplot)
library(castor)
library(ggpubr)

## Make directories if needed
## general figs folder
ifelse(!dir.exists(file.path("figs")), dir.create(file.path("figs")), FALSE)
## ms
ifelse(!dir.exists(file.path("figs", "ms")), dir.create(file.path("figs", "ms")), FALSE)
## extra
ifelse(!dir.exists(file.path("figs", "extra")), dir.create(file.path("figs", "extra")), FALSE)

## Source plotting functions
source("R/plotting-functions.R")

## Define universal variables
cols <- c("#f2706c","#1eb9f1")

mega_phylo <- c("Mammals", "Birds", "Squamates", "Amphibians",
                "Fish", "Chondrichthyans", "Agaricomycetes",
                "Vascular Plants", "Bacteria")

## Read in tree info data
tree_info <- read.csv("data/tree_descriptions.csv")


## Read in and trim output
d_vrm_boot <- "output/vrm/"

## Megaphylogenies
mega_id <- dir(d_vrm_boot)[sapply(mega_phylo, function(x) grep(x, dir(d_vrm_boot)))]
vrm_tt <- vector(mode = "list", length = length(mega_id))

for (i in 1:length(mega_id)){
  vrm_tt[[i]] <- readRDS(paste0(d_vrm_boot, mega_id[i]))
}

names(vrm_tt) <- mega_phylo


vrm_tt <- lapply(vrm_tt, function(x) trim_sims(x))

## All the other trees 
hd_id <- dir(d_vrm_boot)[-which(dir(d_vrm_boot) %in% mega_id)]
vrm_tt_supp <- vector(mode = "list", length = length(hd_id))
clade_names_supp <- vector(length=length(hd_id))

for (i in 1:length(hd_id)){
  vrm_tt_supp[[i]] <- readRDS(paste0(d_vrm_boot, hd_id[i]))
  clade_names_supp[i] <- as.character(strsplit(hd_id[i], ".", fixed=TRUE)[[1]][1])
}

names(vrm_tt_supp) <- clade_names_supp

vrm_tt_supp <- lapply(vrm_tt_supp, function(x) trim_sims(x))



## Get list of megaphylogenies for separate plotting
tree_info_mega <- tree_info[which(tree_info$taxon %in% mega_phylo), ]
tree_list_mega <- as.character(tree_info_mega$tree_name)
names(tree_list_mega) <- as.character(tree_info_mega$taxon)
## reorder
tree_list_mega <- tree_list_mega[mega_phylo]

## Get list of Henao Diaz trees for separate plotting
tree_info_hd <- tree_info[-which(tree_info$taxon %in% mega_phylo), ]
tree_list_hd <- as.character(tree_info_hd$tree_name)
names(tree_list_hd) <- as.character(tree_info_hd$taxon)
## reorder
tree_list_hd <- tree_list_hd[clade_names_supp]



## FIGURE - ltt vs dltt - macro only
figs_lm <- lapply(c(1:length(tree_list_mega)), function(x) {
  tmp <- dltt_wrapper(x, tree_list_mega, vrm_tt)
  fig_dltt_ltt(tmp$tree, dltt=tmp$dltt, clade_name=tmp$clade, cols=cols)
})

ggarrange(figs_lm[[1]], figs_lm[[2]], figs_lm[[3]],
          figs_lm[[4]], figs_lm[[5]], figs_lm[[6]],
          figs_lm[[7]], figs_lm[[8]], figs_lm[[9]],
          ncol=3, nrow=3)

ggsave("figs/ms/ltt_megaphylo.pdf", width=7.5, height=6) 




## FIGURE - ltt vs dltt - henao diaz trees
figs_lhd <- lapply(c(1:20), function(x) {
  tmp <- dltt_wrapper(x, tree_list_hd, vrm_tt_supp)
  fig_dltt_ltt_small(tmp$tree, dltt=tmp$dltt, clade_name=tmp$clade, cols=cols)
})

## Break figure up into two
ggarrange(figs_lhd[[1]], figs_lhd[[2]], figs_lhd[[3]],
          figs_lhd[[4]], figs_lhd[[5]], figs_lhd[[6]],
          figs_lhd[[7]], figs_lhd[[8]], figs_lhd[[9]],
          figs_lhd[[10]],
          figs_lhd[[11]], figs_lhd[[12]], figs_lhd[[13]],
          figs_lhd[[14]], figs_lhd[[15]], figs_lhd[[16]],
          figs_lhd[[17]], figs_lhd[[18]], figs_lhd[[19]],
          figs_lhd[[20]], 
          ncol=5, nrow=4)

ggsave("figs/ms/ltt_henaodiaz1.pdf", width=11, height=8.5) 

figs_lhdb <- lapply(c(21:39), function(x) {
  tmp <- dltt_wrapper(x, tree_list_hd, vrm_tt_supp)
  fig_dltt_ltt_small(tmp$tree, dltt=tmp$dltt, clade_name=tmp$clade, cols=cols)
})

ggarrange(figs_lhdb[[1]], figs_lhdb[[2]], figs_lhdb[[3]],
          figs_lhdb[[4]], figs_lhdb[[5]], figs_lhdb[[6]],
          figs_lhdb[[7]], figs_lhdb[[8]], figs_lhdb[[9]],
          figs_lhdb[[10]],
          figs_lhdb[[11]], figs_lhdb[[12]], figs_lhdb[[13]],
          figs_lhdb[[14]], figs_lhdb[[15]], figs_lhdb[[16]],
          figs_lhdb[[17]], figs_lhdb[[18]], figs_lhdb[[19]],
          ncol=5, nrow=4)

ggsave("figs/ms/ltt_henaodiaz2.pdf", width=11, height=8.5)


## EXTRA FIGURE - PDR through time - megaphylogenies
figs_pdrt <- lapply(1:length(vrm_tt), function(x) 
  pdr_time_boot(vrm_tt[[x]], clade_name=names(vrm_tt)[x], cols=cols))

ggarrange(figs_pdrt[[1]], figs_pdrt[[2]], figs_pdrt[[3]],
          figs_pdrt[[4]], figs_pdrt[[5]], figs_pdrt[[6]],
          figs_pdrt[[7]], figs_pdrt[[8]], figs_pdrt[[9]], 
          ncol=3, nrow=3)

ggsave("figs/extra/pdr_t_megaphylo.pdf", width=7, height=6) ## need to adjust proportions


## EXTRA FIGURE - PDR through time - Henao Diaz trees
figs_pdrt_hd <- lapply(1:20, function(x) 
  pdr_time_boot(vrm_tt_supp[[x]], clade_name=names(vrm_tt_supp)[x], cols=cols))

## Break figure up into two
ggarrange(figs_pdrt_hd[[1]], figs_pdrt_hd[[2]], figs_pdrt_hd[[3]],
          figs_pdrt_hd[[4]], figs_pdrt_hd[[5]], figs_pdrt_hd[[6]],
          figs_pdrt_hd[[7]], figs_pdrt_hd[[8]], figs_pdrt_hd[[9]],
          figs_pdrt_hd[[10]],
          figs_pdrt_hd[[11]], figs_pdrt_hd[[12]], figs_pdrt_hd[[13]],
          figs_pdrt_hd[[14]], figs_pdrt_hd[[15]], figs_pdrt_hd[[16]],
          figs_pdrt_hd[[17]], figs_pdrt_hd[[18]], figs_pdrt_hd[[19]],
          figs_pdrt_hd[[20]], 
          ncol=5, nrow=4)

ggsave("figs/extra/pdr_t_henaodiaz1.pdf", width=11, height=8.5) ## need to adjust proportions

figs_pdrt_hdb <- lapply(21:39, function(x) 
  pdr_time_boot(vrm_tt_supp[[x]], clade_name=names(vrm_tt_supp)[x], cols=cols))

ggarrange(figs_pdrt_hdb[[1]], figs_pdrt_hdb[[2]], figs_pdrt_hdb[[3]],
          figs_pdrt_hdb[[4]], figs_pdrt_hdb[[5]], figs_pdrt_hdb[[6]],
          figs_pdrt_hdb[[7]], figs_pdrt_hdb[[8]], figs_pdrt_hd[[9]],
          figs_pdrt_hdb[[10]],
          figs_pdrt_hdb[[11]], figs_pdrt_hdb[[12]], figs_pdrt_hdb[[13]],
          figs_pdrt_hdb[[14]], figs_pdrt_hdb[[15]], figs_pdrt_hdb[[16]],
          figs_pdrt_hdb[[17]], figs_pdrt_hdb[[18]], figs_pdrt_hdb[[19]],
          ncol=5, nrow=4)

ggsave("figs/extra/pdr_t_henaodiaz2.pdf", width=11, height=8.5) ## need to adjust proportions

## Remove datasets to save memory
rm(vrm_tt)
rm(vrm_tt_supp)




## FIGURE - epsilon0 vs mup0
taxa <- as.character(tree_info$taxon)
rhos <- tree_info$rho

mup_0 <- rep(NA, length=length(taxa))
eps_0 <- rep(NA, length=length(taxa))
## Estimate rp[0] just for additional info
rp_0 <- rep(NA, length=length(taxa))

for (i in 1:length(taxa)){
  file_id <- dir("output/vrm/")[grep(taxa[i], dir("output/vrm/"))]
  tmp <- readRDS(paste0("output/vrm/", file_id))
  mup_0[i] <- get_pulled_ext(tmp, rhos[i], "mle")
  rp_0[i] <- tmp$fit$fitted_PDR[1]
  file_id <- dir("output/srm/")[grep(taxa[i], dir("output/srm/"))]
  tmp <- readRDS(paste0("output/srm/", file_id))
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

mup_eps0(mu_df, cols)

ggsave("figs/ms/epsilon_mup.pdf", width=6, height=6) ## need to fix dimensions




## Calculate numbers for paper
## Number of trees
nrow(tree_info)

## Total number of taxa
sum(tree_info$ntip_tree)

## How many have negative mu_p[0] (raw)
length(which(mu_df$mup < 0)) 

## How mnay have negative mu_p[0] (percentage)
length(which(mu_df$mup < 0)) / nrow(mu_df)

## How many have epsilon approximately = 0
length(which(eps_0 < 0.01))
## How many of these have negative mup?
length(which(mu_df[which(mu_df$eps < 0.01), "mup"] < 0))
## Pvalue from binomial test
bt <- binom.test(n=length(which(eps_0 < 0.01)), 
           x=length(which(mu_df[which(mu_df$eps < 0.01), "mup"] < 0)))
bt$p.value

## How many have epsilon approximately < 0.1
length(which(eps_0 < 0.1))
## How many of these have negative mup?
length(which(mu_df[which(mu_df$eps < 0.1), "mup"] < 0))
## Pvalue from binomial test
bt <- binom.test(n=length(which(eps_0 < 0.1)), 
                 x=length(which(mu_df[which(mu_df$eps < 0.1), "mup"] < 0)))
bt$p.value

## How many have epsilon approximately < 0.25
length(which(eps_0 < 0.25))
## How many of these have negative mup?
length(which(mu_df[which(mu_df$eps < 0.25), "mup"] < 0))
## Pvalue from binomial test
bt <- binom.test(n=length(which(eps_0 < 0.25)), 
                 x=length(which(mu_df[which(mu_df$eps < 0.25), "mup"] < 0)))
bt$p.value
