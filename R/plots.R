## Figures for MS and suppmat
library(ggplot2)
library(cowplot)
library(castor)

## Functions for making plots

## Define universal variables
cols <- c("#5475BA", "#2E3C61")

d <- "output/pdr_boot/"
pdr_data <- vector(mode = "list", length = length(dir(d)))
for (i in 1:length(dir(d))){
  pdr_data[[i]] <- readRDS(paste0(d, dir(d)[i]))
}
clade_names <- as.character(sapply(dir(d), function(x) 
  return(strsplit(x, ".", fixed=TRUE)[[1]][1])))
names(pdr_data) <- clade_names


## Plot 1 : PDR through time (combine later)
pdr_time <- function(res, clade_name, cols=c("#5475BA", "#2E3C61")){
  df <- data.frame(PDR=c(res$pdr_lower$PDR, res$pdr_upper$PDR, res$pdr_mle$PDR), 
                   type=c(rep("lower", length(res$pdr_lower$PDR)), rep("upper", length(res$pdr_upper$PDR)), rep("mle", length(res$pdr_mle$PDR))),
                   ages=c(res$pdr_lower$ages, res$pdr_upper$ages, res$pdr_mle$ages))
  ggplot(df, aes(x=ages, y=PDR, color=type)) + stat_smooth(se=FALSE) + 
    scale_colour_manual(values = c("#5475BA", "#2E3C61", "#5475BA")) + 
    theme_cowplot() + theme(legend.position = "none") + 
    xlab("Time before present (My)") + ylab(expression(r[p]))
  ggsave(paste0("figs/pdr_time/", clade_name, ".pdf"))
}

lapply(1:length(pdr_data), function(x) 
  pdr_time(pdr_data[[x]], clade_name = names(pdr_data)[x]))


## Plot 2: PSR through time
psr_time <- function(res, clade_name, cols=c("#5475BA", "#2E3C61")){
  df <- data.frame(PSR=c(res$pdr_lower$PSR, res$pdr_upper$PSR, res$pdr_mle$PSR), 
                   type=c(rep("lower", length(res$pdr_lower$PSR)), rep("upper", length(res$pdr_upper$PSR)), rep("mle", length(res$pdr_mle$PSR))),
                   ages=c(res$pdr_lower$ages, res$pdr_upper$ages, res$pdr_mle$ages))
  ggplot(df, aes(x=ages, y=PSR, color=type)) + stat_smooth(se=FALSE) + 
    scale_colour_manual(values = c("#5475BA", "#2E3C61", "#5475BA")) + 
    theme_cowplot() + theme(legend.position = "none") + 
    xlab("Time before present (My)") + ylab(expression(lambda[p]))
  ggsave(paste0("figs/psr_time/", clade_name, ".pdf"))
}

lapply(1:length(pdr_data), function(x) 
  psr_time(pdr_data[[x]], clade_name = names(pdr_data)[x]))


## Plot 3: Lambda0 vs. clade age

## grab rho from each dataset from data table
tree_info <- read.csv("data/tree_descriptions.csv")
rhos <- tree_info$rho
names(rhos) <- tree_info$taxon
rhos <- rhos[clade_names]
rhos <- rhos[!is.na(rhos)]

## build a dataframe with all the information
lambda0ml <- sapply(pdr_data, function(x) return(x$fit$fitted_rholambda0))
lambda0lo <- sapply(pdr_data, function(x) return(x$fit$CI95lower$rholambda0))
lambda0up <- sapply(pdr_data, function(x) return(x$fit$CI95upper$rholambda0))
lambda0se <- sapply(pdr_data, function(x) return(x$fit$standard_errors$rholambda0))
rp0ml <- sapply(pdr_data, function(x) return(x$fit$fitted_PDR[1]))
rp0lo <- sapply(pdr_data, function(x) return(x$fit$CI95lower$PDR[1]))
rp0up <- sapply(pdr_data, function(x) return(x$fit$CI95upper$PDR[1]))
rp0se <- sapply(pdr_data, function(x) return(x$fit$standard_errors$PDR[1]))
age <- sapply(pdr_data, function(x) return(max(x$fit$age_grid)))

pdr_df <- data.frame(lambda0rho = lambda0ml, lambda0rho_lo = lambda0lo,
                     lambda0rho_up = lambda0up, lambda0rho_se = lambda0se,
                     rp0 = rp0ml, rp0lo = rp0lo, rp0up = rp0up, rp0se = rp0se, 
                     rho = rhos, age=age, clade=clade_names)
## solve for lambda0 given rho
pdr_df$lambda0 <- pdr_df$lambda0rho / pdr_df$rho
## calculate mu_p0 given lambda0 and rho
pdr_df$mu_p0 <- pdr_df$lambda0 - pdr_df$rp0

## Exclude bacteria since it is so much older
pdr_df_nb <- pdr_df[-which(pdr_df$clade == "Bacteria"), ]

ggplot(pdr_df_nb, aes(x=age, y=lambda0)) +  
  geom_errorbar(aes(ymin=lambda0-lambda0se, ymax=lambda0+lambda0se)) +
  geom_point(colour="#5475BA", size=3, alpha=0.7) + 
  ylab(expression(lambda(0))) + xlab("Age of clade (my)") + 
  theme_cowplot() #+ stat_smooth(se=FALSE, colour="#2E3C61")
ggsave("figs/clade_age/lambda0_cladeage.pdf")


## Plot 4: rp0 vs clade age (standard scale)
## Exclude bacteria since it is so much older than the rest

ggplot(pdr_df_nb, aes(x=age, y=rp0)) +  
  geom_errorbar(aes(ymin=rp0-rp0se, ymax=rp0+rp0se)) +
  geom_point(colour="#5475BA", size=3, alpha=0.7) + 
  ylab(expression(r[p](0))) + xlab("Age of clade (my)") + 
  theme_cowplot() + stat_smooth(se=FALSE, colour="#2E3C61")
ggsave("figs/clade_age/rp0_cladeage.pdf")


## Plot 6: mu_p0 vs. epislon from constant rate model
tmp <- grep(pattern="0_", x=dir("output/srm"))
ds <- "output/srm/"
srm_data <- vector(mode = "list", length = length(dir(ds)[tmp]))
for (i in 1:length(dir(ds)[tmp])){
  srm_data[[i]] <- readRDS(paste0(ds, dir(ds)[i]))
}
clade_names_srm <- as.character(sapply(dir(ds)[tmp], function(x) 
  return(strsplit(x, "_", fixed=TRUE)[[1]][2])))
clade_names_srm <- as.character(sapply(clade_names_srm, function(x) 
  return(strsplit(x, ".", fixed=TRUE)[[1]][1])))
names(srm_data) <- clade_names_srm

## Compute epsilon
lambda0ml_srm <- sapply(srm_data, function(x) return(x$fit$fitted_rholambda0))
#lambda0lo_srm <- sapply(srm_data, function(x) return(x$fit$CI95lower$rholambda0))
#lambda0up_srm <- sapply(srm_data, function(x) return(x$fit$CI95upper$rholambda0))
#lambda0se_srm <- sapply(srm_data, function(x) return(x$fit$standard_errors$rholambda0))
rp0ml_srm <- sapply(srm_data, function(x) return(x$fit$fitted_PDR[1]))
#rp0lo_srm <- sapply(srm_data, function(x) return(x$fit$CI95lower$PDR[1]))
#rp0up_srm <- sapply(srm_data, function(x) return(x$fit$CI95upper$PDR[1]))
#rp0se_srm <- sapply(srm_data, function(x) return(x$fit$standard_errors$PDR[1]))
#age_srm <- sapply(srm_data, function(x) return(max(x$fit$age_grid)))

srm_df <- data.frame(lambda0rho_srm = lambda0ml_srm, rp0_srm = rp0ml_srm, 
                     rho = rhos[clade_names_srm], clade=clade_names_srm)
## compute lambda0
srm_df$lambda0_srm <- srm_df$lambda0rho_srm / srm_df$rho
## compute epsilon
srm_df$mu0 <- srm_df$lambda0_srm - srm_df$rp0_srm
srm_df$eps <- srm_df$mu / srm_df$lambda0_srm

## Plot 7: LTT vs. dLTT plots
ltt_dltt <- function(tree, fit, cols=c("#5475BA", "#2E3C61")){
  tree_age <- get_tree_span(tree)$max_distance
  ltt <- castor::count_lineages_through_time(tree, times=seq(0, (tree_age - 0.1), length.out = 100))
  dltt <- fit$pdr_mle
  min_one_lineage <- which(dltt$LTT >= 1)
  dlineages <- dltt$LTT[min_one_lineage]
  dages <- tree_age - dltt$ages[min_one_lineage]
  ltt_comp <- data.frame(lineages=c(ltt$lineages, dlineages),
                         ages=c(ltt$times, dages),
                         Type=c(rep("LTT", length(ltt$lineages)),
                                rep("dLTT", length(dlineages))))
  ltt_comp$ages <- max(ltt_comp$ages) - ltt_comp$ages
  
  ggplot(ltt_comp, aes(x=ages, y=lineages, colour=Type)) + 
    stat_smooth(se=FALSE) + scale_y_log10() + theme_cowplot() +
    scale_colour_manual(values = cols) 
    ##theme(legend.position = "none", axis.title = element_blank(),
    ##      axis.text = element_text(size=9))
}

tree_info <- read.csv("data/tree_descriptions.csv")

for (i in 1:nrow(tree_info)){
  id <- as.character(tree_info[i,"tree_name"])
  rho <- as.numeric(tree_info[i,"rho"])
  clade <- as.character(tree_info[i, "taxon"])
  tree <- castor::read_tree(file=paste0("data/trees/", id, ".tre"))
  fit <- pdr_data[[clade]]
  ltt_dltt(tree, fit)
  ggsave(paste0("figs/ltt/ltt-dltt-", clade, ".pdf"))
}


## Figure 7: epsilon vs. mu_0  

## 1. compute mu_0 from fitted rp and lambda_0
