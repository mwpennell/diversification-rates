## Figures for MS and suppmat
library(ggplot2)
library(cowplot)

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

lambda0ml <- sapply(pdr_data, function(x) return(x$fit$fitted_rholambda0))
lambda0lo <- sapply(pdr_data, function(x) return(x$fit$CI95lower$rholambda0))
lambda0up <- sapply(pdr_data, function(x) return(x$fit$CI95upper$rholambda0))
lambda0se <- sapply(pdr_data, function(x) return(x$fit$standard_errors$rholambda0))
age <- sapply(pdr_data, function(x) return(max(x$fit$age_grid)))

lambda0_df <- data.frame(lambda0rho = lambda0ml, rho = rhos, 
                         se=lambda0se, age=age, clade=clade_names)
lambda0_df$lambda0 <- lambda0_df$lambda0rho / lambda0_df$rho

## Exclude bacteria since it is so much older
lambda0_df <- lambda0_df[-which(lambda0_df$clade == "Bacteria"), ]

ggplot(lambda0_df, aes(x=age, y=lambda0)) +  
  geom_errorbar(aes(ymin=lambda0-se, ymax=lambda0+se)) +
  geom_point(colour="#5475BA", size=3, alpha=0.7) + 
  ylab(expression(lambda[p](0))) + xlab("Age of clade (my)") + 
  theme_cowplot() #+ stat_smooth(se=FALSE, colour="#2E3C61")
ggsave("figs/clade_age/lambda0_cladeage.pdf")


## Plot 4: rp0 vs clade age (standard scale)
## Exclude bacteria since it is so much older than the rest

rp0 <- sapply(pdr_data, function(x) return(x$fit$fitted_PDR[1]))
se <- sapply(pdr_data, function(x) return(x$fit$standard_errors$PDR[1]))
age <- sapply(pdr_data, function(x) return(max(x$fit$age_grid)))

df <- data.frame(rp0=rp0, se = se, age=age, clade=clade_names)
df <- df[-which(df$clade == "Bacteria"), ]

ggplot(df, aes(x=age, y=rp0)) +  
  geom_errorbar(aes(ymin=rp0-se, ymax=rp0+se)) +
  geom_point(colour="#5475BA", size=3, alpha=0.7) + 
  ylab(expression(r[p](0))) + xlab("Age of clade (my)") + 
  theme_cowplot() + stat_smooth(se=FALSE, colour="#2E3C61")
ggsave("figs/clade_age/rp0_cladeage.pdf")

## Plot 5: rp0 vs clade age (log scale)
## Exclude bacteria since it is so much older than the rest

ggplot(df, aes(x=age, y=rp0)) +  
    #geom_errorbar(aes(ymin=rp0-se, ymax=rp0+se)) +
    geom_point(colour="#5475BA", size=3, alpha=0.7) + 
    ylab(expression(r[p](0))) + xlab("Age of clade (my)") + 
    theme_cowplot() + stat_smooth(method="lm", se=FALSE, colour="#2E3C61") +
    scale_y_log10() + scale_x_log10()
ggsave("figs/clade_age/rp0_cladeage_log.pdf")