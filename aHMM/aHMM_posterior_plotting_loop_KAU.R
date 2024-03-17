#Code for Generating a plot of Ancestry_HMM output
#Julia Harenčár
#

library(readr)
library("ggplot2")
library(dplyr)
theme_set(theme_bw())
library(cowplot)

## clear workspace
dev.off()
rm(list = ls())

# popgen wd
# setwd("/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/aHMM/final_params_0.02er_0.0001pulse")

# QTL wd
setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/full_final_posteriors/")

directory <- "./chr2/tmp/"
chromosome <- "chr2"

files <- list.files(directory, pattern = "*posterior")
IDs <- substr(files, 1, 9)

plot_list = list()
for (i in 1:length(files)) {
  post <- read_tsv(paste0(directory, files[i]))
  # rename cols 
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  # Filter out posteriors with a 0.9 threshold
  post <- post %>% filter(ref_homo > 0.9 | alt_homo > 0.9 | het > 0.9) 
  # Create allele frequency cols for plotting 
  post <- post %>% mutate (ancestry = ref_homo+het*0.5)
  p <- ggplot(post, aes(position, ancestry)) + 
    ggtitle(paste(IDs[i], chromosome)) +
    geom_point() +
    ylim(-0.01,1.01)
  plot_list[[i]] = p
}

pdf("chr2.QTL.final_params_0.02e_0.0001p.pdf", height = 25, width = 12)
# looking at just 15 from mapping pop
plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
          plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
          plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]],
          plot_list[[13]], plot_list[[14]],
          nrow = 9,
          ncol = 2)
dev.off()

# # for popgen (all 66):
# pdf("chr8.popgen.final_params_0.02e_0.0001p.pdf", height = 80, width = 12)
# plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
#           plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
#           plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]],
#           plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]],
#           plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]],
#           plot_list[[21]], plot_list[[22]], plot_list[[23]], plot_list[[24]],
#           plot_list[[25]], plot_list[[26]], plot_list[[27]], plot_list[[28]],
#           plot_list[[29]], plot_list[[30]],plot_list[[31]], plot_list[[32]],
#           plot_list[[33]], plot_list[[34]], plot_list[[35]], plot_list[[36]],
#           plot_list[[37]], plot_list[[38]], plot_list[[39]], plot_list[[40]],
#           plot_list[[41]], plot_list[[42]], plot_list[[43]], plot_list[[44]],
#           plot_list[[45]], plot_list[[46]], plot_list[[47]], plot_list[[48]],
#           plot_list[[49]], plot_list[[50]], plot_list[[51]], plot_list[[52]],
#           plot_list[[53]], plot_list[[54]], plot_list[[55]], plot_list[[56]],
#           plot_list[[57]], plot_list[[58]], plot_list[[59]], plot_list[[60]],
#           plot_list[[60]], plot_list[[62]], plot_list[[63]], plot_list[[64]],
#           nrow = 33,
#           ncol = 2)
# dev.off()
