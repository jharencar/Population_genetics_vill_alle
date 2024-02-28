#Code for Generating a plot of Ancestry_HMM output
#Julia Harenčár
#

library(readr)
library("ggplot2")
library(dplyr)
theme_set(theme_bw())
library(cowplot)

setwd("~/Google Drive/My Drive/Costus/genetic_mapping/local_ancestry/aHMM/ancestryInfer6/")

directory <- "./chr9/"
chromosome <- "scaffold_9"

files <- list.files(directory, pattern = "*posterior")
IDs <- substr(files, 1, 11)

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
    ylim(0,1)
  plot_list[[i]] = p
}  

pdf("chr9.plots.pdf", height = 50)
plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], 
          plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
          plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]],
          plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]],
          plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]],
          plot_list[[21]], plot_list[[22]], plot_list[[23]], plot_list[[24]],
          plot_list[[25]], plot_list[[26]], plot_list[[27]], plot_list[[28]],
          plot_list[[29]], plot_list[[30]],
          nrow = 15,
          ncol = 2)
dev.off()