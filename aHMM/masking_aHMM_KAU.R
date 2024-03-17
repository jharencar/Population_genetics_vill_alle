# Kate's loop to list sites for masking plus my loop to filter the panels
# JGH file setup 

library(readr)
library("ggplot2")
library(dplyr)
theme_set(theme_bw())
library(cowplot)

## clear workspace
rm(list = ls())

## notin operator
`%notin%` <- Negate(`%in%`)

setwd("./aHMM/F1_posteriors_for_masking/")

## list all directories in working directory
directories <- list.dirs(path = ".")
## keep only the directories that contain posterior files
posterior_directories <- directories[grep("chr", directories)]

chrindsite <- c()
for (i in 1:length(posterior_directories)) {
  ## check to see if you have an incorrect number of posterior directories
  if(length(posterior_directories) != 9) {stop("You have more than 9 posterior directories")}
  ## iterate through posterior directories 
  posterior_directory <- posterior_directories[i]
  ## get posterior file names from directory
  posteriors <- list.files(posterior_directory, pattern = "*posterior")
  ## iterate through posterior files
  for (j in 1:length(posteriors)) {
    ## check to see if there are files that don't end in *.posterior
    if (length(grep(".posterior", posteriors)) != length(posteriors)) {stop("one of your posterior files doesn't look like a posterior file")}
    ## get IDs from file name
    ID <- substr(posteriors[j], 1, 11)
    ## get ancestry group from IDs
    group <- substr(ID, 8, nchar(ID))
    # load posterior file
    post <- read_tsv(paste0(posterior_directory, "/", posteriors[j]))
    # rename cols 
    names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
    # Filter out posteriors with a 0.9 threshold
    post <- post %>% filter(ref_homo > 0.9 | alt_homo > 0.9 | het > 0.9)
    # Create allele frequency column
    post <- post %>% mutate (ancestry = ref_homo+het*0.5)
    ## assign misperforming sites in lasius individuals 
    if (group == "ALLE") {
      post <- post %>% mutate(misperforming = case_when(ref_homo <= 0.9 ~ TRUE,
                                                        ref_homo > 0.9 ~ FALSE))
    }
    ## assign misperforming sites in bracteatus individuals
    if (group == "VILL") {
      post <- post %>% mutate(misperforming = case_when(alt_homo <= 0.9 ~ TRUE,
                                                        alt_homo > 0.9 ~ FALSE))
    }
    ## assign misperforming sites in F1 hybrid individual(s)
    if (group == "HYB.") {
      post <- post %>% mutate(misperforming = case_when(het <= 0.9 ~ TRUE,
                                                        het > 0.9 ~ FALSE))
    }
    ## if there is at least one misperforming site, add those sites to the misperforming_sites matrix
    if (length(which(post$misperforming == TRUE)) > 0) {
      post_misperforming <- post %>% 
        filter(misperforming == TRUE) %>%
        select(chrom, position)
      ## if the misperforming_sites matrix doesn't exist, create it
      if(exists("misperforming_sites") == FALSE) {
        misperforming_sites <- post_misperforming
      }
      ## if the misperforming_sites matrix already exists, append new sites to it
      else if (exists("misperforming_sites") == TRUE) {
        misperforming_sites <- rbind(misperforming_sites, post_misperforming)
      }
      ## a vector that provides more information about the misperforming sites
      chrindsite <- c(chrindsite, paste("chromosome:", i, ", ID:", ID, ", site:", post$position[which(post$misperforming == TRUE)]))
    }
  }
}

## remove any duplicated rows
misperforming_sites_unique <- misperforming_sites %>% 
  distinct() %>%
  arrange(chrom, position)

## Count the number of misperforming sites per chromosome
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom1",]) #(KAU-73), 268 (27014, 26746)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom2",]) #0 , 551 (30552, 30001)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom3",]) #0 , 685 (21082, 20397)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom4",]) #0 , 376 (18308, 17932)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom5",]) #0 , 767 (25374, 24607)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom6",]) #533, 0 (19823,19823)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom7",]) #XX, 459 (17887, 17428)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom8",]) #160, 1,751 (19321, 17570)
nrow(misperforming_sites_unique[misperforming_sites_unique$chrom == "Chrom9",]) #92, 496 (17508, 17012)

# ## Plotting functions in case you need to investigate... 
# ggplot(post[post$position<2.5e7,], aes(position, ancestry)) + 
#   #ggtitle(paste(IDs[i], chromosome)) +
#   geom_point() +
#   ylim(0,0.55)
# 
# ggplot(post, aes(position, ancestry)) + 
#   #ggtitle(paste(IDs[i], chromosome)) +
#   geom_point() +
#   ylim(0,0.55)

## use new list of misperforming sites to filter panel files
for (i in 1:length(posterior_directories)) {
  posterior_directory <- posterior_directories[i]
  # filter out rows with values for V1 and V2 matching values of chrom and position in the misperforming sites DF
  panel_file <- list.files(posterior_directory, pattern = "*_QTL.panel") #change for other panel file names, first time was just *panel
  panel <- read.table(paste0(posterior_directory,"/",panel_file))
  panel_name <- paste0("panel_chr", i)
  filtered_panel <- panel %>%
  anti_join(misperforming_sites_unique, by = c("V1" = "chrom", "V2" = "position"))
  assign(panel_name, filtered_panel)
  file_name <- paste0("chr", i, "_QTL_filtered.panel")
  write.table(get(paste0("panel_chr", i)), file = file_name, sep = "\t", row.names = FALSE, col.names=FALSE)
}
