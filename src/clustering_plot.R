# Script for plotting distribution of clustering index for IgG and IgA:

library(ggplot2)
library(cowplot)

results_01107PB <- read.table('../results/test_clustering/summary_stats_clustering_01107PB.csv', header =T, sep = ',')
results_01207PB <- read.table('../results/test_clustering/summary_stats_clustering_01207PB.csv', header =T, sep = ',')
results_SFAPB <- read.table('../results/test_clustering/summary_stats_clustering_SFAPB.csv', header =T, sep = ',')

dataframe_list = list('01107PB' = results_01107PB, '01207PB' = results_01207PB, 'SFAPB' = results_01207PB)

dataset_vector <- c()
isotype_vector <- c()
clustering_idex_vector <- c()

# Condense dataframes into single dataframe for ggplot
for(dataset in c('01107PB', '01207PB', 'SFAPB')){
  
  data_frame <- dataframe_list[[dataset]]
  
  for(isotype in c('IgA', 'IgG')){
    
    obs_subtree_mean_size <- data_frame[,paste('mean_',isotype,'_subtree_size', sep = '')]
    rdm_subtree_mean_size <- data_frame[,paste('mean_',isotype,'_subtree_rdm', sep = '')]
    sd_subtree_size_rdm <- data_frame[,paste('SD_',isotype,'_subtree_rdm', sep = '')]
    
    clustering_index <- (obs_subtree_mean_size - rdm_subtree_mean_size) / sd_subtree_size_rdm
    
    # Remove Inf values caused by clones with no variation in randomizations
    clustering_index <- clustering_index[clustering_index != Inf]
    
    clustering_idex_vector <- c(clustering_idex_vector, clustering_index)
    dataset_vector <- c(dataset_vector, rep(dataset, length(clustering_index)))
    isotype_vector <- c(isotype_vector, rep(isotype, length(clustering_index)))
  }
}

combined_dataframe <- 