# Script for plotting distribution of clustering index for IgG and IgA:

library(ggplot2)
library(cowplot)

results_01107PB <- read.table('../results/test_clustering/summary_stats_clustering_01107PB.csv', header =T, sep = ',')
results_01207PB <- read.table('../results/test_clustering/summary_stats_clustering_01207PB.csv', header =T, sep = ',')
results_SFAPB <- read.table('../results/test_clustering/summary_stats_clustering_SFAPB.csv', header =T, sep = ',')

dataframe_list = list('01107PB' = results_01107PB, '01207PB' = results_01207PB, 'SFAPB' = results_SFAPB)

dataset_vector <- c()
isotype_vector <- c()
clustering_index_vector <- c()
p_value_vector <- c()

# Condense dataframes into single dataframe for ggplot
for(dataset in c('01107PB', '01207PB', 'SFAPB')){
  
  data_frame <- dataframe_list[[dataset]]
  
  for(isotype in c('IgA', 'IgG')){
    
    obs_subtree_mean_size <- data_frame[,paste('mean_',isotype,'_subtree_size', sep = '')]
    rdm_subtree_mean_size <- data_frame[,paste('mean_',isotype,'_subtree_rdm', sep = '')]
    sd_subtree_size_rdm <- data_frame[,paste('SD_',isotype,'_subtree_rdm', sep = '')]
    p_value <- data_frame[, paste('p_', isotype, '_subtree_rdm', sep = '')]
    
    clustering_index <- (obs_subtree_mean_size - rdm_subtree_mean_size) / sd_subtree_size_rdm
    
    
    # Remove Inf values caused by clones with no variation in randomizations
    #clustering_index <- clustering_index[clustering_index != Inf]
    
    clustering_index_vector <- c(clustering_index_vector, clustering_index)
    dataset_vector <- c(dataset_vector, rep(dataset, length(clustering_index)))
    isotype_vector <- c(isotype_vector, rep(isotype, length(clustering_index)))
    p_value_vector <- c(p_value_vector, p_value)
  }
}

dataset_vector <- factor(dataset_vector, levels = c('01107PB', '01207PB', 'SFAPB'))
isotype_vector <- factor(isotype_vector, levels = c('IgA','IgG'))

combined_dataframe <- data.frame(dataset = dataset_vector, isotype = isotype_vector,
                                 clustering_index = clustering_index_vector, p_value = p_value_vector)


# ================ DEFINING GGPLOT2 PARAMETERS ===============
title_size <- 7
axis_title_size <- 11
y_axis_text_size <- 11
x_axis_text_size <- 11

axis_title_size_subplot <- 5
y_axis_text_size_subplot <- 5
x_axis_text_size_subplot <- 5

ylab_distance <- 7
xlab_distance <- 7
plot_margins <- c(0.1, 0.1, 0.1, 0.1)

lineage_names <- c(
  'CH103' = "CH103 (H)",
  'CH103L' = "CH103 (L)",
  'VRC26' = "VRC26 (H)",
  'VRC26L' = "VRC26 (L)",
  'VRC01_13' = 'VRC01-13 (H)',
  'VRC01_01' = 'VRC01-01 (H)',
  'VRC01_19' = 'VRC01-19 (H)'
)

ggplot_theme <- theme_bw() +
  theme(axis.title.y = element_text(size = axis_title_size,
                                    margin = margin(0,ylab_distance,0,0)),
        axis.title.x = element_text(size = axis_title_size,
                                    margin = margin(xlab_distance,0,0,0)),
        axis.text.x = element_text(size = x_axis_text_size, vjust = 0.6),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_line(size = 0.3),
        axis.ticks.length = unit(0.2, "cm"),
        #plot.margin = unit(c(1, 1, 1, 1),"cm"),
        title = element_text(size = title_size),
        strip.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = 'top',
        legend.text=element_text(size=16),
        legend.title = element_text(size = 16),
        panel.grid.minor = element_line(size = 0.1),
        panel.grid.major = element_line(size = 0.2)
  )

###### PLOTS #####

c_index_pl <- ggplot(combined_dataframe, aes(x = clustering_index, fill = isotype)) + 
    geom_density(alpha = 0.6) + ggplot_theme +
    facet_wrap(~dataset, scales = "free") + 
    scale_x_continuous(limits = c(-150,150)) + 
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = 0, linetype=2) + 
    xlab('Clustering index') + ylab('Density')

pdf('../figures/clustering_index_dist.pdf', height = 5, width = 12)
plot(c_index_pl)
dev.off()

p_value_pl <- ggplot(combined_dataframe, aes(x = p_value, fill = isotype)) + 
  geom_histogram(colour = 'black') + ggplot_theme +
  facet_grid(isotype~dataset) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,65)) +
  geom_vline(xintercept = 0.05, linetype=2) + 
  xlab('p value') + ylab('Frequency')



