# Script for plotting the relationships between divergence and IgG/IgA composition:

library(ggplot2)
library(cowplot)

results_01107PB <- read.table('../results/divergence_vs_composition_01107PB.csv', header =T, sep = ',')
results_01207PB <- read.table('../results/divergence_vs_composition_01207PB.csv', header =T, sep = ',')
results_SFAPB <- read.table('../results/divergence_vs_composition_SFAPB.csv', header =T, sep = ',')

dataframe_list = list('01107PB' = results_01107PB, '01207PB' = results_01207PB, 'SFAPB' = results_SFAPB)

dataset_vector <- c()
isotype_vector <- c()
clustering_index_vector <- c()
p_value_vector <- c()

# Condense dataframes into single dataframe for ggplot
combined_dataframe <- rbind(results_01107PB, results_01207PB, results_SFAPB)

dataset <- c(rep('01107PB', nrow(results_01107PB)), rep('01207PB', nrow(results_01207PB)), 
                    rep('SFAPB', nrow(results_SFAPB)))

dataset <- factor(dataset, levels = c('01107PB','01207PB','SFAPB'))

combined_dataframe <- cbind(dataset, combined_dataframe)

rm(dataset)

# Classify clone into IgA but not IgG, IgG but not IgA, or both IgG and IgA:
clone_type <- c()

for(row in 1:nrow(combined_dataframe)){
  n_iga = combined_dataframe$n_iga[row]
  n_igg = combined_dataframe$n_igg[row]
  
  if(n_iga > 0 & n_igg == 0){
    type <- 'IgA_no_IgG'
  }
  if(n_igg > 0 & n_iga == 0){
    type <- 'IgG_no_IgA'
  }
  if(n_igg > 0 & n_iga > 0){
    type <- 'IgG_and_IgA'
  }
  if(n_igg == 0 & n_iga == 0){
    type <- 'no_IgG_no_IgA'
  }
  
  clone_type <- c(clone_type, type)
}

clone_type <- factor(clone_type, levels = c('IgA_no_IgG', 'IgG_no_IgA', 'IgG_and_IgA', 'no_IgG_no_IgA'))
combined_dataframe <- cbind(combined_dataframe, clone_type)
rm(clone_type)


#================ DEFINING GGPLOT2 PARAMETERS ===============
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
# Mean divergence vs. effective fraction of IgA (i.e., IgA/(IgG + IgA))
mean_div_vs_frac_IgA <- ggplot(combined_dataframe, aes(x = mean_divergence, y = n_iga/(n_iga + n_igg))) +
  geom_point() + ggplot_theme + facet_wrap(~dataset, scales = "free") +
  ylab('IgA / (IgG + IgA)') + xlab ('Mean divergence from germline')

pdf('../figures/mean_divergence_vs_fraction_iga.pdf', height = 5, width = 12)
plot(mean_div_vs_frac_IgA)
dev.off()

# Divergence boxplot (as function of clone type IgG-no-IgA, IgA-no-IgG, IgG-and-IgA, no-IgG-no-IgA)
divergence_vs_clone_type <- ggplot(combined_dataframe, aes(x = clone_type, y = mean_divergence)) +
  geom_boxplot() + ggplot_theme + facet_wrap(~dataset, scales = "free") +
  xlab('Clone type') + ylab('Mean divergence from germline')

pdf('../figures/divergence_vs_clone_type.pdf', height = 7, width = 15)
plot(divergence_vs_clone_type)
dev.off()

aov_results <- aov(data=combined_dataframe, formula = mean_divergence~clone_type + dataset + dataset*clone_type)
TukeyHSD(aov_results, "clone_type")


# Divergence vs. size separately for IgG-no-IgA, IgA-no-IgG, IgG-and-IgA, no-IgG-no-IgA
div_vs_size <- ggplot(subset(combined_dataframe, clone_type != 'XXXXXX'), 
                      aes(x = log10(n_sequences), y = mean_divergence, colour = clone_type)) +
  geom_point(shape = 1) + 
  ggplot_theme +
  facet_wrap(~dataset, scales = "fixed") + 
  geom_smooth(method = 'lm', se = FALSE) +
  scale_x_continuous(limits=c(0.5,4)) +
  scale_color_brewer(type = 'qual',palette=6) +
  xlab("Number of sequences (log10)") + ylab("Mean divergence from germline")

pdf('../figures/divergence_vs_size.pdf', height = 5, width = 12)
plot(div_vs_size)
dev.off()




