library(ggplot2)
library(cowplot)

results_01107PB <- read.table('../results/divergence_vs_composition_01107PB.csv', header =T, sep = ',')
results_01207PB <- read.table('../results/divergence_vs_composition_01207PB.csv', header =T, sep = ',')
results_SFAPB <- read.table('../results/divergence_vs_composition_SFAPB.csv', header =T, sep = ',')

dataframe_list = list('01107PB' = results_01107PB, '01207PB' = results_01207PB, 'SFAPB' = results_SFAPB)


dataset <- c(rep('01107PB', nrow(results_01107PB)), rep('01207PB', nrow(results_01207PB)), 
             rep('SFAPB', nrow(results_SFAPB)))

dataset <- factor(dataset, levels = c('01107PB','01207PB','SFAPB'))

combined_dataframe <- rbind(results_01107PB, results_01207PB, results_SFAPB)
combined_dataframe <- cbind(dataset, combined_dataframe)

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

### PLOTS
# Clone size vs. mean divergence from UCA
size_vs_mean_divergence_pl <- ggplot(subset(combined_dataframe, fraction_igg >= 0.9), 
                                     aes(x = mean_divergence, y = n_sequences)) +
  geom_point() + facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  xlab('Mean divergence from UCA') + ylab('Clone size')

pdf('../figures/exploratory_plots/size_vs_mean_divergence.pdf', height = 5, width = 12)
plot(size_vs_mean_divergence_pl)
dev.off()


# Mean divergence vs mean S5F of UCA
mean_divergence_vs_UCA_S5F_pl <- ggplot(subset(combined_dataframe, fraction_igg >= 0.9), 
                                        aes(x = naive_S5F,y = mean_divergence)) +
  geom_point() + facet_wrap(~dataset) + ggplot_theme +
  xlab('Mean S5F-mutability of UCA') + ylab('Mean divergence from UCA')

pdf('../figures/exploratory_plots/mean_divergence_vs_UCA_S5F.pdf', height = 5, width = 12)
plot(mean_divergence_vs_UCA_S5F_pl)
dev.off()

# Mean divergence vs number of hotspots in UCA
mean_divergence_vs_UCA_HS_pl <- ggplot(subset(combined_dataframe, fraction_igg >= 0.9), 
                                        aes(x = naive_HS,y = mean_divergence)) +
  geom_point() + facet_wrap(~dataset) + ggplot_theme +
  xlab('Number of hotspots in UCA') + ylab('Mean divergence from UCA')

pdf('../figures/exploratory_plots/mean_divergence_vs_UCA_HS.pdf', height = 5, width = 12)
plot(mean_divergence_vs_UCA_HS_pl)
dev.off()


# Clone size vs. mean S5F of UCA
size_vs_UCA_S5F_pl <- ggplot(subset(combined_dataframe, fraction_igg >= 0.9), 
                             aes(x = naive_S5F,y = n_sequences)) +
  geom_point() + facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  xlab('Mean S5F-mutability of UCA') + ylab('Clone size')

pdf('../figures/exploratory_plots/size_vs_UCA_S5F.pdf', height = 5, width = 12)
plot(size_vs_UCA_S5F_pl)
dev.off()

# Clone size vs. number of hotspots in UCA
size_vs_UCA_HS_pl <- ggplot(subset(combined_dataframe, fraction_igg >= 0.9), 
                             aes(x = naive_HS,y = n_sequences)) +
  geom_point() + facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  xlab('Number of hotspots in UCA') + ylab('Clone size')

pdf('../figures/exploratory_plots/size_vs_UCA_HS.pdf', height = 5, width = 12)
plot(size_vs_UCA_HS_pl)
dev.off()


# Distribution of clone sizes
clone_size_distribution_pl <- ggplot(subset(combined_dataframe), 
                                          aes(x = n_sequences)) +
  geom_histogram(fill = 'grey70', colour = 'black') + 
  #geom_density() +
  facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  #scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab ('Clone size') +
  ylab ('Frequency')

pdf('../figures/exploratory_plots/clone_size_distribution.pdf', height = 5, width = 12)
plot(clone_size_distribution_pl)
dev.off()

# Distribution of clone sizes (log-log)
clone_size_distribution_loglogscale_pl <- ggplot(subset(combined_dataframe), 
                                     aes(x = n_sequences)) +
  geom_histogram(fill = 'grey70', colour = 'black') + 
  facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  scale_y_log10(expand = c(0,0)) +
  scale_x_log10() +
  xlab ('Clone size') +
  ylab ('Frequency')

pdf('../figures/exploratory_plots/clone_size_distribution_loglogscale_pl.pdf', height = 5, width = 12)
plot(clone_size_distribution_loglogscale_pl)
dev.off()

# Distribution of mean divergence from UCA across clones:
mean_divergence_distribution_pl <- ggplot(subset(combined_dataframe), 
                                      aes(x = mean_divergence)) +
  geom_histogram(fill = 'grey70', colour = 'black') + 
  facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab ('Mean divergence from UCA') +
  ylab ('Frequency')
  
pdf('../figures/exploratory_plots/mean_divergence_distribution.pdf', height = 5, width = 12)
plot(mean_divergence_distribution_pl)
dev.off()

# Distribution of mean divergence from UCA across clones (log x axis):
mean_divergence_distribution_logscale_pl <- ggplot(subset(combined_dataframe), 
                                          aes(x = mean_divergence)) +
  geom_histogram(fill = 'grey70', colour = 'black') + 
  facet_wrap(~dataset, scales = 'free') + ggplot_theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_log10() +
  xlab ('Mean divergence from UCA') +
  ylab ('Frequency')

pdf('../figures/exploratory_plots/mean_divergence_distribution_logscale.pdf', height = 5, width = 12)
plot(mean_divergence_distribution_logscale_pl)
dev.off()

  









