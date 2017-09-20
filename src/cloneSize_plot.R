#setwd("C:/Users/user-1/Dropbox/Isotype switching/isotypeFrac_cloneSize")
library(ggplot2)
data = read.table("isotypeFrac_cloneSize.txt", header = TRUE, sep = ',')

#parameters
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

#ggplot theme
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

#plots of clone size vs fraction IgA
numGplot = ggplot(data, aes(x=numIgG, y=effFracIgA)) + geom_point() + facet_wrap(~dataset, scale="free") + ylab('IgA/(IgG+IgA)') + xlab('number of IgG')
pdf("numG_vs_fracA.pdf", height=5, width=12)
plot(numGplot)
dev.off()

effSizePlot = ggplot(data, aes(x=effSize, y=effFracIgA)) + geom_point() + facet_wrap(~dataset, scale="free") + ylab('IgA/(IgG+IgA)') + xlab('IgG+IgA')
pdf("effSize_vs_fracA.pdf", height = 5, width = 12)
plot(effSizePlot)
dev.off()

#pearson correlation tests for clone size and fraction IgA
cor_list = c()
cor_list = rbind(cor_list, cor.test(data011$effSize, data011$effFracIgA))
cor_list = rbind(cor_list, cor.test(data012$effSize, data012$effFracIgA))
cor_list = rbind(cor_list, cor.test(dataSFA$effSize, dataSFA$effFracIgA))
cor_list = rbind(cor_list, cor.test(data011$numIgG, data011$effFracIgA))
cor_list = rbind(cor_list, cor.test(data012$numIgG, data012$effFracIgA))
cor_list = rbind(cor_list, cor.test(dataSFA$numIgG, dataSFA$effFracIgA))


cor_list[1,][c("statistic", "p.value", "estimate", "method")]


statFile = file("cor_tests.txt", "a")
for (row in 1:nrow(cor_list)){
write(paste(cor_list[row,][c("estimate", "p.value", "statistic", "method")], collapse=","), statFile, append=TRUE)
}
close(statFile)

#spearman correlation tests for clone size and fraction IgA
s1 = cor.test(data011$effSize, data011$effFracIgA, method="spearman")
s2 = cor.test(data012$effSize, data012$effFracIgA, method="spearman")
s3 = cor.test(dataSFA$effSize, dataSFA$effFracIgA, method="spearman")
s4 = cor.test(data011$effSize, data011$effFracIgA, method="spearman")
s4 = cor.test(data011$numIgG, data011$effFracIgA, method="spearman")
s5 = cor.test(data012$numIgG, data012$effFracIgA, method="spearman")
s6 = cor.test(dataSFA$numIgG, dataSFA$effFracIgA, method="spearman")

spear_list = c()
spaer_list = rbind(s1,s2,s3,s4,s5,s6)

spearFile = file("spear_stats.txt", "a")
for(row in 1:nrow(spear_list)){
write(paste(spear_list[row,][c("estimate", "p.value", "statistic", "method")], collapse=','), spearFile, append=TRUE)
}
close(spearFile)



#Read frequency distribution of clone sizes
sizeData = read.table("frac_clone_with_effSize.txt", header=TRUE, sep=',')

#
sizePlot = ggplot() + geom_point(data=sizeData, aes(x=log2(effSize), y=log2(onlyIgA)), color='red') + geom_point(data=sizeData,aes(x= log2(effSize), y=log2(onlyIgG)), color='blue') + geom_point(data=sizeData, aes(x=log2(effSize), y=log2(onlyIgA_or_onlyIgG)), color = 'green') + facet_wrap(~dataset, scale="free") + ggplot_theme + scale_y_continuous(limits=c(-12.5, 0)) + scale_x_continuous(limits=c(0, 10)) + xlab("IgG+IgA") + ylab("fraction of clones")
pdf("sizeplot.pdf", height=5, width=12)
plot(sizePlot)
dev.off()


#binomial test for excess of singleton IgA
bn011 = binom.test(101, 379, 0.2028)
bn012 = binom.test(151, 464, 0.2456)
bnSFA = binom.test(153, 642, 0.1116)

bnFile = file("bn_test.txt", "a")
write(paste(bn011["p.value"]), bnFile, append=TRUE)
write(paste(bn012["p.value"]), bnFile, append=TRUE)
write(paste(bnSFA["p.value"]), bnFile, append=TRUE)
close(bnFile)