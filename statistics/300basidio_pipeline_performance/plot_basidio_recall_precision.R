# Run batch version of cutting_diagnostics.py first to obtain the CSV

setwd("~/PycharmProjects/mycointrons/statistics/300basidio_pipeline_performance")
recall_precision <- read.csv("300basidio_recall_precision.csv", sep=';')

ordered_r <- recall_precision[order(recall_precision$recall), ]
ordered_r['portion'] <- c(1:length(ordered_r$recall)/length(ordered_r$recall))
plot(portion ~ recall, data=ordered_r)

recall_precision$exon_precision <- 1 - recall_precision$exon_breaking_fpr
ordered_p = recall_precision[order(recall_precision$exon_precision), ]
ordered_p['portion'] <- c(1:length(ordered_p$exon_precision)/length(ordered_p$exon_precision))
plot(portion ~ exon_precision, data=ordered_p)

plot(exon_precision ~ recall, data=recall_precision)
