setwd("C:/Users/vuleq/Desktop/mykointron")
library(tidyverse)
library(RColorBrewer)
library(reshape2)

fungi_phyla <- c('ascomycota', 
                 'basidiomycota', 
                 'blastocladiomycota', 
                 'cryptomycota', 
                 'chytridiomycota', 
                 'microsporidia', 
                 'mucoromycota',
                 'zoopagomycota')

file_suff <- '_metrics.csv'

phyla_metrics_list <- lapply(fungi_phyla, function(x){
  metrics_file <- paste0(x, file_suff)
  print(metrics_file)
  
  metrics_df <- read.csv(metrics_file, sep=';', 
                         colClasses=c('character','numeric','numeric','numeric','numeric','numeric', 'numeric'))
  metrics_df %>% mutate(phylum = x)
})

# Function for merging DF of different phyla together
my_merge <- function(df1, df2){
  bind_rows(df1, df2)
}

# Data frame with columns ## v1 | phylum ##
metrics_data_all <- Reduce(my_merge, phyla_metrics_list) %>%
  mutate()

# Add metrics such as Recall and Exon precision
metrics_data_all <- metrics_data_all %>% 
  mutate(recall_all = 100 * true_cuts / all_introns,
         recall_detectable = 100 * true_cuts / detectable_introns,
         dmg_exons = 100 * exon_cuts / exons)

################
# RECALL PLOTS #
################
metrics_recall <- metrics_data_all %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(fungi, phylum, recall_all, recall_detectable) %>%
  melt(id.vars = c('fungi', 'phylum'), measure.vars = c('recall_all', 'recall_detectable'))

give.n <- function(x){
  return(c(y = median(x) + 2, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

metrics_recall %>% ggplot(aes(x=phylum, y=value, fill=phylum)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  facet_grid(cols = vars(variable),
             labeller = as_labeller(c("recall_all" = "All introns",
                                      "recall_detectable" = "Detectable introns (range 40 - 100nt)"))) +
  ylab("Recall (%)") +
  ggtitle("Proportion of correctly cut introns") +
  scale_color_brewer(palette = "Dark2")

###################
# PRECISION PLOTS #
###################
metrics_exon_dmg <- metrics_data_all %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(fungi, phylum, dmg_exons)


metrics_exon_dmg %>% ggplot(aes(x=phylum, y=dmg_exons, fill=phylum)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ylab("Damaged exons (%)") +
  scale_color_brewer(palette = "Set3") +
  ggtitle("Damaged exons due to cuts")

####################
# Proportions plot #
####################
library(ggrepel)

metrics_gain <- metrics_data_all %>% 
  mutate(gain = recall_detectable / dmg_exons) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(gain = ifelse(gain != 0, gain, 0.025))

is_outlier <- function(x) {
  print(quantile(x, 0.75))
  return(x > quantile(x, 0.75) + 5 * IQR(x) | x < quantile(x, 0.25) - 1.8 * IQR(x))
}

metrics_gain <- metrics_gain %>% group_by(phylum) %>%
  mutate(outlier = ifelse (is_outlier(gain), fungi, NA)) 

GRAND_MEAN <- mean(metrics_gain$gain)

metrics_gain %>% ggplot(aes(x=phylum, y=gain, fill = phylum)) +
  geom_boxplot() +
  scale_y_continuous(trans = 'log2', breaks = c(0.125, 1, round(GRAND_MEAN,3), 8, 64)) +
  geom_hline(yintercept=GRAND_MEAN, color = "black") +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
  geom_text_repel(aes(label = outlier), na.rm = TRUE, hjust = -0.3) +
  ylab("Gain ratio") +
  scale_color_brewer(palette = "Set3") +
  ggtitle("Gain ratio (intron recall/exon damage)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "none")
