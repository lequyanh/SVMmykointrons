setwd('/home/anhvu/PycharmProjects/mycointrons/statistics/intron-lenghts-statistics/')
library(dplyr)
library(ggplot2)

my_merge <- function(df1, df2){
  rbind(df1, df2, by = "species")
}

get.intron.lens <- function(data_filename) {
  fungi <- gsub('_introns.csv', '', data_filename)
  
  file <- paste0('basi_introns_dataset','/', data_filename)
  introns_dataset <- read.csv(file[1], header=T, sep=";")
  
  print(fungi)
  lens <- apply(introns_dataset, MARGIN = 1, function(x){nchar(x[4])})
  introns_lens_dataset <- introns_dataset %>% 
    mutate(len = lens,
           species = fungi) %>%
    select(species, len, label)
}

basidio_i_lens <- read.table('./phyla_lengths/basidiomycota-intron-lens.txt')
asco_i_lens <- read.table('./phyla_lengths/ascomycota-intron-lens.txt')

fungi_files <- read.table('fungi_list.txt')

# Get list of labeled intron data for each fungi in file
# Each DF has # scaffold # start # end # sequence # label #
# Remove unused columns and compute length of each sequence
intron_lens <- lapply(fungi_files$V1, get.intron.lens)

# Transform the list into long-form DF (for plotting)
# Columns # species # label # len #
intron_data_long <- Reduce(my_merge, intron_lens)
intron_data_long$len <- as.numeric(intron_data_long$len) 

# Plot distribution of positive introns
i_distr_plot <- ggplot(intron_data_long %>% filter(label == 1), aes(x=len, color=species)) +
  geom_density() +
  geom_density(aes(V1), data = basidio_i_lens, size=1, color='#050517') + 
  geom_density(aes(V1), data = asco_i_lens, size=1, color='#38040E') + 
  xlim(0, 200) +
  theme_minimal() +
  xlab("Candidate length") + labs(title = "Positive intron candidates length distribution ")
i_distr_plot

ggsave(filename = 'i_distr_plot.png', plot = i_distr_plot)

# Print quantiles of the distribution
quantile((intron_data_long %>% filter(label == 1))$len, c(0.025, 0.05, 0.95, 0.975), na.rm=T)
ecdf((intron_data_long %>% filter(label == 1))$len)(40)
ecdf((intron_data_long %>% filter(label == 1))$len)(100)

###################################
# Process phyla_intron_gene_stats #
###################################

fungi_phyla <- c('ascomycota', 
                 'basidiomycota',
                 'blastocladiomycota', 
                 'cryptomycota', 
                 'chytridiomycota', 
                 'microsporidia', 
                 'mucoromycota',
                 'zoopagomycota')

summary_table <- sapply(fungi_phyla, function(phylum) {
  stats_df <- read_csv(sprintf('./phyla_introns_genes_stats/%s_intron_gene_stats.csv', phylum))
  # hist(stats_df$intron_count, main = phylum)
  return(list(median_gene_count = median(stats_df$gene_count),
              mean_introns_per_gene = mean(stats_df$median_introns_per_gene),
              median_intron_count = median(stats_df$intron_count)))
}) %>% t()

summary_table
