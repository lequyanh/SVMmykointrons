setwd('/home/anhvu/PycharmProjects/mycointrons/statistics/intron-lenghts-statistics/')
#setwd('C:/Users/AnhVu/Desktop/myko')
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
  
  # tp_vs_fp <- sapply(20:200, function(x) {
  #   introns_dataset %>% filter(len < x) %>% 
  #     group_by(label) %>% 
  #     summarise(count = n()) %>%
  #     select(count)
  # })
  
  # temp <- as.data.frame(tp_vs_fp) %>% t() %>% as.data.frame()
  # temp <- temp %>% mutate(candidate_lens = 20:200)
  # 
  # ggplot(temp, aes(x=candidate_lens)) +
  #   geom_line(aes(y=V1)) +
  #   geom_line(aes(y=V2))
}

basidio_i_lens <- read.table('../../pipeline/basidiomycota-intron-lens.txt')
asco_i_lens <- read.table('../../pipeline/ascomycota-intron-lens.txt')

fungi_files <- read.table('fungi_list.txt')

# Get list of labeled intron data for each fungi in file
# Each DF has # scaffold # start # end # sequence # label #
# Remove unused columns and compute lenght of each sequence
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
quantile((intron_data_long %>% filter(label == 1))$len, c(0.05, 0.10, 0.9, 0.95), na.rm=T)
ecdf((intron_data_long %>% filter(label == 1))$len)(40)
ecdf((intron_data_long %>% filter(label == 1))$len)(100)
