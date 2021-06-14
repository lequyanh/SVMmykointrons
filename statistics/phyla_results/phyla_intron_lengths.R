setwd("C:/Users/vuleq/Desktop/mykointron")
library(tidyverse)

fungi_phyla <- c('ascomycota', 
                 'basidiomycota', 
                 'blastocladiomycota', 
                 'cryptomycota', 
                 'chytridiomycota', 
                 'microsporidia', 
                 'mucoromycota', 
                 'zoopagomycota')

lens_file_suff <- '-intron-lens.txt'

# Load intron lens for each phyla into lists
intrn_lens_list <- lapply(fungi_phyla, function(x){
  lens_file <- paste0(x, lens_file_suff)
  
  lens_df <- read.table(lens_file)
  print(dim(lens_df))
  lens_df %>% mutate(phylum = x)
})

# Function for merging DF of different phyla together
my_merge <- function(df1, df2){
  rbind(df1, df2, by='phylum')
}

# Data frame with columns ## v1 | phylum ##
lens_data_all <- Reduce(my_merge, intrn_lens_list) %>%
  filter(phylum != 'phylum') %>%
  mutate(V1 = as.numeric(V1))

# Keep aside data for basidiomycota and ascomycota (for highlighting)
lens_data_asco_basi <- lens_data_all %>%
  filter(phylum %in% c('ascomycota', 'basidiomycota'))

lens_data_all %>% ggplot(aes(x=phylum, y=V1, fill=phylum)) +
  geom_boxplot() +
  ylim(0,300) +
  # geom_density(size=0.75) +
  # geom_density(aes(V1), data = lens_data_asco_basi, size=1.5) +
  # xlim(0, 400) +
  # xlab("Intron length") +
  # labs(title = "Intron length distribution across phyla") +
  theme_minimal()