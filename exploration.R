
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)


# read data ---------------------------------------------------------------

final_hits = read_csv("final_hits.csv")

seqs_hits_extended = read_csv("seqs_hits_extended.csv")


seqs_hits_extended %>% 
  count(Species) %>%
  ggplot(aes(x = n)) +
  geom_histogram(stat = 'bin') +
  theme_light()



seqs_hits_extended %>% 
  count(Species) %>% view
seqs_hits_extended %>%
  distinct(Species, .keep_all = TRUE) %>% 
  write_csv('unique_seq_hits.csv')

