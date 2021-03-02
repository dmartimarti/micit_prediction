
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)
library(glue)
library(heatmaply)

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







# read freq table ---------------------------------------------------------

full_df = read_excel("species_freqs.xlsx", 
                            sheet = "FULL") %>% 
  select(-(`...1`:`Prevalence in Paraf. Controls`)) %>% 
  unite(ID, phylum,class,order,family,genus,species, sep='_', remove = F) %>% 
  rename(index = `Unnamed: 2`)


hits_df = read_excel("species_freqs.xlsx", 
                     sheet = "Hits") %>% 
  select(-(`...1`:`Prevalence in Paraf. Controls`)) %>% 
  unite(ID, phylum,class,order,family,genus,species, sep='_', remove = F) %>% 
  rename(index = `Unnamed: 2`)

hits_matrix = read_excel("aay9189_TableS4.xlsx", 
                              sheet = "Hits", skip = 2) %>% 
  rename(index = `...1`)


full_df
hits_df

# name of samples
samples = unique(str_split(names(full_df),pattern = '[.]',simplify = T)[,1])[10:24]
 

# test sum
var = str_replace(samples[1],pattern = ' ',replacement = '_')

full_df %>% 
  select(index,ID, starts_with(samples[1])) %>% 
  rowwise(index,ID) %>% 
  mutate(cosa = sum(c_across(starts_with(samples[1]))),.before = samples[1])



# helper function
sum_by = function(data, by, var) {
  suffix = var
  data %>%
    select(index,ID, starts_with({{ var }})) %>% 
    rowwise({{ by }}) %>%
    mutate('{suffix}_sum' := sum(c_across(starts_with({{ var }}))), # walrus op
           .before = {{ var }}) %>% 
    ungroup %>% 
    select(index,ID,contains('_sum'))
}

# iterate over full df
end_df = full_df %>% select(index,ID)
for (sample in samples){
  print(glue('Processing sample {sample}'))
  end_df = end_df %>% 
    left_join(full_df %>% 
                sum_by(ID,sample)) 
  
}

# drop NA values (1 row only)
end_df = end_df %>%
  drop_na(`Breast (NAT)_sum`)

### correct by hits ####
## VERY IMPORTANT
# now we need to multiply the values by the "hit" matrix from
# supplementary table 4. 1 is a hit and the bug passed the quality
# filter, and 0 it's a contaminant

## calculate the number of hits per tissue
n_hits = hits_matrix %>% 
  summarise(across(contains('('), ~sum(.x))) %>% 
  as.data.frame 

n_hits[2,] = colnames(n_hits)

colnames(n_hits) = NULL

n_hits = t(n_hits)

colnames(n_hits) = c('hits', 'tissue')

n_hits %>% 
  as_tibble %>%
  mutate(hits=as.numeric(hits)) %>% 
  ggplot(aes(x = hits, y = tissue)) +
  geom_point() +
  theme_light()


n_hits['tissue'] = rownames(n_hits)




fix_end_df = end_df %>% 
  left_join(hits_matrix)

names(hits_matrix)

# multiply by hits
fix_end_df = fix_end_df %>% 
  mutate(`Breast (NAT)_sum` = `Breast (NAT)_sum`*`Breast (NAT)`,
         `Breast (N)_sum` = `Breast (N)_sum` * `Breast (N)`,
         `Breast (N-FA)_sum` =`Breast (N-FA)_sum` * `Breast (N)`,
         `Breast (T)_sum` = `Breast (T)_sum` * `Breast (T)`,
         `Lung (NAT)_sum` = `Lung (NAT)_sum` * `Lung (NAT)`,
         `Lung (T)_sum` = `Lung (T)_sum` * `Lung (T)`,
         `Melanoma (T)_sum` = `Melanoma (T)_sum` * `Melanoma (T)`,
         `Pancreas (T)_sum` = `Pancreas (T)_sum` * `Pancreas (T)`,
         `Ovary (NAT)_sum` = `Ovary (NAT)_sum` * `Ovary (N+NAT)`,
         `Ovary (N- fallop)_sum` = `Ovary (N- fallop)_sum` * `Ovary (N+NAT)`,
         `Ovary (T)_sum` = `Ovary (T)_sum` * `Ovary (T)`,
         `Colon (NAT)_sum` = `Colon (NAT)_sum` * `Colon (NAT)`,
         `Colon (T)_sum` = `Colon (T)_sum` * `Colon (T)`,
         `Bone (T)_sum` = `Bone (T)_sum` * `Bone (T)`,
         `GBM (T)_sum` = `GBM (T)_sum` * `GBM (T)`)


fix_end_df = fix_end_df[,1:17]

fix_end_df = fix_end_df %>%
  drop_na(`Breast (NAT)_sum`)


# write 
write_csv(fix_end_df, 'total_tissue_reads_sum.csv')


# calculate freqs in the big table

total_freqs = fix_end_df

# calculate frequencies 
# helping function
scale2 = function(x, na.rm = TRUE) (x) / sum(x)

# 
total_freqs = total_freqs %>% 
  mutate(across(contains('('),scale2))

# this is the total freqs without 
global_freqs = total_freqs %>% filter(index %in% hits_df$index)

colSums(global_freqs[,3:17]) %>% 
  tibble(sample = names(colSums(global_freqs[,3:17]))) %>% 
  ggplot(aes(x = ., y = sample)) + geom_point() +
  labs(y = 'Sample',
       x = 'Total frequency of hits') +
  theme_light()


ggsave('hits_global_freqs_v2.pdf', height = 7,width = 7)

# and now the freqs within the hits

hits_freqs = fix_end_df %>% filter(index %in% hits_df$index)
# calculate frequencies 

hits_freqs = hits_freqs %>% 
  mutate(across(contains('('),scale2))

# add some useful columms
hits_freqs = hits_df %>% select(index:species) %>% 
  left_join(hits_freqs)


list_of_datasets = list(
  raw_values = end_df,
  total_freqs = total_freqs,
  hits_global_freqs = global_freqs,
  hits_freqs = hits_freqs
)

write.xlsx(list_of_datasets,'freqs_hits.xlsx')





# micit prod by tissue ----------------------------------------------------

library(readxl)
micit = read_excel("MICIT2_selected_strains.xlsx", 
                                      sheet = "sim_table_ngm", skip = 1)
names(micit)

micit = micit %>% 
  select(species = Species, NGM = NGM_MicitOptimization, serum = serum_MicitOptimization)

micit_sp = micit$species



hits_freqs_micit = hits_freqs %>% 
  filter(species %in% micit_sp) %>% 
  select(-(index:genus)) 



A = as.matrix(hits_freqs_micit[,2:16])

B = as.matrix(micit[,2:3])

micit_prod = t(A) %*% B


micit_prod = micit_prod %>% 
  as.data.frame %>% 
  tibble(tissue = rownames(micit_prod),validate = NULL)


micit_prod %>% 
  pivot_longer(NGM:serum, names_to = 'Media', values_to = 'Micit_prod') %>% 
  ggplot(aes(x = Micit_prod, y = reorder(tissue, Micit_prod, mean), fill=Media, group=Media)) +
  geom_bar(stat = 'identity') + 
  labs(x = 'Micit production',
       y = 'Tissue') +
  facet_wrap(~Media, scales = 'free_x') +
  theme_light() 

ggsave('Micit_prod_byTissue.pdf', height = 8, width = 10)

micit_prod %>% 
  filter(str_detect(tissue, '\\(T\\)')) %>% 
  pivot_longer(NGM:serum, names_to = 'Media', values_to = 'Micit_prod') %>% 
  ggplot(aes(x = Micit_prod, y = reorder(tissue, Micit_prod, mean), fill=Media, group=Media)) +
  geom_bar(stat = 'identity') + 
  labs(x = 'Micit production',
       y = 'Tissue') +
  facet_wrap(~Media, scales = 'free_x') +
  theme_light() 

ggsave('Micit_prod_byTissue_cancer.pdf', height = 8, width = 10)




# correlation 
# plot
library(broom)
micit %>% 
  filter(!(NGM == 0 & serum == 0)) %>% 
  ggplot(aes(NGM, serum)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  theme_light()

model = micit %>% 
  filter(!(NGM == 0 & serum == 0)) %>% 
  nest(everything()) %>% 
  mutate(model = map(data, lm, formula = 'NGM ~ serum'))  %>% 
  select(model)

glance((model$model[[1]]))
