## Capture the output and error logs
output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

library(dplyr)
library(vroom)
library(purrr)
library(tidyr)

BEST_FILES <- snakemake@input[[1]]

map_dfr(1:length(BEST_FILES), function(i) {
  
  df <- vroom(BEST_FILES[i])
  df$pool <- gsub(".best", "", basename(BEST_FILES[i]))
  df
  
}) %>%
  mutate(prediction = sapply(strsplit(BEST, split = "-"), "[[", 1)) -> df

# Deal with single pool prediction by demuxlet. For pools with unique samples, demuxlet predicts all barcodes as ambiguous. 
# In that case, set all barcodes as singlets and let the preprocessing pipeline handle the homotypic doublets.
single_pools <- df %>% 
  select(pool, sample = SNG.1ST) %>% 
  unique() %>% 
  group_by(pool) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) %>% 
  pull(pool)

if(length(single_pools) > 0) {
  for(p in single_pools) {
    df <- df %>% 
      mutate(prediction = ifelse(pool == p, "SNG", prediction))
  }
}

singlets <- df %>% 
  filter(prediction == "SNG") %>% 
  select(BARCODE, sample = SNG.1ST, pool) %>% 
  mutate(prediction = "singlet")

doublets <- df %>% 
  filter(prediction == "DBL") %>% 
  mutate(sample = paste(DBL.1ST, DBL.2ND, sep = "&")) %>%
  select(BARCODE, sample, pool) %>% 
  separate_rows(sample, sep = "&") %>% 
  mutate(prediction = "doublet")

all_df <- singlets %>% 
  bind_rows(doublets)

write.csv(all_df, file = snakemake@output[[1]], row.names = F, quote = F)

sink()
sink()