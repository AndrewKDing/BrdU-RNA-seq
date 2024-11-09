library(readxl)
library(tidyverse)

animals <- read_xlsx("./data/library_info.xlsx") %>%
  mutate(sample = `Sample no.`, 
         animal = factor(Animal), 
         timepoint = factor(Timepoint, levels = c("21dpi","nec")),
         cd8_depleted = factor(cd8_depleted, levels = c("N","Y")),
         infected = factor(infected, levels = c("N","Y")),
         aids = factor(aids, levels = c("N","Y")),
         sive = factor(sive, levels = c("N","Y")),
         plus_minus = factor(`plus/minus`, levels = c("minus","plus")),
         library_date,
         .keep = "none") %>%
  #relocate(sample, .before = animal) %>%
  
  #let's fix the order
  select(sample, 
         animal, 
         timepoint, 
         cd8_depleted, 
         infected, 
         aids, 
         sive, 
         plus_minus, 
         library_date) %>% 
  
  #next steps are to get the sample names to match
  mutate(sample = ifelse(sample <= 9, paste0("0",sample), sample)) %>%
  mutate(sample = paste0("sample_",sample))

#converting animal data into a matrix for differential expression
metadata <- as.matrix(animals)
rownames(metadata) <- metadata[, 1]
