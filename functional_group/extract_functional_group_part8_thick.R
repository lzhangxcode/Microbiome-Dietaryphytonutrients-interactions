# the aim of this script is to extract which one is unique and which one is overalp compounds : same as part6, but here rename the thick ones 
# date: 2024.07.23 
# author: lu zhang

library(tidyverse)
setwd("./minisp_70ECs/functional_group")
filtered_functional_group <- read_tsv("filtered_compounds.tsv")
np_database_overall_fullEC <- read_tsv("nutrichem_md_manual_with_non_microbiotaEC_fullEC.tsv")
system("mkdir -p images_rename_thick/")

# the figures are all in images_rename

metabolite_group <- read_tsv("../secondary_compounds_overlap.tsv")
np_database_overall_fullEC_sel <- np_database_overall_fullEC %>% 
  filter(CompoundID %in% metabolite_group$compounds) %>% 
  select(CompoundID, CompoundTag) %>% 
  unique() %>% 
  group_by(CompoundID) %>% 
  summarise(CompoundTag = paste(CompoundTag, collapse = ","))
write_tsv(np_database_overall_fullEC_sel, str_c("./minisp_70ECs/functional_group/", "np_database_overall_fullEC_sel_thick.tsv"))

np_database_top_tag <- np_database_overall_fullEC %>% 
  filter(CompoundID %in% metabolite_group$compounds) %>% 
  select(CompoundID, CompoundTag) %>%
  group_by(CompoundID) %>%
  count(CompoundTag, sort = TRUE) %>%
  slice_max(n, n = 1) %>%
  ungroup() %>%
  select(CompoundID, CompoundTag) %>% 
  unique() %>% 
  mutate(CompoundTag_new = str_replace_all(CompoundTag, " ", "_"))


filtered_functional_group_mark <- filtered_functional_group %>% 
  left_join(., metabolite_group, by = c("CompoundID" = "compounds")) %>%
  left_join(., np_database_top_tag) %>% 
  mutate(Index_new = str_c(Index, CompoundTag_new, group, sep = "__")) 

files <- list.files(path = "images_thick", pattern = "png")
i <- 1 
for (file in files){
  print(i)
  print(file)
  index <- str_split(file, pattern = "_fr")[[1]][1] 
  print(index)
  index_new <- filtered_functional_group_mark %>% 
    filter(Index == index) %>% 
    pull(Index_new)
  print(index_new)
  file_new <- str_c(index_new, file, sep = "_")
  file.copy(str_c("images_thick/", file), str_c("images_rename_thick/", file_new))
  print(file_new)
  i <- i + 1 
}

output_dir <- "./minisp_70ECs/functional_group/"
save.image(str_c(output_dir, "/extract_functional_group_part8.RData"))
