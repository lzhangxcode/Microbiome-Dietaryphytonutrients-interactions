# summary of the sig.ECs that are sig.diff between strawberry and non-strawberry group at day6 or day18 in RNA level 
# author: lu.zhang
# date: 2024.july.31 

################################################################################

#  load the libraries 

################################################################################

library(tidyverse)

################################################################################

#  load the input files 

################################################################################

# the sig.EC between G2 and G3 at day6 

G2G3_day6_sig_EC <- read_tsv("./sigec_day6/D6_sig_EC_overlap.tsv") %>% 
  filter(p < 0.05)


# the sig.EC between G2 and G3 at day18

G2G3_day18_sig_EC <- read_tsv("./sigec_day18/D18_sig_EC_overlap.tsv") %>% 
  filter(p < 0.05)


# correlation with DAI and HE score both G2 and G3 

daiscore_metag_straw_clr_G2G3 <- read_tsv("./G2G3_correlation/EC_daiscore_cor_tibble.tsv") %>%  # variable called metag,it's a typo，input from rna
  filter(P < 0.05) %>%
  mutate(R_direction = case_when(R > 0 ~ "positive", R < 0 ~ "negative")) 


hescore_metag_straw_clr_G2G3 <- read_tsv("./G2G3_correlation/HE_correlation/Total_score_EC_HEscore_cor_tibble.tsv") %>% # variable called metag,it's a typo，input from rna
  filter(P < 0.05) %>%
  mutate(R_direction = case_when(R > 0 ~ "positive", R < 0 ~ "negative")) 

################################################################################

#  calculate which are both sig.diff 

################################################################################

overlap_EC <- bind_rows(G2G3_day6_sig_EC %>% mutate(day = "day6"),
                        G2G3_day18_sig_EC %>% mutate(day = "day18"))

overlap_EC_list <- overlap_EC %>% 
  group_by(EC) %>% 
  summarise(n = n()) %>% 
  filter(n >= 1)

overlap_EC_sel <- overlap_EC %>% 
  filter(EC %in% overlap_EC_list$EC)

intersect(daiscore_metag_straw_clr_G2G3$EC, overlap_EC_sel$EC) %>% length() 
intersect(hescore_metag_straw_clr_G2G3$EC, overlap_EC_sel$EC) %>% length() 


################################################################################

#  write the summary output 

################################################################################

# check the directions of the abundance 

dai_overlap <- intersect(daiscore_metag_straw_clr_G2G3$EC,  overlap_EC_sel$EC) 
he_overlap <- intersect(hescore_metag_straw_clr_G2G3$EC, overlap_EC_sel$EC)

# input data 
metag_clr <- read_tsv("metatranscriptomics_join_enzyme_clr.tsv")
colnames(metag_clr)[1] <- "EC"
colnames(metag_clr) <- str_remove(colnames(metag_clr), ".clean_Abundance-CPM")

# select day6 and day18 data separately 
sample_overlap <- read_tsv("mice_metatable_G2G3.txt") 
sample_day6 <- read_tsv("mice_metatable_G2G3.txt") %>% 
  filter(day %in% "D6")
sample_day18 <- read_tsv("mice_metatable_G2G3.txt") %>% 
  filter(day %in% "D18")

metag_day6 <- metag_clr %>% select("EC", sample_day6$sample) 
metag_day18 <- metag_clr %>% select("EC", sample_day18$sample)

input_ec <- dai_overlap[1]
input_metag_sel <- metag_day6

check_which_one_is_higher <- function(input_ec = NA, input_metag_sel = NA, meta = sample_overlap){
  input_metag_sel_long <- input_metag_sel %>% 
    pivot_longer(!EC) %>% 
    filter(EC %in% input_ec) %>% 
    left_join(., meta, by = c("name" = "sample"))
  
  print(input_ec)
  g2_median <- input_metag_sel_long %>% 
    filter(group == "G2") %>% 
    pull(value) %>% 
    median()
  g3_median <- input_metag_sel_long %>% 
    filter(group == "G3") %>% 
    pull(value) %>% 
    median()
  
  if (g2_median > g3_median){
    higher_value <- "g2"
  }else if (g2_median < g3_median){
    higher_value <- "g3"
  }else if (g2_median == g3_median){
    higher_value <- "unknown"
  }
  
  return(higher_value)
}



overlap_direction <- tibble(ECs = overlap_EC_sel$EC %>% unique()) %>% 
  mutate(dai_cor = ifelse(ECs %in% dai_overlap, "Y", "N")) %>% 
  mutate(he_cor = ifelse(ECs %in% he_overlap, "Y", "N")) %>% 
  rowwise() %>% 
  mutate(higher_value_d6 = check_which_one_is_higher(input_ec = ECs, input_metag_sel = metag_day6)) %>% 
  mutate(higher_value_d18 = check_which_one_is_higher(input_ec = ECs, input_metag_sel = metag_day18))


# same trend 
overlap_direciton_same_trend <- overlap_direction %>% filter(higher_value_d6 == "g3") %>% filter(higher_value_d18 == "g3")  

save.image("strawberry_summary_sig.EC_day6andday18_metatranscriptomics_only_one_group.RData")

write_tsv(overlap_direciton_same_trend, "./strawberry_summary_sig.EC_day6andday18_metatranscriptomcis_overlap_direction_only_one_groups.tsv")


