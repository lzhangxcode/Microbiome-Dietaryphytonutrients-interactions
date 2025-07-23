# this script is used to compare ECs among different groups 
# author: lu.zhang 
# date: 2024.04.04

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 7) {
  stop("Rscript script.R cpm_file target_date target_g1 target_g2 metatable target_EC output_dir.n", call.=FALSE)
} 


library(tidyverse)

###################################

# input file 

###################################
print("reading input file... ... ")

EC_mice <- args[1] %>% read_tsv()

target_date <- args[2]

target_group <- c(args[3], args[4])

mice_metatable <-  args[5] %>% read_tsv()

overlap_ec <- args[6] %>% read_tsv()

output_dir <- args[7]

print(str_c("target date is :", target_date))
print(str_c("target groups are :", target_group))
print(str_c("output dir is:", output_dir))

###################################

# process input file 

###################################

print("processing input file ... ...")
colnames(EC_mice)[1] <- "EC"
EC_mice_unstratified <- EC_mice %>% 
  filter(., !grepl("\\|", EC))
colnames(EC_mice_unstratified) <- colnames(EC_mice_unstratified) %>% str_remove(".clean_Abundance-CPM")


# extract the date & the group 

target_samples <- tibble(sample = colnames(EC_mice_unstratified)[2:ncol(EC_mice_unstratified)]) %>% 
  left_join(., mice_metatable) %>% 
  filter(day %in% target_date) %>% 
  filter(group %in% target_group)
  

EC_mice_unstratified_sel <- EC_mice_unstratified %>% 
  dplyr::select("EC", all_of(target_samples$sample))

# produce long format 
EC_mice_unstratified_sel_long <- EC_mice_unstratified_sel %>% 
  pivot_longer(!EC) %>% 
  left_join(., target_samples, by = c("name" = "sample"))

# produce the final EC long format 
EC_mice_unstratified_sel_long_final <- EC_mice_unstratified_sel_long 

###################################

# wilcox test  

###################################
print("wilcox test ... ...")

EC_mice_unstratified_sel_long_final_statistics <- EC_mice_unstratified_sel_long_final %>%
  group_by(EC) %>%
  rstatix::wilcox_test(value ~ group) %>%
  ungroup()

EC_mice_unstratified_sel_long_final_statistics_sel <- EC_mice_unstratified_sel_long_final_statistics %>%
  mutate(q = p.adjust(p, method = "fdr")) %>%
  filter(p < 0.05)


write_tsv(EC_mice_unstratified_sel_long_final_statistics, str_c(output_dir, "all_ecs_wilcox_statistics.tsv"))

###################################

# filtered the targeted ECs 
# eg. strawberry compounds

###################################
print("filtered targeted EC... ...")
library(rstatix)
# filter the EC abundance table 
EC_mice_unstratified_sel_long_final_overlap <- EC_mice_unstratified_sel_long_final %>%
  dplyr::filter(EC %in% intersect(EC_mice_unstratified_sel_long_final_statistics_sel$EC, overlap_ec$ec))

# filter the statistics table 
EC_mice_unstratified_sel_long_final_statistics_sel_overlap <- EC_mice_unstratified_sel_long_final_statistics_sel %>%
  dplyr::filter(EC %in% intersect(EC_mice_unstratified_sel_long_final_statistics_sel$EC, overlap_ec$ec)) %>%
  add_y_position(scales = "free_y")

# write the statistic result into tsv table 
write_tsv(EC_mice_unstratified_sel_long_final_statistics_sel_overlap, str_c(output_dir, str_c(target_date, "_sig_EC_overlap.tsv"), sep = "/"))


###################################

# save the RData 

###################################

print("saving RData... ")
save.image(str_c(output_dir, str_c(target_date, "_sig_EC.RData")))


