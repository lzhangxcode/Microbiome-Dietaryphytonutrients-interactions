# this script is used for the correlation between species abundance/prevalence and compounds number 
# author: Lu.Zhang
# date:2024.04.26 

setwd("./correlation_species_prevalenceandabundance/")


input_dir <- "./correlation_species_prevalenceandabundance/"
output_dir <- "./correlation_species_prevalenceandabundance/result/"

library(tidyverse)
#####################################################

# prepare the input file 

#####################################################

# ec raw profile 
ec_raw_profile <- str_c(input_dir,"merge_Cohort.tsv") %>% read_tsv()
colnames(ec_raw_profile)[1] <- "ec"
colnames(ec_raw_profile) <- str_remove_all(colnames(ec_raw_profile), "_Abundance-CPM") %>% str_remove_all(., "_T1") %>% str_remove_all("_Abundance-RPKs")

# modified ec profile - np
modified_ec_profile_sep <- readRDS(str_c(input_dir, "modified_ec_profile_after_np.rds"))
# metatable 
combine_metatable <- read_tsv(str_c(input_dir, "combine_metadata_sel.tsv"), guess_max = 3000)
# species table 
metaphlan_file <- read_tsv(str_c(input_dir, "merge_Cohort_metaphlan_species.tsv"))
colnames(metaphlan_file) <- str_remove_all(colnames(metaphlan_file), "_metaphlan_bowtie2") %>% str_remove_all(., "_T1")
metaphlan_file_selv2 <- metaphlan_file %>%
  dplyr::select(., c("clade_name", combine_metatable$internal_sample_id))
# tree tip 
tree_tip <- read_tsv(str_c(input_dir, "tree_tips.tsv"))


#####################################################

# calculate the species-compounds link 

#####################################################

ec_compounds_link <- modified_ec_profile_sep %>% select(ec_simplied, Compounds) %>% unique()
 
ec_raw_profile_sel <- ec_raw_profile %>%
  select(., c("ec", combine_metatable$internal_sample_id)) %>% #3068samples 
  add_column(ec_simplied = .$ec %>% str_split(., "\\|", simplify = T) %>% .[,1]) 

ec_raw_profile_combine_sel <- ec_raw_profile_sel %>% 
  left_join(., ec_compounds_link) %>% 
  select(ec, ec_simplied, Compounds) %>% 
  unique() %>% 
  filter(!is.na(Compounds)) 

# remove all zero rows 
keep_items <- ec_raw_profile_sel %>% 
  select(-ec_simplied) %>% 
  column_to_rownames(., var = "ec") %>% 
  mutate(rowsum = rowSums(.)) %>% 
  filter(rowsum > 0)

# filtered out those without any linked compounds 
ec_raw_profile_combine_sel_compounds <- ec_raw_profile_combine_sel %>% 
  filter(ec %in% rownames(keep_items)) 

# linked the compounds with species name 
ec_raw_profile_combine_sel_compounds_sp <- ec_raw_profile_combine_sel_compounds %>% 
  add_column(species = str_split(.$ec, pattern = "\\|", simplify = TRUE)[,2]) %>%
  filter(species != "")

# only species and compounds and ec link 
ec_sp_compounds <- ec_raw_profile_combine_sel_compounds_sp %>% 
  select(ec_simplied, species, Compounds) %>% # unclassified
  unique()

# calcualte each specie's related compounds number 
sp_compounds_num <- ec_sp_compounds %>% 
  select(species, Compounds) %>% 
  unique() %>% 
  group_by(species) %>% 
  summarise(n = n(), compounds = paste(Compounds, collapse = "and")) %>% 
  ungroup()

#####################################################

# calculate the species abundance/prevalence 

#####################################################

# species mean abundance 
sp_abundance <- metaphlan_file_selv2 %>% 
  pivot_longer(!clade_name) %>% 
  group_by(clade_name) %>% 
  summarise(mean_sp = mean(value)) %>% 
  filter(mean_sp > 0) %>% 
  ungroup()
sp_abundance_each_sample <- metaphlan_file_selv2 %>% 
  pivot_longer(!clade_name) %>% 
  group_by(name) %>% 
  summarise(sum_sp = sum(value)) %>% 
  ungroup() # 

# species mean prevalence 
sp_prevalence <- metaphlan_file_selv2 %>% 
  pivot_longer(!clade_name) %>% 
  filter(value > 0) %>% 
  group_by(clade_name) %>% 
  summarise(sample_n = n()) %>% 
  ungroup() %>% 
  mutate(prevalence_sp = sample_n/3068)

# filter only the bacteria ones 
sp_abundance_prevalence <- sp_abundance %>% 
  left_join(., sp_prevalence) %>% 
  select(clade_name, mean_sp, prevalence_sp) %>% 
  mutate(species = paste0(str_split(clade_name, pattern = "\\|", simplify = T) %>% .[,6], ".", 
                               str_split(clade_name, pattern = "\\|", simplify = T) %>% .[,7]))


# matched the name 
sp_abundance_prevalence_sel <- sp_abundance_prevalence %>% 
  left_join(sp_compounds_num, .) %>% 
  filter(!is.na(clade_name)) 

#####################################################

# correlation analysis - spearman correlation

#####################################################
# further filter the bacteria species only 

sp_abundance_prevalence_sel_bac <- sp_abundance_prevalence_sel %>% 
  filter(grepl("k__Bacteria", clade_name)) 

corr_abundance <- sp_abundance_prevalence_sel_bac %>% rstatix::cor_test(n, mean_sp, method = "spearman")
#var1  var2      cor statistic        p method  
#<chr> <chr>   <dbl>     <dbl>    <dbl> <chr>   
#  1 n     mean_sp  0.36 53126005. 2.06e-25 Spearman

corr_prevalence <- sp_abundance_prevalence_sel_bac %>% rstatix::cor_test(n, prevalence_sp, method = "spearman")
#var1  var2            cor statistic        p method  
#<chr> <chr>         <dbl>     <dbl>    <dbl> <chr>   
#  1 n     prevalence_sp  0.28 59855149. 1.99e-15 Spearman


sp_abundance_prevalence_sel_bac$log_n <- log2(sp_abundance_prevalence_sel_bac$n)

#####################################################

# correlation analysis - plot 

#####################################################


library(ggplot2)

p_abundance <- ggplot(sp_abundance_prevalence_sel, aes(x = n, y = mean_sp)) +
  geom_point(size = 2, shape = 1) +
  geom_smooth(method=lm)
  
p_prevalence <- ggplot(sp_abundance_prevalence_sel, aes(x = n, y = mean_sp)) +
  geom_point(size = 2, shape = 1) 

ggpubr::ggscatter(sp_abundance_prevalence_sel, x = "n", y = "mean_sp",
          color = "black", shape = 1, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


sp_abundance_prevalence_sel_bac$log_ab <- log2(sp_abundance_prevalence_sel_bac$mean_sp)
abundance_compounds_plot <- ggpubr::ggscatter(sp_abundance_prevalence_sel_bac, x = "log_ab", y = "n",
                                               color = "white", fill = "#5ab4ac", shape = 21, size = 6, # Points color, shape and size
                                               add = "reg.line",  # Add regressin line
                                               add.params = list(color = "white", fill = "lightgray"), # Customize reg. line
                                               conf.int = TRUE, # Add confidence interval
                                               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                                               cor.coeff.args = list(method = "spearman", label.x = 0.1, label.sep = "\n")) + 
  xlab("log2(Abundance of bacterial species)") + ylab("N of phytonutrients")

ggsave(str_c(output_dir, "/abundance_compounds_correlation.pdf"), abundance_compounds_plot, width = 5.5, height = 5)


prevalence_compounds_plot <- ggpubr::ggscatter(sp_abundance_prevalence_sel_bac, x = "prevalence_sp", y = "n",
                  color = "white", fill = "#d8b365", shape = 21, size = 6, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "white", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.x = 0.8, label.sep = "\n")) + 
                  xlab("Prevalence of bacterial species") + ylab("N of phytonutrients")

ggsave(str_c(output_dir, "/prevalence_compounds_correlation.pdf"), prevalence_compounds_plot, width = 5.5, height = 5)


abundance_compounds_plotv2 <- abundance_compounds_plot + xlab("log2(Abundance of bacterial species)") + ylab("N of phytonutrients") 

ggsave(str_c(output_dir, "/abundance_compounds_correlationv2.pdf"), abundance_compounds_plotv2, width = 5.5, height = 5)

prevalence_compounds_plotv2 <- prevalence_compounds_plot + xlab("Prevalence of bacterial species") + ylab("N of phytonutrients")

ggsave(str_c(output_dir, "/prevalence_compounds_correlationv2.pdf"), prevalence_compounds_plotv2, width = 5.5, height = 5)


# change again the figure 

abundance_compounds_plotv3 <- abundance_compounds_plotv2 + 
  annotate("text", x = -10, y = 550, label = as.character(expression(paste(italic("R ="), " 0.36, ", italic("P ="), " 2.06e-25"))), parse = TRUE, size = 6) + 
  theme(axis.title=element_text(size = 20)) + 
  xlab("log2(abundance of bacterial species)") 

ggsave(str_c(output_dir, "/abundance_compounds_correlationv3.pdf"), abundance_compounds_plotv3, width = 5.5, height = 5)

prevalence_compounds_plotv3 <- prevalence_compounds_plotv2 + 
  annotate("text", x = 0.75, y = 550, label = as.character(expression(paste(italic("R ="), " 0.28, ", italic("P ="), " 1.99e-15"))), parse = TRUE, size = 6) +
  theme(axis.title=element_text(size = 20)) 

ggsave(str_c(output_dir, "/prevalence_compounds_correlationv3.pdf"), prevalence_compounds_plotv3, width = 5.5, height = 5)

#####################################################

# further analysis on statistiscs of compounds in each phylum

#####################################################


sp_abundance_prevalence_sel_bac_phylum <- sp_abundance_prevalence_sel_bac %>%
  mutate(phylum = str_split(clade_name, pattern = "\\|", simplify = T) %>% .[,2])


sp_abundance_prevalence_sel_bac_phylum_above200 <- sp_abundance_prevalence_sel_bac_phylum %>%
  select(species, n, compounds, phylum) %>%
  filter(n > 200) %>%
  group_by(phylum) %>%
  summarise(target_n = n()) %>%
  ungroup() %>%
  left_join(., phylum_number_nofilter) %>%
  mutate(percent = target_n/total_n) %>%
  arrange(desc(percent))



write_tsv(sp_abundance_prevalence_sel_bac_phylum_above200, str_c(output_dir, "sp_abundance_prevalence_sel_bac_phylum_above200.tsv"))


