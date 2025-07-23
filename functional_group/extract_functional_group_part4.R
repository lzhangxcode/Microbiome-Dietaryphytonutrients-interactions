# the aim of this script is to fisher test enrichment analysis 
# date: 2024.may.24 

library(tidyverse)
functional_group <- read_tsv("functional_group_result.tsv")
metabolite_group <- read_tsv("../secondary_compounds_overlap.tsv")

# make the 4 cells table 
tibble_fisher_p <- tibble(functional_group = colnames(functional_group)[-c(1,2)], p = -1, q = -1, direction = -1)
unique_group <- metabolite_group$compounds[metabolite_group$group %in% "unique"]
fisher_test_result <- list()
contingency_table <- list()
for (group in colnames(functional_group)[-c(1,2)]){
  
  print(group)
  in_unique_group <- functional_group %>% 
    select(CompoundID, all_of(group)) %>% 
    filter(CompoundID %in% unique_group) 
  
  # Count the number of compounds in unique group containing the functional group
  in_unique_group_contain_function <- sum(in_unique_group[[group]] > 0)
  print(str_c("unique group contain the function:", in_unique_group_contain_function))
  # Count the total number of compounds in the unique group
  total_in_unique_group <- nrow(in_unique_group)
  
  # Count the number of compounds not in unique group containing the functional group
  not_in_unique_group <- functional_group %>%
    filter(!CompoundID %in% unique_group) %>%
    select(CompoundID, all_of(group))
  
  not_in_unique_group_contain_function <- sum(not_in_unique_group[[group]] > 0) 
  print(str_c("overlap group contain the function:", not_in_unique_group_contain_function))
  
  # Count the total number of compounds not in unique group
  total_not_in_unique_group <- nrow(not_in_unique_group)
  
  # Create a 2x2 contingency table
  contingency_table[[group]] <- matrix(c(in_unique_group_contain_function, 
                                total_in_unique_group - in_unique_group_contain_function,
                                not_in_unique_group_contain_function,
                                total_not_in_unique_group - not_in_unique_group_contain_function), 
                              nrow = 2, byrow = TRUE)
  print(contingency_table[[group]])
  # Perform Fisher's exact test
  fisher_test_result[[group]] <- fisher.test(contingency_table[[group]])
  
  # Extract p-value
  p_value <- fisher_test_result[[group]]$p.value
  tibble_fisher_p$p[tibble_fisher_p$functional_group == group] <- p_value

  # Calculate proportion of compounds containing the functional group in each group
  prop_in_unique_group <- in_unique_group_contain_function / total_in_unique_group
  prop_in_overlap_group <- not_in_unique_group_contain_function / total_not_in_unique_group
  
  # Determine direction of enrichment
  direction <- ifelse(prop_in_unique_group > prop_in_overlap_group, "Unique Group", "Overlap Group")
  tibble_fisher_p$direction[tibble_fisher_p$functional_group == group] <- direction
  
  
}

tibble_fisher_p$q <- p.adjust(tibble_fisher_p$p, method = "fdr")

tibble_fisher_p_sig <- tibble_fisher_p %>% 
  filter(p < 0.05) %>% 
  filter(q < 0.2)
# https://rdkit.org/docs/source/rdkit.Chem.Fragments.html
output_dir <- "./minisp_70ECs/functional_group/"
write_tsv(tibble_fisher_p, str_c(output_dir, "functional_group_fisher_test.tsv"))
write_tsv(tibble_fisher_p_sig, str_c(output_dir, "functional_group_fisher_test_sig_q0.2.tsv"))

save.image(str_c(output_dir, "/extract_functional_group_part4.RData"))



