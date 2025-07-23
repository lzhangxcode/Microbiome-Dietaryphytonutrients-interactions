# the aim of the script is to get the pvalue of adjusted age, gender, BMI 
# author: lu.zhang
# date: 2024.juli.29 

#####################################################################

#   input files 

#####################################################################

library(tidyverse)
library(vegan)

modified_ec_profile_beta_input <- readRDS(str_c("./adonis/", "modified_ec_profile_beta_input.rds"))


####################################################################

#   calculate the pvalue 

#####################################################################

# it has been calculated before from function : distance_produce 
load("./adonis/step16_manual_check_adonis_env_compoundclass.RData")
> modified_ec_profile_distance_result[["bray_adonis"]]
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = distance.bray ~ age + gender + BMI + continent, data = metadat_sel)
            Df SumOfSqs      R2        F Pr(>F)    
continent    3    8.414 0.13876 106.9326  0.001 ***
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> modified_ec_profile_distance_result[["aitchison_adonis"]]
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = distance.aitchison ~ age + gender + BMI + continent, data = metadat_sel)
            Df SumOfSqs      R2       F Pr(>F)    
continent    3  1647997 0.11770 88.6464  0.001 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# don't need to save RData here 
