# this script is to correlation with HE score : similar code as DAI score as in script correlation_target_ec.R 
# author: lu.zhang 
# date: 2024.04.09 

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 7) {
  stop("Rscript script.R abudance_file target_ec mice_metatable hescore hescore_type output_dir group_label.n", call.=FALSE)
} 


# 
library(tidyverse)
###################################

# input file 

###################################

# ec profile 
EC_mice <- args[1] %>% read_tsv()

# ec list 
ec_list <- args[2] %>% read_tsv()

# mice metatable 
mice_metatable <- args[3] %>% read_tsv()

# HE score 
HE_score <- readxl::read_xlsx(args[4], sheet = "Score") %>% 
  fill(., Group) %>% 
  filter(!`Animal No.` %in% c("SEM", "Mean"))

# score 
#score <- "Total,score" # "IBD severity level" "Hyperplasia" "Ulcers" "Lesion area" "Total score"
score <- args[5]
score <- str_replace(score, ",", " ")

# output_dir 
output_dir <- args[6]

# label 
label <- args[7]
label <- unlist(strsplit(label, ",", " "))

###################################

# extract the target ECs  

###################################

print("extract target ec profile ... ... ")

colnames(EC_mice)[1] <- "EC"
EC_mice_unstratified <- EC_mice %>% 
  filter(., !grepl("\\|", EC))
colnames(EC_mice_unstratified) <- colnames(EC_mice_unstratified) %>% str_remove(".clean_Abundance-CPM")

# filtered by the list 
EC_mice_unstratified_sel <- EC_mice_unstratified %>% 
  filter(EC %in% ec_list$EC) %>% 
  select(EC, mice_metatable$sample)

###################################

# proceed HE score 

###################################
print("proceed HE score... ...")
HE_score_sel <- HE_score %>% select("Group", "Animal No.", all_of(score))
colnames(HE_score_sel)[2] <- c("AnimalID")

HE_score_sel_d18 <- HE_score_sel %>% 
  mutate(sample = str_c("S", AnimalID))


###################################

# spearman correlation table 

###################################
print("build spearman correlation table... ...")

Group_d18 <- mice_metatable %>% filter(group %in% label & day == "D18")

HE_score_d18_order <- HE_score_sel_d18 %>% 
  left_join(Group_d18,., by = "sample") 
HE_score_d18_order$HE_score <- as.numeric(HE_score_d18_order$`Total score`)


###################################

# further filter the ec profile 

###################################

print("ec profile further filtered only G4 & day18")
EC_mice_unstratified_d18 <- EC_mice_unstratified_sel %>% 
  select(EC, HE_score_d18_order$sample)

print("the EC table dim")
print(dim(EC_mice_unstratified_d18))



###################################

# correlation  

###################################

print("calculate the correlation...")


EC_HEscore_cor_list <- list()
EC_HEscore_cor_tibble <- tibble(EC = EC_mice_unstratified_d18$EC, P = -1, R = -100)
EC_mice_unstratified_d18_temp <- EC_mice_unstratified_d18 %>% 
  column_to_rownames(., var = "EC")
for (i in EC_mice_unstratified_d18$EC){
  print(i)
  ec_x <- EC_mice_unstratified_d18_temp[i,] %>% as.numeric() %>% as.vector()
  dai_y <- HE_score_d18_order$HE_score %>% as.numeric()
  EC_HEscore_cor_list[[i]] <- cor.test(x = ec_x, y = dai_y, method = 'spearman')
  EC_HEscore_cor_tibble[EC_HEscore_cor_tibble$EC %in% i,]$P <- EC_HEscore_cor_list[[i]]$p.value
  EC_HEscore_cor_tibble[EC_HEscore_cor_tibble$EC %in% i,]$R <- EC_HEscore_cor_list[[i]]$estimate %>% as.vector()
  
}

EC_HEscore_cor_tibble$Q <- p.adjust(EC_HEscore_cor_tibble$P, method = "fdr")

score <- str_replace(score, " ", "_")
write_tsv(EC_HEscore_cor_tibble, str_c(output_dir, str_c(score, "_EC_HEscore_cor_tibble.tsv")))


print("RData saving ... ...")

save.image(str_c(output_dir, str_c(score, "_correlation_HEscore.RData")))


