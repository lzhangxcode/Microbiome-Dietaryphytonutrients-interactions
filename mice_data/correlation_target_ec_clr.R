# this script is used to correlate the target(strawberry) ECs with the DAI score 
# Lu.Zhang 
# 2024.04.08 

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
  stop("Rscript script.R abundance_file target_ec mice_metatable daiscore output_dir.n", call.=FALSE)
} 


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

# DAI score 
DAI_score <- readxl::read_xlsx(args[4], sheet = "DAI_simple") %>% 
  janitor::row_to_names(., row_number = 1)

# output_dir 
output_dir <- args[5]

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

# proceed DAI score 

###################################
print("proceed DAI score... ...")
colnames(DAI_score)[1:2] <- c("Group", "AnimalID")

days_dai <- colnames(DAI_score)[c(3,8,13,18,23,28,33,38,43,48,53,58)] %>% str_remove("Days post grouping: ")
DAI_score_simple <- DAI_score %>% 
  .[,c(1,2,6,11,16,21,26,31,36,41,46,51,56,61)]
colnames(DAI_score_simple)[3:14] <- str_c(days_dai, "DAIscore", sep = "_")

DAI_score_final <- DAI_score_simple %>% 
  slice(-1) %>%
  tidyr::fill(., Group) 

DAI_score_d18 <- DAI_score_final %>% 
  select(AnimalID, D18_DAIscore)%>% 
  mutate(sample = str_c("S", AnimalID))


###################################

# spearman correlation table 

###################################
print("build spearman correlation table... ...")

Group_d18 <- mice_metatable %>% filter(day == "D18") 

DAI_score_d18_order <- DAI_score_d18 %>% 
  left_join(Group_d18,., by = "sample")


###################################

# further filter the ec profile 

###################################

print("ec profile further filtered day18")
EC_mice_unstratified_d18 <- EC_mice_unstratified_sel %>% 
  select(EC, DAI_score_d18_order$sample)

print("the EC table dim")
print(dim(EC_mice_unstratified_d18))



###################################

# correlation  

###################################

print("calculate the correlation...")


EC_daiscore_cor_list <- list()
EC_daiscore_cor_tibble <- tibble(EC = EC_mice_unstratified_d18$EC, P = -1, R = -100)
EC_mice_unstratified_d18_temp <- EC_mice_unstratified_d18 %>% 
  column_to_rownames(., var = "EC")
for (i in EC_mice_unstratified_d18$EC){
  print(i)
  ec_x <- EC_mice_unstratified_d18_temp[i,] %>% as.numeric() %>% as.vector()
  dai_y <- DAI_score_d18_order$D18_DAIscore %>% as.numeric()
  EC_daiscore_cor_list[[i]] <- cor.test(x = ec_x, y = dai_y, method = 'spearman')
  EC_daiscore_cor_tibble[EC_daiscore_cor_tibble$EC %in% i,]$P <- EC_daiscore_cor_list[[i]]$p.value
  EC_daiscore_cor_tibble[EC_daiscore_cor_tibble$EC %in% i,]$R <- EC_daiscore_cor_list[[i]]$estimate %>% as.vector()
  
}

EC_daiscore_cor_tibble$Q <- p.adjust(EC_daiscore_cor_tibble$P, method = "fdr")


write_tsv(EC_daiscore_cor_tibble, str_c(output_dir, "EC_daiscore_cor_tibble.tsv"))



print("RData saving ... ...")

save.image(str_c(output_dir, "correlation_strawberry.RData"))

