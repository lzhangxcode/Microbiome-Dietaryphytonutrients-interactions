# the aim of this script is to make figure for the mice RNA data
# author: lu.zhang 
# date: 2024.08.04 

library(tidyverse)
library(rstatix)

################################################################################

#      calculate the wilcox effect size 

################################################################################

# calculate day18 wilcox effectsize 

load("./sigec_day18/D18_sig_EC.RData")
EC_mice_unstratified_sel_long_final_statistics_effectsize <- EC_mice_unstratified_sel_long_final %>%
  group_by(EC) %>%
  rstatix::wilcox_effsize(value ~ group) %>%
  ungroup()

write_tsv(EC_mice_unstratified_sel_long_final_statistics_effectsize, "./sigec_day18/all_ecs_wilcox_statistics_effect_size.tsv")

rm(list = ls())


load("./sigec_day6/D6_sig_EC.RData")

EC_mice_unstratified_sel_long_final_statistics_effectsize <- EC_mice_unstratified_sel_long_final %>%
  group_by(EC) %>%
  rstatix::wilcox_effsize(value ~ group) %>%
  ungroup()

write_tsv(EC_mice_unstratified_sel_long_final_statistics_effectsize, "./sigec_day6/all_ecs_wilcox_statistics_effect_size.tsv")
rm(list = ls())

################################################################################

#       input file 

################################################################################

# the sig.EC between G2 and G3 at day6 

G2G3_day6_sig_EC <- read_tsv("./sigec_day6/D6_sig_EC_overlap.tsv") %>% 
  filter(p < 0.05)
G2G3_day6_sig_EC_full <- read_tsv("./sigec_day6/all_ecs_wilcox_statistics.tsv") 


# the sig.EC between G2 and G3 at day18

G2G3_day18_sig_EC <- read_tsv("./sigec_day18/D18_sig_EC_overlap.tsv") %>% 
  filter(p < 0.05)
G2G3_day18_sig_EC_full <- read_tsv("./sigec_day18/all_ecs_wilcox_statistics.tsv") 


# correlation with DAI and HE score both G2 and G3 

daiscore_metag_straw_clr_G2G3 <- read_tsv("./G2G3_correlation/EC_daiscore_cor_tibble.tsv") %>% 
  filter(P < 0.05) %>%
  mutate(R_direction = case_when(R > 0 ~ "positive", R < 0 ~ "negative")) 

daiscore_metag_straw_clr_G2G3_full <- read_tsv("./G2G3_correlation/EC_daiscore_cor_tibble.tsv") 


hescore_metag_straw_clr_G2G3 <- read_tsv("./G2G3_correlation/HE_correlation/Total_score_EC_HEscore_cor_tibble.tsv") %>% 
  filter(P < 0.05) %>%
  mutate(R_direction = case_when(R > 0 ~ "positive", R < 0 ~ "negative")) 

hescore_metag_straw_clr_G2G3_full <- read_tsv("./G2G3_correlation/HE_correlation/Total_score_EC_HEscore_cor_tibble.tsv")

# abundance data 

metag_clr <- read_tsv("/metatranscriptomics_join_enzyme_clr.tsv")
colnames(metag_clr)[1] <- "EC"
colnames(metag_clr) <- str_remove(colnames(metag_clr), ".clean_Abundance-CPM")

sample_day6 <- read_tsv("mice_metatable_G2G3.txt") %>% 
  filter(day %in% "D6")
sample_day18 <- read_tsv("mice_metatable_G2G3.txt") %>% 
  filter(day %in% "D18")

metag_day6 <- metag_clr %>% select("EC", sample_day6$sample) 
metag_day18 <- metag_clr %>% select("EC", sample_day18$sample)

# ec list : sig.EC, correlated with at least one of the clinical parameters
overlap_direciton_same_trend <- read_tsv("strawberry_summary_sig.EC_day6andday18_metatranscriptomcis_overlap_direction_only_one_groups.tsv")

# all meta 
sample_overlap <- read_tsv("./mice_metatable_G2G3.txt") 

# effect size 
day18_effectsize <- read_tsv("./sigec_day18/all_ecs_wilcox_statistics_effect_size.tsv")
day6_effectsize <- read_tsv("./sigec_day6/all_ecs_wilcox_statistics_effect_size.tsv")


################################################################################

#       extract the EC list of either sig.EC at day6 or sig.EC at day18 

################################################################################

overlap_EC <- bind_rows(G2G3_day6_sig_EC %>% mutate(day = "day6"),
                        G2G3_day18_sig_EC %>% mutate(day = "day18"))

dai_overlap <- intersect(daiscore_metag_straw_clr_G2G3$EC,  overlap_EC$EC) 
he_overlap <- intersect(hescore_metag_straw_clr_G2G3$EC, overlap_EC$EC)

ec_list <- c(dai_overlap, he_overlap) %>% unique() 

overlap_direciton_same_trend_correlated <- overlap_direciton_same_trend %>% 
  filter(dai_cor == "Y" | he_cor == "Y") 

################################################################################

#       extract the correlation value and the p value 

################################################################################

calculate_the_fc_for_each_ec <- function(input_ec = NA, input_metag_sel = NA, meta = sample_overlap){
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
  
  fc <- g3_median/g2_median 
  result_list <- list("g2_median" = g2_median, 
                      "g3_median" = g3_median,
                      "fc" = fc)
  return(result_list)
}


overlap_direciton_same_trend_correlated_fc <- overlap_direciton_same_trend_correlated %>% 
  rowwise() %>% 
  mutate(g2_median_d6 = calculate_the_fc_for_each_ec(input_ec = ECs, input_metag_sel = metag_day6) %>% .[[1]]) %>% 
  mutate(g3_median_d6 = calculate_the_fc_for_each_ec(input_ec = ECs, input_metag_sel = metag_day6) %>% .[[2]]) %>% 
  mutate(fc_d6 = calculate_the_fc_for_each_ec(input_ec = ECs, input_metag_sel = metag_day6) %>% .[[3]]) %>% 
  mutate(g2_median_d18 = calculate_the_fc_for_each_ec(input_ec = ECs, input_metag_sel = metag_day18) %>% .[[1]]) %>% 
  mutate(g3_median_d18 = calculate_the_fc_for_each_ec(input_ec = ECs, input_metag_sel = metag_day18) %>% .[[2]]) %>% 
  mutate(fc_d18 = calculate_the_fc_for_each_ec(input_ec = ECs, input_metag_sel = metag_day18) %>% .[[3]])

overlap_direciton_same_trend_correlated_fc_pvalues <- overlap_direciton_same_trend_correlated_fc %>% 
  left_join(., G2G3_day6_sig_EC_full %>% select(EC, p) %>% rename(ECs = EC, d6_p = p)) %>% 
  left_join(., G2G3_day18_sig_EC_full %>% select(EC, p) %>% rename(ECs = EC, d18_p = p)) %>% 
  left_join(., daiscore_metag_straw_clr_G2G3_full %>% select(EC, P, R) %>% rename(ECs = EC, dai_r = R, dai_p = P)) %>% 
  left_join(., hescore_metag_straw_clr_G2G3_full %>% select(EC, P, R) %>% rename(ECs = EC, he_r = R, he_p = P)) %>% 
  left_join(., day6_effectsize %>% select(EC, effsize) %>% rename(ECs = EC, day6_effsize = effsize)) %>% 
  left_join(., day18_effectsize %>% select(EC, effsize) %>% rename(ECs = EC, day18_effsize = effsize)) 

################################################################################

#       Complex Heatmap for the correlation figure 

################################################################################

library(ComplexHeatmap)
library("circlize")

ec_compounds_metag_g2g3 <- read_tsv("EC_strawberry_both_metat.tsv") %>% 
  select(EC, Compounds) %>% 
  unique()
keep_ec <- ec_compounds_metag_g2g3 %>% # this compound is a false record after manual check. So I remove it
  filter(!Compounds %in% "CDBNO:14035;CHEBI:15366;CHEBI:30089;CHEBI:32029;CHEBI:32954;CHEBI:62947;CHEBI:62984;CHEBI:63045;CHEM000606;CHEMLIST:4013384;CHEMLIST:4013767;CHEMLIST:4014426;CHEMLIST:4038148;CHEMLIST:4085560;CHEMLIST:4250360;CHEMLIST:4253661;CHEMLIST:4254351;CHEMLIST:4257788;CHEMLIST:4263892;CHEMLIST:4276491;CHEMLIST:4277629;CHEMLIST:4277651;CHEMLIST:4278482;CHEMLIST:4278541;CID:10129899;CID:10915;CID:11192;CID:12432;CID:15337;CID:175;CID:17789862;CID:206;CID:22833492;CID:2514;CID:517044;CID:65313;CID:8895;CID:9317")


overlap_direciton_same_trend_correlated_fc_pvalues_sel <- overlap_direciton_same_trend_correlated_fc_pvalues %>% 
  filter(ECs %in% unique(keep_ec$EC))

heatmap_input <- overlap_direciton_same_trend_correlated_fc_pvalues_sel %>% 
  mutate(day6_change = g3_median_d6 - g2_median_d6, day18_change = g3_median_d18 - g2_median_d18) %>% 
  dplyr::select(ECs, day6_change, day18_change, dai_r, he_r)

heatmap_input_pvalue <- overlap_direciton_same_trend_correlated_fc_pvalues_sel %>% 
  dplyr::select(ECs, d6_p, d18_p, dai_p, he_p)

# Convert the data to matrix format
mat <- as.matrix(heatmap_input[, -1])
rownames(mat) <- heatmap_input$ECs

# Define a function to map p-values to colors (gray for p < 0.05, white otherwise)
pval_to_color <- function(pval) {
  ifelse(pval < 0.05, "gray", "white")
}


# order mat and heatmap_input_pvalue
mat_order <- mat[order(rownames(mat)),]
heatmap_input_pvalue_order <- heatmap_input_pvalue[order(heatmap_input_pvalue$ECs), ]

d6_18_p <- heatmap_input_pvalue_order[, c("d6_p", "d18_p")]
dai_he_p <- heatmap_input_pvalue_order[, c("dai_p", "he_p")]

# Draw the heatmap with the row annotation
col_fun = colorRamp2(c(0, 2.5, 5), c( "white", "#2d7ecc","#195c9c"))

ht <- Heatmap(mat_order[, 1:2], name = "Change", 
              row_title = "Enzyme",
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              column_labels = c("Day6", "Day18"),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              #heatmap_legend_param = list(title = "Expression", at = c(min(mat_order), 0, max(mat_order))),
              #right_annotation = row_anno, 
              width = unit(1, "cm"),
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(d6_18_p[i, j] < 0.05)
                  grid.text("*", x, y, gp = gpar(fontsize = 10))
              }
)

#cor_col_fun = colorRamp2(c(0, -1), c( "white", "#199c8b"))
cor_col_fun = colorRamp2(c(0, -0.4, -0.7, -1), c( "white", "#cef5f0", "#88f2e4", "#199c8b"))

ht_cor <- Heatmap(mat_order[, 3:4], name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE, width = unit(1, "cm"),
                  column_labels = c("DAI score", "Histological score"),
                  column_names_gp = gpar(fontsize = 10),
                  col = cor_col_fun,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(dai_he_p[i, j] < 0.05)
                      grid.text("*", x, y, gp = gpar(fontsize = 10))
                  }
)

ht + ht_cor


P_sig = Legend(pch = "*", type = "points", labels = "P < 0.05")

pdf("metranscriptomics_RNAv2.pdf", height = 7*1.5, width = 7*1.5)
draw(ht + ht_cor, annotation_legend_list = list(P_sig), merge_legend = TRUE)
dev.off()


#############################################################################################

# correlation figure for the 6.3.1.1 

##############################################################################################
# 2024.08.06 

# Histological score correlation figure 

rm(list = ls())

load("./HE_correlation/Total_score_correlation_HEscore.RData")
library(ggpubr)
target_ec <- "6.3.1.1"
combine_tibble_for_cor_plot <- tibble(EC_Abundance = EC_mice_unstratified_d18_temp[target_ec,] %>% as.numeric() %>% as.vector(),                                                                                    DAI_score = HE_score_d18_order$HE_score %>% as.numeric()) %>% 
  mutate(sample = colnames(EC_mice_unstratified_d18_temp)) %>% 
  left_join(., Group_d18)
p_value <- EC_HEscore_cor_tibble_q$P[EC_HEscore_cor_tibble_q$EC == target_ec] %>% round(., digits=5)
r_value <- EC_HEscore_cor_tibble_q$R[EC_HEscore_cor_tibble_q$EC == target_ec] %>% round(., digits=5)
combine_tibble_for_cor_plot$Group <- combine_tibble_for_cor_plot$group
combine_tibble_for_cor_plot$Group[combine_tibble_for_cor_plot$Group == "G2"] <- "DSS+Vehicle"
combine_tibble_for_cor_plot$Group[combine_tibble_for_cor_plot$Group == "G3"] <- "DSS+Strawberry"

p_cor <- ggscatter(combine_tibble_for_cor_plot, title = target_ec, x = "EC_Abundance", y = "DAI_score", size = 4,                                             
                   add = "reg.line", conf.int = TRUE, add.params = list(color = "white", fill = "lightgray"),cor.coeff.args = list(label.y = 2.5), 
                   cor.coef = TRUE, cor.method = "spearman", cor.coef.size = 8, color = "Group", 
                   #xlab = "EC Abundance", ylab = "DAI score")
                   xlab = "", ylab = "") +
  theme(text = element_text(size=20),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle(target_ec) +
  scale_color_manual(labels = c("DSS+Strawberry", "DSS+Vehicle"), values=c("#1b5d9d", "#00695d")) + geom_point(size = 4, aes(color = Group)) + 
  xlab("Expression") + ylab("Histological score")

ggsave("./HE_correlation/total_score_6.3.1.1_RNA.pdf", p_cor, width = 5.5)



