# this script is used to draw the PCoA plot for the targeted ecs between different groups
# name : lzhang
# date : 2024.04.05 


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 7) {
  stop("Rscript script.R cpm_file target_date target_g1 target_g2 metatable overlap_EC output_dir.n", call.=FALSE)
} 

############################################

# input profile 

############################################

library(tidyverse)

# mice metatable 
mice_metatable <- args[1] %>% read_tsv()

# ec profile 
EC_mice <- args[2] %>% read_tsv()

# ec list 
ec_list <- args[3] %>% read_tsv()

# day comparison or group comparison 
target_comparison <- args[4]

# target groups
target_groups <- args[5] # G2- DSS + Vehicle, G3- DSS + strawberry

# filename 
filename <- args[6]

# output_dir 
output_dir <- args[7]


############################################

# make comparison table for target groups

############################################
print("extract target ec profile ... ... ")

target_groups_vector <- str_split(target_groups, pattern = ",") %>% unlist()

mice_metatable_sel <- mice_metatable %>% 
  filter(!!sym(target_comparison) %in% target_groups_vector)

mice_metatable_sel_unique <- mice_metatable_sel %>% 
  select(day, group) %>% select(-!!sym(target_comparison)) %>%
  distinct() %>% 
  mutate(compar1 = target_groups_vector[1], compar2 = target_groups_vector[2]) 
mice_metatable_sel_unique$separate_G <- pull(mice_metatable_sel_unique, 1)

#comparison_tibble <- tibble(!!target_comparison := target_groups_vector) # !! makes target_comparison is recognized as a variable 


############################################

# filter the target EC profile 

############################################
print("extract target ec profile ... ... ")

colnames(EC_mice)[1] <- "EC"
EC_mice_unstratified <- EC_mice %>% 
  filter(., !grepl("\\|", EC))
colnames(EC_mice_unstratified) <- colnames(EC_mice_unstratified) %>% str_remove(".clean_Abundance-CPM")

# filtered by the list 
EC_mice_unstratified_sel <- EC_mice_unstratified %>% 
  filter(EC %in% ec_list$EC)

# filtered the target profile depends on the comparisons 

extract_targeted_groups <- function(input_ec = EC_mice_unstratified_sel, target_col = NA, target = NA, meta = mice_metatable_sel){ # target_col = colnames(mice_metatable_sel_unique)[1]
  
  target_samples <- meta %>% filter(!!sym(target_col) %in% target)
  print(target_samples)
  input_ec_sel <- input_ec %>% dplyr::select(EC, all_of(target_samples$sample))
  input_ec_sel_nozero <- input_ec_sel %>% 
    column_to_rownames(., var = "EC") %>% 
    mutate(rowsum = rowSums(.)) %>% 
    filter(rowsum > 0) %>% 
    select(-rowsum) %>%
    rownames_to_column(., var = "EC") 
  
  result <- list(input_ec_sel_nozero, target_samples)
  return(result)
  
}

mice_metatable_sel_unique_nest <- mice_metatable_sel_unique %>% 
  rowwise() %>% 
  mutate(ec_profile = list(extract_targeted_groups(target_col = colnames(mice_metatable_sel_unique)[1],
                                                   target = separate_G) %>% .[[1]])) %>% 
  mutate(meta_sel = list(extract_targeted_groups(target_col = colnames(mice_metatable_sel_unique)[1],
                                                 target = separate_G) %>% .[[2]]))


############################################

# PCoA distance calculation  

############################################
print("calculate distance and statistics ... ... ")

input_ec_profile <- mice_metatable_sel_unique_nest$ec_profile[[1]]
input_meta <- mice_metatable_sel_unique_nest$meta_sel[[1]]

PCoA_distance_calculation <- function(input_ec_profile = NA, input_meta = NA, label = NA){ # label = target_comparison
  # ec profile 
  profile <- input_ec_profile %>% 
    column_to_rownames(., var = "EC") %>% 
    t()

  # distance calculation 
  set.seed(1234)
  distance.bray <- profile %>%
    vegan::vegdist(.,method = 'bray') %>%
    as.matrix() 
  
  set.seed(1234)
  aitchson_input <- profile
  
  distance.aitchison <- profile %>%
    vegan::vegdist(.,method = 'robust.aitchison') %>%
    as.matrix() 
  
  if (label == "group"){
    print(label)
    # adonis 
    permanova_result_aitchison<- vegan::adonis2(distance.aitchison ~ group, data = input_meta)
    
  }else if (label == "day"){
    print(label)
    # adonis 
    permanova_result_aitchison<- vegan::adonis2(distance.aitchison ~ day, data = input_meta)
  }

  # statistical result 
  statistic_result <- list()
  statistic_result[["aitchison"]] <- distance.aitchison
  
  statistic_result[["aitchison_adonis"]] <- permanova_result_aitchison
  
  statistic_result[["meta"]] <- input_meta
  
  return(statistic_result)
  
}

mice_metatable_sel_unique_nest_statistic <- mice_metatable_sel_unique_nest %>% 
  rowwise() %>% 
  mutate(statistic_out = list(PCoA_distance_calculation(input_ec_profile = ec_profile,
                                                        input_meta = meta_sel,
                                                        label = target_comparison))) 

############################################

# statistical table write  

############################################
print("producing rsquare table ... ... ")

extract_r_square <- function(input = NA){
  
  x <- input
  
  aitchison_r2 <- x[["aitchison_adonis"]]$R2[1]
  aitchison_p <- x[["aitchison_adonis"]]$`Pr(>F)`[1]
  
  summary_tb <- tibble(aitchison_r2 = aitchison_r2, aitchison_p = aitchison_p)

  return(summary_tb)
  
}



mice_metatable_sel_unique_nest_statistic_rsquare <- mice_metatable_sel_unique_nest_statistic %>% 
  rowwise() %>% 
  mutate(r_square_p_table = list(extract_r_square(input = statistic_out)))


############################################

# visualization the pcoa plot 

############################################
# function for calculate the points for pcoa plot 
print("drawing pcoa plot ... ... ")

PcoA_plot <- function(points_df = NA, eig = NA, p_value = NA, rsqaure = NA, tag = NA){
    p <- ggpubr::ggscatter(points_df, x = "x", y = "y",
                           color = "Group", shape = "Group", fill = "Group", ellipse = TRUE) + # )
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
      annotate(geom="text", x = max(points_df$x), y = (max(points_df$y) + 1.5*sd(points_df$y)), label = list(bquote(P == .(as.vector(p_value))))) + # x = 0.2, y = 0.3, 
      annotate(geom="text", x = max(points_df$x), y = (max(points_df$y) + sd(points_df$y)), label = list(bquote(R^2 == .(as.vector(rsqaure))))) + # x = 0.2, y = 0.3, 
      ggtitle(tag) + 
      # ggforce::geom_mark_ellipse(aes(fill = Group, color = Group)) + 
      theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  
  return(p)
}

PCoA_plots_all <- function(distance.bray = NA, permanova_result = NA, output_dir = NA, metadat_sel = NA, filetag = NA, g_or_d = NA){
  
  set.seed(1234)
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  
  metadat_sel <- metadat_sel %>% dplyr::filter(sample %in% rownames(points)) %>% 
    mutate(samples = sample) %>% 
    column_to_rownames(., var = "sample")
  metadat_sel <- metadat_sel[rownames(points),]

  points <- cbind(points, metadat_sel) 
  points <- points %>% mutate(Group = !!sym(g_or_d))
  
  group_pcoa <- PcoA_plot(points_df = points, eig = eig,
                               p_value = permanova_result$`Pr(>F)`[1],
                               rsqaure = round(permanova_result$R2[1], 3),
                               tag = filetag)
  
  
  
  system(str_c("mkdir -p ", output_dir, "/beta_diversity_plots/tmp/"))
  saveRDS(points, str_c(output_dir, "/beta_diversity_plots/tmp/pcoa_points_", filetag, ".rds"))
  saveRDS(eig, str_c(output_dir, "/beta_diversity_plots/tmp/pcoa_eig_", filetag, ".rds"))
  return(group_pcoa)    
  
}


# function for exceute function: PCoA_plots_all
produce_pcoa_plots <- function(distance_result = NA, output_dir = NA, metadat_sel = NA, category = NA, g_or_d = NA){
  library(ggpubr)
  
  pcoa_result_aitchison <- list()
  pcoa_result_aitchison <- PCoA_plots_all(distance.bray = distance_result[["aitchison"]], 
                                          permanova_result = distance_result[["aitchison_adonis"]], 
                                          output_dir = output_dir, 
                                          metadat_sel = metadat_sel,
                                          filetag = paste(category, "aitchison", sep = "_"),
                                          g_or_d = g_or_d)
  pcoa_result_all <- list()
  pcoa_result_all[["aitchison"]] <- pcoa_result_aitchison
  
  
  return(pcoa_result_all)
  
}



mice_metatable_sel_unique_nest_statistic_rsquare_plots <- mice_metatable_sel_unique_nest_statistic_rsquare %>% 
  rowwise() %>% 
  mutate(plots_out = list(produce_pcoa_plots(distance_result = statistic_out,
                                             output_dir = output_dir,
                                             metadat_sel = mice_metatable_sel,
                                             category = filename,
                                             g_or_d = target_comparison)))

saveRDS(mice_metatable_sel_unique_nest_statistic_rsquare_plots, str_c(output_dir, filename, "_mice_metatable_rsquare_plots.rds"))

write_tsv(mice_metatable_sel_unique_nest_statistic_rsquare_plots %>% dplyr::select(separate_G, compar1, compar2, r_square_p_table) %>% 
            unnest(cols = c(r_square_p_table)), 
          str_c(output_dir, filename, "_mice_adonis_rsquare_table.tsv"))

for (k in 1:length(mice_metatable_sel_unique_nest_statistic_rsquare_plots$separate_G)){
  ggsave(str_c(output_dir,str_c(mice_metatable_sel_unique_nest_statistic_rsquare_plots$separate_G[[k]], "_PCoA_plots_aitchison.pdf")),
         mice_metatable_sel_unique_nest_statistic_rsquare_plots$plots_out[[k]]$aitchison)
}


############################################

# save RData  

############################################

print("saving RData... ...")
save.image(str_c(output_dir, filename, "_PCoA_", target_comparison, ".RData"))


