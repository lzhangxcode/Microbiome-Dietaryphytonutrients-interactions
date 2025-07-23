# this script is aim to draw sample accumulation curve based on random sampling result 
# author: l.zhang 
# date: 29,jan,2024 - 31,jan,2024 

########################################

#       packages needed                # 

########################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 8) {
	  stop("Rscript sample_accumulation_curve_plot.R metatable random_sampling_result_dietary_ec modified_profile all_ec_sep speceis_sep(these '_sep' files are from random sampling scripts) output_dir output_plot_dir depth_file", call.=FALSE)
} 


library(tidyverse)
library(ggplot2)

########################################


#       input                          # 


########################################

# combined metatable 
combine_metatable <- args[1] %>% read_tsv(., guess_max = 3000)

# modified ec profile : random sampling result
modified_ec <- readRDS(args[2])

# modified ec profile - per ec per line
modified_ec_profile_simple <- args[3] %>% readRDS()

# all_ec_sep 
all_ec_sep <- args[4] %>% readRDS() 

# species_sep 
species_sep <- args[5] %>% readRDS()

# output dir 
output_dir <- args[6]

# output_dir for plot 
output_dir_plot <- args[7] 

# sequencing depth file 
sequencing_depth_file <- args[8] %>% readRDS()

########################################################################################################################

#                                   save matrix                                                                        #

########################################################################################################################

save_matrix <- function(input_matrix = NA){

        enzyme_melt_Af <- reshape2::melt(input_matrix[["Africa"]]) # Var1: repeat time, Var2 : sampling_size
        enzyme_melt_Af$Var1 <- "Africa"

        enzyme_melt_Am <- reshape2::melt(input_matrix[["America"]])
        enzyme_melt_Am$Var1 <- "America"

        enzyme_melt_As <- reshape2::melt(input_matrix[["Asia"]])
        enzyme_melt_As$Var1 <- "Asia"

        enzyme_melt_Eu <- reshape2::melt(input_matrix[["Europe"]])
        enzyme_melt_Eu$Var1 <- "Europe"

        enzyme_melt_Oc <- reshape2::melt(input_matrix[["Oceania"]])
        enzyme_melt_Oc$Var1 <- "Oceania"


        all_continents_melt <- rbind(enzyme_melt_Af, enzyme_melt_Am, enzyme_melt_As, enzyme_melt_Eu, enzyme_melt_Oc)
        all_continents_melt$Var2 <- as.factor(all_continents_melt$Var2)

        return(all_continents_melt)
}


# all ec 
all_ec_matrix <- readRDS(paste(output_dir, "/enzyme_matrix_for_plot_allec.rds", sep = "/"))

# species 
sp <- readRDS(paste(output_dir, "/enzyme_matrix_for_plot_species.rds", sep = "/")) 

sp_melt <- save_matrix(input_matrix = sp)
all_continents_melt_ec <- save_matrix(input_matrix = all_ec_matrix)



########################################################################################################################

#                                   draw plots                                                                         #

########################################################################################################################

all_continents_melt_ec$feature_group <- "EC"
modified_ec$feature_group <- "Dietary EC"
sp_melt$feature_group <- "Species"

combined_input <- rbind(all_continents_melt_ec, modified_ec, sp_melt)

saveRDS(sp_melt, str_c(output_dir, "/sp_melt.rds"))
saveRDS(combined_input, str_c(output_dir, "/combined_input.rds"))
saveRDS(all_continents_melt_ec, str_c(output_dir, "/all_ec_melt.rds"))
saveRDS(modified_ec, str_c(output_dir, "/modified_ec_melt.rds"))

draw_boxplot <- function(dataset = NA, color = NA, continent = NA){

  max_table <- dataset %>% group_by(feature_group) %>% top_n(value, n = 1) %>% select(value, feature_group) %>% unique()
  modified_ec_max <- ggplot(data = dataset, aes(x = Var2, y = value, color = Var1, fill = feature_group)) +
    geom_boxplot(lwd = 0.1, outlier.size = 0.2, outlier.shape = NA)+
    labs(x = "Number of Samples", y = "Number of Features") +
    xlim(unique(dataset$Var2)) +
    ylim(1,3000) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 26, color = "black"),
          axis.text.x = element_text(size = 24, color = rep(c(rep("transparent", 49), "black"), 29)),
          axis.text.y = element_text(size = 24, color = "black"),
          axis.title.y = element_text(size = 26, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22),
          legend.key.size = unit(1, "cm")
    ) +
    scale_fill_manual(name = "", values = color, limits = continent) +
    scale_color_manual(name = "", values = color, limits = continent) +
    theme(legend.position = "none") + 
    annotate("text", x = c(length(unique(dataset$Var2)) - 50), y = max_table$value[max_table$feature_group == "EC"], label = "EC", size = 8) + 
    annotate("text", x = c(length(unique(dataset$Var2)) - 50), y = max_table$value[max_table$feature_group == "Dietary EC"], label = "Dietary EC", size = 8) + 
    annotate("text", x = c(length(unique(dataset$Var2)) - 50), y = max_table$value[max_table$feature_group == "Species"], label = "Species", size = 8)
 
   return(modified_ec_max)

}

africa <- filter(combined_input, Var1 == "Africa")
america <- filter(combined_input, Var1 == "America")
asia <- filter(combined_input, Var1 == "Asia")
europe <- filter(combined_input, Var1 == "Europe")
oceania <- filter(combined_input, Var1 == "Oceania")

        
p_draw_africa <- draw_boxplot(africa, color = "#E64B35FF", continent = "Africa")
p_draw_america <- draw_boxplot(america, color = "#4DBBD5FF", continent = "America")
p_draw_asia <- draw_boxplot(asia, "#00A087FF", "Asia")
p_draw_europe <- draw_boxplot(europe, "#3C5488FF", "Europe")
p_draw_oceania <- draw_boxplot(oceania, "#F39B7FFF", "Oceania")


pdf(paste(output_dir_plot, "/africa_accumulation_curves.pdf", sep = "/"), height = 8, width = 25)
p_draw_africa
dev.off()

pdf(paste(output_dir_plot, "/america_accumulation_curves.pdf", sep = "/"), height = 8, width = 25)
p_draw_america
dev.off()

pdf(paste(output_dir_plot, "/asia_accumulation_curves.pdf", sep = "/"), height = 8, width = 25)
p_draw_asia
dev.off()

pdf(paste(output_dir_plot, "/europe_accumulation_curves.pdf", sep = "/"), height = 8, width = 25)
p_draw_europe
dev.off()

pdf(paste(output_dir_plot, "/oceania_accumulation_curves.pdf", sep = "/"), height = 8, width = 25)
p_draw_oceania
dev.off()


pdf(paste(output_dir_plot, "continents_accumulation_curves.pdf", sep = "/"), height = 30, width = 20)
gridExtra::grid.arrange(p_draw_africa, p_draw_america, p_draw_asia, p_draw_europe, p_draw_oceania, ncol = 1)
dev.off()

# for each continent draw boxplots 

p_draw <- draw_boxplot(all_continents_melt_ec, 
		       color = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF"), 
		       continent = c("Africa", "America", "Asia", "Europe", "Oceania"))
pdf(paste(output_dir_plot, "continents_combine_modified_enzyme_numbers.pdf", sep = "/"), height = 20, width = 40)
p_draw
dev.off()


########################################################################################


#                       calculate metatable age,gender,bmi median                      # 


########################################################################################
# this is another part of the figure1a 

# calculate the median value of age & BMI & depth for each continent 

combine_metadata_statistic_age <- combine_metatable %>% 
        filter(!is.na(age)) %>% 
        group_by(continent) %>% 
        summarise(median_age = median(age), sd_age = sd(age))

combine_metadata_statistic_BMI <- combine_metatable %>%      
        filter(!is.na(BMI)) %>%
        group_by(continent) %>%
        summarise(median_BMI = median(BMI), sd_BMI = sd(BMI))


# depth file 

combine_metadata_statistic_depth <- combine_metatable %>%      
        left_join(., sequencing_depth_file %>% select(internal_sample_id, nonhost)) %>%
        group_by(continent) %>%
        summarise(median_depth = median(nonhost), sd_depth = sd(nonhost))

# gender 
combine_metadata_statistic_gender <- combine_metatable %>% 
        filter(!is.na(gender)) %>% 
	group_by(continent, gender) %>% summarise(num = n()) %>% 
	ungroup()
combine_metadata_statistic_gender_sel <- combine_metatable %>%
        filter(!is.na(gender)) %>%
        group_by(continent) %>% summarise(total_num = n()) %>%
        ungroup() 

combine_metadata_statistic_gender_combine <- combine_metadata_statistic_gender %>%
	left_join(., combine_metadata_statistic_gender_sel)


write_tsv(combine_metadata_statistic_age, str_c(output_dir_plot, "/combine_metadata_statistic_median_age.tsv"))
write_tsv(combine_metadata_statistic_BMI, str_c(output_dir_plot, "/combine_metadata_statistic_median_BMI.tsv"))
write_tsv(combine_metadata_statistic_depth, str_c(output_dir_plot, "/combine_metadata_statistic_median_depth.tsv"))
write_tsv(combine_metadata_statistic_gender_combine, str_c(output_dir_plot, "/combine_metadata_statistic_median_gender.tsv"))


########################################################################################


#                       ribbon plot                                                    # 


########################################################################################


draw_ribbonplot <- function(dataset = NA, color = NA, continent = NA, original_dataset = NA, break_vector = NA){
  
  dataset$Var2 <- as.numeric(dataset$Var2)
  max_table <- original_dataset %>% group_by(feature_group) %>% top_n(value, n = 1) %>% select(value, feature_group) %>% unique()
  print(color)
  p <- ggplot(data = dataset, aes(x = Var2, y = median_feature, ymin = lower_quantile, ymax = upper_quantile, fill = color, group = feature_group)) + 
    geom_line() + 
    geom_ribbon(alpha = 0.5) +
    labs(x = "Number of Samples", y = "Number of Features") +
    ylim(1,3000) +
    theme_bw() +
    scale_x_continuous(breaks= break_vector) +
    theme(axis.title.x = element_text(size = 26, color = "black"), 
          axis.text.x = element_text(size = 24, color = "black"),  
          axis.text.y = element_text(size = 24, color = "black"), 
          axis.title.y = element_text(size = 26, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22),
          legend.key.size = unit(1, "cm"),
          axis.ticks = element_blank()
    ) +
    scale_fill_manual(name = "", values = color) +
    #scale_color_manual(name = "", values = color, limits = continent) +
    theme(legend.position = "none") + 
    annotate("text", x = c(length(unique(original_dataset$Var2)) - 50), y = max_table$value[max_table$feature_group == "EC"], label = "EC", size = 12) + 
    annotate("text", x = c(length(unique(original_dataset$Var2)) - 50), y = max_table$value[max_table$feature_group == "Dietary EC"], label = "Dietary EC", size = 12) + 
    annotate("text", x = c(length(unique(original_dataset$Var2)) - 50), y = max_table$value[max_table$feature_group == "Species"], label = "Species", size = 12)
    return(p)
}

africa_median_quantile <- africa %>%
  group_by(Var1, Var2, feature_group) %>%
  summarise(median_feature = median(value), value_quantile = quantile(value, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

africa_median_quantile_md <- africa_median_quantile %>%
  filter(q != 0.5) %>%
  pivot_wider(names_from = q, values_from = value_quantile) %>%
  rename(lower_quantile = `0.25`, upper_quantile = `0.75`)


america_median_quantile <- america %>%
  group_by(Var1, Var2, feature_group) %>%
  summarise(median_feature = median(value), value_quantile = quantile(value, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

america_median_quantile_md <- america_median_quantile %>%
  filter(q != 0.5) %>%
  pivot_wider(names_from = q, values_from = value_quantile) %>%
  rename(lower_quantile = `0.25`, upper_quantile = `0.75`)


asia_median_quantile <- asia %>%
  group_by(Var1, Var2, feature_group) %>%
  summarise(median_feature = median(value), value_quantile = quantile(value, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

asia_median_quantile_md <- asia_median_quantile %>%
  filter(q != 0.5) %>%
  pivot_wider(names_from = q, values_from = value_quantile) %>%
  rename(lower_quantile = `0.25`, upper_quantile = `0.75`)


europe_median_quantile <- europe %>%
  group_by(Var1, Var2, feature_group) %>%
  summarise(median_feature = median(value), value_quantile = quantile(value, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

europe_median_quantile_md <- europe_median_quantile %>%
  filter(q != 0.5) %>%
  pivot_wider(names_from = q, values_from = value_quantile) %>%
  rename(lower_quantile = `0.25`, upper_quantile = `0.75`)


oceania_median_quantile <- oceania %>%
  group_by(Var1, Var2, feature_group) %>%
  summarise(median_feature = median(value), value_quantile = quantile(value, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))

oceania_median_quantile_md <- oceania_median_quantile %>%
  filter(q != 0.5) %>%
  pivot_wider(names_from = q, values_from = value_quantile) %>%
  rename(lower_quantile = `0.25`, upper_quantile = `0.75`)

africa_p <- draw_ribbonplot(dataset = africa_median_quantile_md, color = "#E64B35FF", continent = "Africa", original_dataset = africa, break_vector = c(0, 150, 300))
america_p <- draw_ribbonplot(dataset = america_median_quantile_md, color = "#4DBBD5FF", continent = "America", original_dataset = america, break_vector = c(0, 300, 600))
asia_p <- draw_ribbonplot(dataset = asia_median_quantile_md, color = "#00A087FF", continent = "Asia", original_dataset = asia, break_vector = c(0, 200, 400))
europe_p <- draw_ribbonplot(dataset = europe_median_quantile_md, color = "#3C5488FF", continent = "Europe", original_dataset = europe, break_vector = c(0, 700, 1400))
oceania_p <- draw_ribbonplot(dataset = oceania_median_quantile_md, color = "#F39B7FFF", continent = "Oceania", original_dataset = oceania, break_vector = c(0, 50, 100))


output_dir <- output_dir_plot
ggsave(str_c(output_dir, "/africa_accumulation_curves_ribbon.pdf"), africa_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/africa_accumulation_curves_ribbon.png"), africa_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/america_accumulation_curves_ribbon.pdf"), america_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/america_accumulation_curves_ribbon.png"), america_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/asia_accumulation_curves_ribbon.pdf"), asia_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/asia_accumulation_curves_ribbon.png"), asia_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/europe_accumulation_curves_ribbon.pdf"), europe_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/europe_accumulation_curves_ribbon.png"), europe_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/oceania_accumulation_curves_ribbon.pdf"), oceania_p, width = 9.49, height = 5.76)
ggsave(str_c(output_dir, "/oceania_accumulation_curves_ribbon.png"), oceania_p, width = 9.49, height = 5.76)


save.image("sample_accumulation_curves_plots.RData")


