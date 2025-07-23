# the script is used to calculate per person compounds number after 1% prevalence filtering 
# food numbers will be in another script (by using the full EC database)
# author:lu zhang 
# date: 2024.05.31

###################################################

# packages 

###################################################


# the idea of the whole script 

# 1. first I need to get the prevalence filtered ECs list from the previous result 



# 2. second I calculate the foods number per samples after prevalence filtering ECs 



# 3. third I need the distribution plot for the foods number 


# 4. Kolmogorov–Smirnov test for distribution profile
 

# 5. normalized by sequencing depth and redo the ks test 


# 6. ks test for the normalized profile for foods

# this script is run by terminal 

############################################################

# loading package                                          #

############################################################

library(ggplot2)
library(rstatix)
library(ggridges)
library(tidyverse)
library(dplyr)


############################################################

# input file                                               #

############################################################

# the modified ec profile per sample 
#modified_ec_profile_after_np <- readRDS(args[1])
modified_ec_profile_after_np <- readRDS("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result1_rarefaction_curve/modified_ec_profile/modified_ec_profile_after_np.rds") # this file, column ec_simplied Compounds are the needed ones.

# the input list of the prevalence filtering from other script : 
#keep_list_1_percent_prevalence <- readRDS(args[2])
keep_list_1_percent_prevalence <- readRDS("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result8_related_smallmolecules/1_percent/food/keep_ec_list_for_each_continents_food.rds")

# metatable
#combine_metadata <- read_tsv(args[3], guess_max = 3000)
combine_metadata <- read_tsv("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result11_the_smallest_species_the_mostEC/compounds/combine_metadata_sel.tsv", guess_max = 3000)


# depth file
#depth_file <- readRDS(args[4])
depth_file <- readRDS("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result5_alpha_diversity_species/non_rarefy/sequencing_depth_revised.rds")

# output_dir 
#output_dir <- args[5]
output_dir <- "/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result8_related_smallmolecules/1_percent/food/"

# full db ec 
#full_ec_db <- read_tsv(args[6])
full_ec_db <- read_tsv("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result12_beneficial_food/nutrichem_md_manual_with_non_microbiotaEC_fullEC.tsv")

############################################################

#  separate the groups into diff.continents                #

############################################################

# different continents 
combine_metadata_Af <- combine_metadata %>%
    dplyr::filter(., continent == "Africa")
combine_metadata_Am <- combine_metadata %>%
    dplyr::filter(., continent == "America")
combine_metadata_As <- combine_metadata %>%
    dplyr::filter(., continent == "Asia")
combine_metadata_Eu <- combine_metadata %>%
    dplyr::filter(., continent == "Europe")
combine_metadata_Oc <- combine_metadata %>%
    dplyr::filter(., continent == "Oceania")

############################################################

#  separate the profile for each continent                 #

############################################################

keep_ec_list_for_continent <- keep_list_1_percent_prevalence

# africa
ec_profile_Af_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Af$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["af"]]$ec_simplied)
ec_profile_Af_con <- ec_profile_Af_con[rowSums(ec_profile_Af_con[-(1:2)])>0,]

# america
ec_profile_Am_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Am$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["am"]]$ec_simplied)
ec_profile_Am_con <- ec_profile_Am_con[rowSums(ec_profile_Am_con[-(1:2)])>0,]

# asia
ec_profile_As_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_As$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["as"]]$ec_simplied)
ec_profile_As_con <- ec_profile_As_con[rowSums(ec_profile_As_con[-(1:2)])>0,]

# europe
ec_profile_Eu_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Eu$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["eu"]]$ec_simplied)
ec_profile_Eu_con <- ec_profile_Eu_con[rowSums(ec_profile_Eu_con[-(1:2)])>0,]

# oceania
ec_profile_Oc_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Oc$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["oc"]]$ec_simplied)
ec_profile_Oc_con <- ec_profile_Oc_con[rowSums(ec_profile_Oc_con[-(1:2)])>0,]


#############################################################

# calculate the compounds number 

#############################################################


calculate_compounds_number <- function(input = NA, sample_names = NA, full_db = full_ec_db){
    
    input_food <- input %>% 
	    left_join(., full_db %>% select(CompoundID, PlantID), by = c("Compounds" = "CompoundID")) %>% 
	    separate_rows(., PlantID, sep = ";")
    
    print(input_food %>% filter(is.na(PlantID)) %>% nrow())

    input_long <- input_food %>%
	dplyr::select(-Compounds) %>% 
        pivot_longer(!c(ec_simplied, PlantID), names_to = "samples", values_to = "cpm") %>%
        unique() # it doesn't matter unique or not, cause i have a unique step in extraction of compound & food

    summary_tibble <- tibble(samples = sample_names, foods_num = -1)
    summary_compounds_food_list <- list()

    for (k in sample_names){
        input_long_sample <- filter(input_long, samples == k) %>% filter(., cpm > 0)
        # compounds number
        foods_sample <- input_long_sample$PlantID %>% unique()
        foods_sample_num <- length(foods_sample)
        summary_compounds_food_list[[k]][["foods"]] <- foods_sample

        summary_tibble[summary_tibble$samples == k,]$foods_num <- foods_sample_num

    }

    summary_list = list(summary_compounds_food_list, summary_tibble) 
    return(summary_list)

}

# calculate each sample corresponding compounds & the related food number depends on the compounds
modified_ec_profile <- list("africa" = ec_profile_Af_con, "america" = ec_profile_Am_con, "asia" = ec_profile_As_con, "europe" = ec_profile_Eu_con, "oceania" = ec_profile_Oc_con)
summary_compounds_table <- list()
for (i in names(modified_ec_profile)){
	print(i)
	sample_names <- colnames(modified_ec_profile[[i]])[-1] #%>% .[!. %in% c("Compounds")]
	summary_compounds_table[[i]] <- calculate_compounds_number(input = modified_ec_profile[[i]], sample_names = sample_names)

}

#############################################################

# combine each continent compounds table 

#############################################################

# summary list 

summary_compounds_overall_table <- do.call(rbind, lapply(summary_compounds_table, `[[`, 2))
summary_compounds_overall_table$group <- "group" 
summary_compounds_overall_table <- summary_compounds_overall_table %>% 
	left_join(., combine_metadata, by = c("samples" = "internal_sample_id"))

summary_compounds_overall_table_sel <- summary_compounds_overall_table %>% 
	filter(samples != "Compounds")

saveRDS(summary_compounds_overall_table_sel, paste(output_dir, "summary_table_continents.rds", sep = "/"))

############################################################

# plot overall distribution of the compounds/foods number
# based on continents

############################################################

print("ploting the distribution compounds number and foods number for each continent... ...")


plot_each_continents_distribution_update_abbreviation <- function(summary_table_prevalence = NA, filetag = NA){
  # for compounds number: color table 
  color_food_tibble <- tibble(continent = c("AF", "AM", "AS", "EU", "OC"),
                                        color = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF"))
  summary_table_prevalence <- summary_table_prevalence %>% 
          mutate(continent = case_when(continent == "Africa" ~ "AF", 
                                       continent == "America" ~ "AM", 
                                       continent == "Asia" ~ "AS", 
                                       continent == "Europe" ~ "EU", 
                                       continent == "Oceania" ~ "OC"))

  summary_table_prevalence_original <- summary_table_prevalence
  statistic <- summary_table_prevalence %>% group_by(continent) %>% summarise(mean = mean(foods_num), sd = sd(foods_num)) # calculate the mean 
  order_statistic <- statistic %>% arrange(desc(mean), continent) %>% left_join(., color_food_tibble) # from the biggest to lowest to order the mean and match color 
  summary_table_prevalence$continent <- factor(summary_table_prevalence$continent, levels = order_statistic$continent) # make continent as factors 
  order_statistic$y <- c(1.5, 2.5, 3.5, 4.5, 5.5) # give the y to the continents 
  print(order_statistic)
  plot_ridges_compounds <- ggplot(summary_table_prevalence, aes(x = foods_num, y = continent, fill = continent)) +
    geom_density_ridges() +
    theme_ridges() +
    scale_x_continuous(name = "N of foods biotransf. by microbiome") + #, breaks = c(0.7,0.72,0.74,0.76,0.78,0.80,0.82), limits = c(0.7,0.83)) +
    theme(legend.position = c(.06, .85)) +
    scale_colour_manual(values=order_statistic$color) +
    scale_fill_manual(values=order_statistic$color) +
    annotate("label", x = 750, y = order_statistic$y[order_statistic$continent == "AF"],
             label = paste("AF:", round(statistic$mean[statistic$continent == "AF"], 0), " \u00B1 ", round(statistic$sd[statistic$continent == "AF"], 0), sep = ""),
             color = "#E64B35FF", size = 4) +
    annotate("label", x = 750, y = order_statistic$y[order_statistic$continent == "AM"],
             label = paste("AM:", round(statistic$mean[statistic$continent == "AM"], 0), " \u00B1 ", round(statistic$sd[statistic$continent == "AM"], 0), sep = ""),
             color = "#4DBBD5FF", size = 4) +
    annotate("label", x = 750, y = order_statistic$y[order_statistic$continent == "AS"],
             label = paste("AS:", round(statistic$mean[statistic$continent == "AS"], 0), " \u00B1 ", round(statistic$sd[statistic$continent == "AS"], 0), sep = ""),
             color = "#00A087FF", size = 4) +
    annotate("label", x = 750, y = order_statistic$y[order_statistic$continent == "EU"],
             label = paste("EU:", round(statistic$mean[statistic$continent == "EU"], 0), " \u00B1 ", round(statistic$sd[statistic$continent == "EU"], 0), sep = ""),
             color = "#3C5488FF", size = 4) +
    annotate("label", x = 750, y = order_statistic$y[order_statistic$continent == "OC"],
             label = paste("OC:", round(statistic$mean[statistic$continent == "OC"], 0), " \u00B1 ",round(statistic$sd[statistic$continent == "OC"], 0), sep = ""),
             color = "#F39B7FFF", size = 4) + ylab(NULL) 

  saveRDS(plot_ridges_compounds, paste(output_dir, paste0(filetag, "_plot_continents_foods_ridges.rds"), sep = "/"))
  pdf(paste(output_dir, paste0(filetag, "_continents_foods_number_distribution_plot.pdf"), sep = "/"))
  print(plot_ridges_compounds)
  dev.off()
  statistics_compounds <- summary_table_prevalence %>%
    wilcox_test(data =., foods_num ~ continent) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")

  result <- list(plot_ridges_compounds, statistics_compounds)
  return(result)

}

overall_distribution_continents <- list()
overall_distribution_continents <- plot_each_continents_distribution_update_abbreviation(summary_table_prevalence = summary_compounds_overall_table_sel, 
												       filetag = "1percent")
############################################################

# the compounds/foods range for the distribution profile  

############################################################

# the range of the compounds and foods of all continents 

summary_table_min_max_0.1 <- summary_compounds_overall_table_sel %>% 
        group_by(continent) %>% 
        summarise(min_food = min(foods_num), max_food = max(foods_num))


min_food <- min(summary_table_min_max_0.1$min_food)
max_food <- max(summary_table_min_max_0.1$max_food)

print(str_c("the min compound number is ", min_food, ",the max compound number is ", max_food))
# [1] "the min compound number is 672,the max compound number is 1035"


# calculate the mode and median for the profile 


summary_0.1_table <- summary_compounds_overall_table_sel %>% dplyr::select(samples, foods_num, continent)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# from this link: https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode

mode_median_summary <- summary_0.1_table %>%
        group_by(continent) %>%
        summarise(foods_num_median = median(foods_num), foods_num_mean = mean(foods_num), foods_num_sd = sd(foods_num), foods_num_mode = Mode(foods_num)) %>%
        ungroup()

write_tsv(mode_median_summary, str_c(output_dir, "/1_percent_mode_median.foods.tsv"))

# asia and oceania median a bit difference 

############################################################

# Kolmogorov–Smirnov test for distribution profile

############################################################

# ks distribution test for the related compounds number 

# reference of doing ks distribution test https://www.statology.org/kolmogorov-smirnov-test-r/

combination_continent <- combn(unique(summary_0.1_table$continent), 2) %>% as.data.frame() %>% t() %>% as.data.frame()  
colnames(combination_continent)[1:2] <- c("compar1", "compar2") 
combination_continent$ks_p <- -1

ks_0.1_table <- list()

for (i in 1:nrow(combination_continent)){

        label <- paste(combination_continent$compar1[i], combination_continent$compar2[i], sep = "_")
        summary_0.1_table_sel_1 <- summary_0.1_table %>%
                filter(continent %in% c(combination_continent$compar1[i]))
        summary_0.1_table_sel_2 <- summary_0.1_table %>%
                filter(continent %in% c(combination_continent$compar2[i]))

        ks_0.1_table[[label]] <- ks.test(summary_0.1_table_sel_1$foods_num, summary_0.1_table_sel_2$foods_num)
        combination_continent$ks_p[i] <- ks_0.1_table[[label]]$p.value
}


write_tsv(combination_continent, str_c(output_dir, "/1_percent_ks_distribution_test.tsv"))

# about ks pvalue
#From the output we can see that the test statistic is 0.99 and the corresponding p-value is 1.299e-14. Since the p-value is less than .05, we reject the null hypothesis. We have sufficient evidence to say that the two sample datasets do not come from the same distribution.


p_no_sig_ks <- combination_continent %>% filter(ks_p > 0.05)
if (nrow(p_no_sig_ks) == 0){

        print("all comparison in ks test are significant different")

}else{
        print("the following are not sig. in ks test")
        print(p_no_sig_ks)

}

print("ks distribution calculation finished")

#[1] "the following are not sig. in ks test"
#   compar1 compar2      ks_p
#V9    Asia Oceania 0.2604606


############################################################

# normalized by sequencing depth & compare again 

############################################################

# normalized the number by the sequencing depth 
new_table <- summary_compounds_overall_table_sel
new_table_combined <- new_table %>% 
        select(samples, foods_num, continent) %>% 
        left_join(., depth_file %>% select(-continent), by = c("samples" = "internal_sample_id"))

new_table_combined$foods_num <- (new_table_combined$foods_num/new_table_combined$nonhost)*1000000


summary_0.1_table_normalized <- new_table_combined

ks_0.1_table_normalized <- list()

combination_continent$ks_p_normalized <- -1

for (i in 1:nrow(combination_continent)){

        label <- paste(combination_continent$compar1[i], combination_continent$compar2[i], sep = "_")
        summary_0.1_table_sel_1 <- summary_0.1_table_normalized %>%
                filter(continent %in% c(combination_continent$compar1[i]))
        summary_0.1_table_sel_2 <- summary_0.1_table_normalized %>%
                filter(continent %in% c(combination_continent$compar2[i]))

        ks_0.1_table[[label]] <- ks.test(summary_0.1_table_sel_1$foods_num, summary_0.1_table_sel_2$foods_num)
        combination_continent$ks_p_normalized[i] <- ks_0.1_table[[label]]$p.value
}



write_tsv(combination_continent, str_c(output_dir, "1_percent_ks_distribution_test_normalized.tsv"))
saveRDS(new_table_combined, str_c(output_dir, "1_percent_normalized_foods_table.rds"))


p_no_sig_ks_normalized <- combination_continent %>% filter(ks_p_normalized > 0.05)
if (nrow(p_no_sig_ks_normalized) == 0){

        print("all comparison in ks test are significant different")

}else{
        print("the following are not sig. in ks test")
        print(p_no_sig_ks_normalized)

}

print("ks distribution calculation for normalized profile finished")
# [1] "all comparison in ks test are significant different"


############################################################

# density plot 

############################################################


# draw density plot for each continent
print("drawing the density plot for each continent...")

plot_each_continents_distribution_update_normalized <- function(summary_table_prevalence = NA, filetag = NA, x_position = NA, round_type = NA){
  # for compounds number
  # colortable
  color_food_tibble <- tibble(continent = c("AF", "AM", "AS", "EU", "OC"),
                                        color = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF"))

  summary_table_prevalence <- summary_table_prevalence %>%
          mutate(continent = case_when(continent == "Africa" ~ "AF",
                                       continent == "America" ~ "AM",
                                       continent == "Asia" ~ "AS",
                                       continent == "Europe" ~ "EU",
                                       continent == "Oceania" ~ "OC"))

  summary_table_prevalence_original <- summary_table_prevalence
  statistic <- summary_table_prevalence %>% group_by(continent) %>% summarise(mean = mean(foods_num), sd = sd(foods_num))
  order_statistic <- statistic %>% arrange(desc(mean), continent) %>% left_join(., color_food_tibble)
  summary_table_prevalence$continent <- factor(summary_table_prevalence$continent, levels = order_statistic$continent)
  order_statistic$yaxis <- c(1.5, 2.5, 3.5, 4.5, 5.5)
  pdf(paste(output_dir, paste0(filetag, "_continents_foods_number_distribution_plot_update_normalized.pdf"), sep = "/"))
  plot_ridges_compounds <- ggplot(summary_table_prevalence, aes(x = foods_num, y = continent, fill = continent)) +
    geom_density_ridges() +
    theme_ridges() +
    scale_x_continuous(name = "N of foods biotransf. by microbiome/Million reads") + #, breaks = c(0.7,0.72,0.74,0.76,0.78,0.80,0.82), limits = c(0.7,0.83)) +
    theme(legend.position = c(.85, .85)) +
    scale_colour_manual(values=order_statistic$color) +
    scale_fill_manual(values=order_statistic$color) +
    annotate("label", x = x_position, y = order_statistic$yaxis[order_statistic$continent == "AF"],
             label = paste("AF:", round(statistic$mean[statistic$continent == "AF"], round_type), " \u00B1 ", 
                           round(statistic$sd[statistic$continent == "AF"], round_type), sep = ""),
             color = "#E64B35FF", size = 4) +
    annotate("label", x = x_position, y = order_statistic$yaxis[order_statistic$continent == "AM"],
             label = paste("AM:", round(statistic$mean[statistic$continent == "AM"], round_type), " \u00B1 ", 
                           round(statistic$sd[statistic$continent == "AM"], round_type), sep = ""),
             color = "#4DBBD5FF", size = 4) +
    annotate("label", x = x_position, y = order_statistic$yaxis[order_statistic$continent == "AS"],
             label = paste("AS:", round(statistic$mean[statistic$continent == "AS"], round_type), " \u00B1 ", 
                           round(statistic$sd[statistic$continent == "AS"], round_type), sep = ""),
             color = "#00A087FF", size = 4) +
    annotate("label", x = x_position, y = order_statistic$yaxis[order_statistic$continent == "EU"],
             label = paste("EU:", round(statistic$mean[statistic$continent == "EU"], round_type), " \u00B1 ", 
                           round(statistic$sd[statistic$continent == "EU"], round_type), sep = ""),
             color = "#3C5488FF", size = 4) +
    annotate("label", x = x_position, y = order_statistic$yaxis[order_statistic$continent == "OC"],
             label = paste("OC:", round(statistic$mean[statistic$continent == "OC"], round_type), " \u00B1 ",
                           round(statistic$sd[statistic$continent == "OC"], round_type), sep = ""),
             color = "#F39B7FFF", size = 4) + ylab(NULL)

  saveRDS(plot_ridges_compounds, paste(output_dir, paste0(filetag, "_plot_continents_foods_ridges_normalized.rds"), sep = "/"))
  print(plot_ridges_compounds)
  dev.off()
  statistics_compounds <- summary_table_prevalence %>%
    wilcox_test(data =., foods_num ~ continent) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  write_tsv(statistics_compounds, paste(output_dir, paste0(filetag, "statistics_foods_continents_update_normalized.tsv"), sep = "/"))

  result <- list(plot_ridges_compounds, statistics_compounds)
  return(result)

}


overall_distribution_continents_update_normalized <- plot_each_continents_distribution_update_normalized(summary_table_prevalence = new_table_combined, 
														       filetag = "1_percent", x_position = 100, round_type = 0)


print("saving RData...")
save.image(str_c(output_dir, "/step19_figuremodification_manual_check_1percent_food.RData"))


print("foods number finished.")





