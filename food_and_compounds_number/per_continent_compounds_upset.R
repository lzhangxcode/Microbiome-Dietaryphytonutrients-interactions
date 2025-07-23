# the aim of the script is to do 1% prevalence filtering for the samples & draw upset plot for compounds 
# author: lu.zhang 
# date: 2024.Mai.21 

######################################################################

# Data input 

######################################################################

library(tidyverse)
# modified ec profile 
modified_ec_profile_after_np <- readRDS(args[1])


# metatable
combine_metadata <- read_tsv(args[2], guess_max = 3000)

# output 
output_dir <- args[3]

######################################################################

# acquire the 1% prevalence filtering EC 

######################################################################

prevalence_filtering <- function(input_modified_ec = NA, prevalence = NA){
  colnames(input_modified_ec)[1] <- "Genefamilies" # it is actually ec, i just too lazy to change names
  sample_num <- input_modified_ec %>% dplyr::select(-Compounds, -Genefamilies, -Tag) %>% ncol()
  input_modified_ec_long <- input_modified_ec %>%
    dplyr::select(-Compounds,-Tag) %>%
    unique() %>%
    pivot_longer(!Genefamilies) %>%
    add_column(num = ifelse(.$value>0,1,0)) %>%
    group_by(Genefamilies) %>%
    summarize(sum = sum(num)) %>%
    mutate(freq = sum / sample_num) %>%
    ungroup()


  keep_ec <- input_modified_ec_long %>%
    filter(freq >= prevalence)

  input_modified_ec_sel <- input_modified_ec %>% filter(Genefamilies %in% keep_ec$Genefamilies)

  return(input_modified_ec_sel)

}


modified_ec_profile_after_np_1_percent <- prevalence_filtering(input_modified_ec = modified_ec_profile_after_np, prevalence = 0.01)
saveRDS(modified_ec_profile_after_np_1_percent, str_c(output_dir, "keep_list_1_percent_allsamples.rds"))


######################################################################

# filter each continent the ECs 

######################################################################

# for each continent separate the profile 

table(combine_metadata$continent)
#Africa America    Asia  Europe Oceania 
#    326     784     476    1379     103 


ec_profile <- modified_ec_profile_after_np_1_percent
colnames(ec_profile)[1] <- "ec_simplied"

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

######################################################################

# 1% prevalence per continent 

######################################################################

prevalence_filtering_percontinent <- function(input_modified_ec = modified_ec_profile_after_np, metatable = combine_metadata, prevalence = 0.01){

	# continent number 
	continent_number <- data.frame(table(metatable$continent))
	# filter each continent 
	input_modified_ec_long <- input_modified_ec %>%
    		dplyr::select(-Compounds,-Tag) %>%
    		unique() %>%
    		pivot_longer(!ec_simplied) %>%
    		add_column(num = ifelse(.$value>0,1,0)) %>%
		left_join(., metatable, by = c("name" = "internal_sample_id")) %>% 
    		group_by(continent, ec_simplied) %>%
    		summarize(sum = sum(num)) %>%
		left_join(., continent_number, by = c("continent" = "Var1")) %>%
    		mutate(freq = sum / Freq) %>%
    		ungroup()
	
	# filter ec > prevalence 
        keep_ec <- input_modified_ec_long %>%
    		filter(freq >= prevalence)

	
	# produce 
	keep_ec_list <- list()
	keep_ec_list[["af"]] <- keep_ec %>% filter(continent == "Africa")
	keep_ec_list[["am"]] <- keep_ec %>% filter(continent == "America")
	keep_ec_list[["as"]] <- keep_ec %>% filter(continent == "Asia")
	keep_ec_list[["oc"]] <- keep_ec %>% filter(continent == "Oceania")
	keep_ec_list[["eu"]] <- keep_ec %>% filter(continent == "Europe")

	# return
	return(keep_ec_list)

}



# 1% prevalence per continent 

keep_ec_list_for_continent <- prevalence_filtering_percontinent(input_modified_ec = modified_ec_profile_after_np, prevalence = 0.01)


######################################################################

# filter each continent ec  

######################################################################
# 1% prevalence per continent 

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


######################################################################

# acquire the compound list - percontinent  

######################################################################
# 1 percent 
listInput_con <- list(AF = ec_profile_Af_con$Compounds %>% unique(),
                  AM = ec_profile_Am_con$Compounds %>% unique(),
                  AS = ec_profile_As_con$Compounds %>% unique(),
                  EU = ec_profile_Eu_con$Compounds %>% unique(),
                  OC = ec_profile_Oc_con$Compounds %>% unique())

# calculate the compounds number (for the order in the upset plot)

africa_comp_num_con <- length(ec_profile_Af_con$Compounds %>% unique())
america_comp_num_con <- length(ec_profile_Am_con$Compounds %>% unique())
asia_comp_num_con <- length(ec_profile_As_con$Compounds %>% unique())
europe_comp_num_con <- length(ec_profile_Eu_con$Compounds %>% unique())
oceania_comp_num_con <- length(ec_profile_Oc_con$Compounds %>% unique())


# draw plot 
library(UpSetR)

color_compound_tibble_con <- tibble(continent = c("AF", "AM", "AS", "EU", "OC"),
                                        color = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF"),
                                        number = c(africa_comp_num_con, america_comp_num_con, asia_comp_num_con, europe_comp_num_con, oceania_comp_num_con)) %>%
                arrange(desc(number))



pdf(paste(output_dir, "upset_compound_1percent_con.pdf", sep = "/"))#, width = 6, height = 5)
text_scale_options3 <- c(2.5, 2, 2, 2, 2.5, 3)
compounds_upset <- UpSetR::upset(fromList(listInput_con), sets = c("AF", "AM", "AS", "EU", "OC"),
                                   sets.x.label = "LMWMs",
                                   mainbar.y.label = "N of phytonutrients biotransf. \nby microbiome",
                                   mb.ratio = c(0.65, 0.35), point.size = 4,
                                   order.by = "freq", #empty.intersections = "on",
                                   text.scale = text_scale_options3, sets.bar.color = color_compound_tibble_con$color)
print(compounds_upset)
dev.off()



######################################################################

# save RData 

######################################################################

save.image(str_c(output_dir, "/compounds_number_1_percent.RData"))


# save the list into rds
saveRDS(keep_ec_list_for_continent, str_c(output_dir, "/keep_ec_list_for_each_continents.rds"))
saveRDS(keep_ec_list_for_continent_5percent, str_c(output_dir, "/keep_ec_list_for_each_continents_5percent.rds"))


