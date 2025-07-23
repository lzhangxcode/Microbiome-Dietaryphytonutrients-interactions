# the aim of the script is to do beta diversity for the phytonutrients-associated ec profile 
# date: 22,Dec,2023 
# author: l.zhang 


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 6) {
                  stop("Rscript ../step16_adonis_modified_profile_manual_check.R  modified_ec_profile depth_file sequencing_platform metatable nutrichem output_dir", call.=FALSE)
}



############################################################

# loading package                                          #

############################################################

library(vegan)
library(tidyverse)


############################################################

# input file                                               #

############################################################

# modified ec profile after natural product check 
modified_ec_profile_after_np <- readRDS(args[1])

# sequencing depth file 
depth_file <- readRDS(args[2])

# sequencing platform file 
sequencing_platform <- read_tsv(args[3])

# combine metatable 
combine_metadata <- read_tsv(args[4], guess_max = 3000)

# nutrichem class file 
nutrichem_md <- read_tsv(args[5])

# output dir 
output_dir <- args[6]


############################################################

# update metatable sequencing platform                     #

############################################################

# produce the metatable with age, gender & BMI information 
produce_new_metadata <- function(metadata = NA, sequencing_depth_file = NA, sequencing_platform_manually_check = NA){

	metadata_input <- metadata %>%
        	select(internal_sample_id, age, BMI, gender, continent, study_name, sequencing_platform, country, subject_id)
	sequencing_depth_sel <- sequencing_depth_file %>% 
        	select(internal_sample_id, input, kept, fwd_only, rev_only, dropped, host, nonhost, library)
	metadat <- metadata_input

    	metadat_sel <-  metadat %>%
        	filter(., !is.na(age) & !is.na(BMI) & !is.na(gender)) %>% 
        	left_join(., sequencing_depth_sel, by = "internal_sample_id")
	
	# update the sequencing platform
        colnames(sequencing_platform_manually_check)[3] <- "manually_checked_platform"

    	for (study in sequencing_platform_manually_check$study_name){
        	print(study)
        	if (study == "LeChatelierE_2013"){
            		print("don't change")
        	}else{
            		platform <- sequencing_platform_manually_check[sequencing_platform_manually_check$study_name == study,]$manually_checked_platform %>% unique()
            		print(platform)
            		platform_original <- metadat_sel[metadat_sel$study_name == study,]$sequencing_platform %>% unique()
            		print(platform_original)
            		metadat_sel[metadat_sel$study_name == study,]$sequencing_platform <- platform
	        }
   
    	}

	return(metadat_sel)

}

metadata_sel_revised <- produce_new_metadata(metadata = combine_metadata,
                                            sequencing_depth_file = depth_file,
                                            sequencing_platform_manually_check = sequencing_platform)

print(str_c("in total,", nrow(metadata_sel_revised), "samples with age, gender & BMI information"))



############################################################

# input profile for beta diversity                         #

############################################################

prepare_beta_diversity_input <- function(input = NA, metadat = NA){
 
    beta_diversity_input <- input %>% 
	dplyr::select(-Compounds,-Tag) %>% 
        dplyr::select(., c("ec_simplied", metadat$internal_sample_id)) %>%  
        unique() %>%
	column_to_rownames(., var = "ec_simplied")

    beta_diversity_input <- beta_diversity_input[rowSums(abs(beta_diversity_input)) != 0,]
    

    beta_diversity_input_t <- beta_diversity_input %>% t(.)
     

    input_list <- list()
    metadat$BMI <- as.numeric(metadat$BMI)
    input_list[[1]] <- metadat
    input_list[[2]] <- beta_diversity_input_t

    return(input_list)
}

modified_ec_profile_beta_input  <- prepare_beta_diversity_input(input = modified_ec_profile_after_np, metadat = metadata_sel_revised)

print(str_c("the modified ec profile has ", ncol(modified_ec_profile_beta_input[[2]]), " ECs"))
print(str_c("the modified ec profile has ", nrow(modified_ec_profile_beta_input[[2]]), " samples"))


############################################################

# adonis statistical results                               #

############################################################

# produce adonis result 
distance_produce <- function(beta_diversity = NA){ # 1st list metadata; 2nd list ec_profile

    metadat_sel <- beta_diversity[[1]]
    metadat_sel$depth <- metadat_sel$nonhost
    profile <- beta_diversity[[2]]


    # distance calculation
    set.seed(1234)
    distance.bray <- profile %>%
        vegdist(.,method = 'bray') %>%
        as.matrix()

    set.seed(1234)
    aitchson_input <- profile
    min_tmp <- min(aitchson_input[aitchson_input>0])
    aitchson_input[aitchson_input == 0] <- min_tmp/2
    distance.aitchison <- as.matrix(robCompositions::aDist(aitchson_input))



    set.seed(1234)
    permanova_result_bray <- vegan::adonis2(distance.bray ~ age + gender + BMI + continent, data = metadat_sel)
    set.seed(1234)
    permanova_result_aitchison <- vegan::adonis2(distance.aitchison ~ age + gender + BMI + continent, data = metadat_sel)

     # separate adonis for each factor - bray curtis 
    permanova_result_bray_continent <- vegan::adonis2(distance.bray ~ continent, data = metadat_sel)
    permanova_result_bray_country <- vegan::adonis2(distance.bray ~ country, data = metadat_sel)
    permanova_result_bray_gender <- vegan::adonis2(distance.bray ~ gender, data = metadat_sel)
    permanova_result_bray_age <- vegan::adonis2(distance.bray ~ age, data = metadat_sel)
    permanova_result_bray_BMI <- vegan::adonis2(distance.bray ~ BMI, data = metadat_sel)
    permanova_result_bray_depth <- vegan::adonis2(distance.bray ~ depth, data = metadat_sel)
    permanova_result_bray_sequencingplatform <- vegan::adonis2(distance.bray ~ sequencing_platform, data = metadat_sel)
    permanova_result_bray_study <- vegan::adonis2(distance.bray ~ study_name, data = metadat_sel)

    # separate adonis for each factor - aitchison
    permanova_result_aitchison_continent <- vegan::adonis2(distance.aitchison ~ continent, data = metadat_sel)
    permanova_result_aitchison_country <- vegan::adonis2(distance.aitchison ~ country, data = metadat_sel)
    permanova_result_aitchison_gender <- vegan::adonis2(distance.aitchison ~ gender, data = metadat_sel)
    permanova_result_aitchison_age <- vegan::adonis2(distance.aitchison ~ age, data = metadat_sel)
    permanova_result_aitchison_BMI <- vegan::adonis2(distance.aitchison ~ BMI, data = metadat_sel)
    permanova_result_aitchison_depth <- vegan::adonis2(distance.aitchison ~ depth, data = metadat_sel)
    permanova_result_aitchison_sequencingplatform <- vegan::adonis2(distance.aitchison ~ sequencing_platform, data = metadat_sel)
    permanova_result_aitchison_study <- vegan::adonis2(distance.aitchison ~ study_name, data = metadat_sel)

    statistic_result <- list()
    statistic_result[["bray_distance"]] <- distance.bray
    statistic_result[["aitchison"]] <- distance.aitchison
    statistic_result[["bray_adonis"]] <- permanova_result_bray
    statistic_result[["aitchison_adonis"]] <- permanova_result_aitchison

    statistic_result[["bray_adonis_continent"]] <- permanova_result_bray_continent
    statistic_result[["bray_adonis_country"]] <- permanova_result_bray_country
    statistic_result[["bray_adonis_gender"]] <- permanova_result_bray_gender
    statistic_result[["bray_adonis_age"]] <- permanova_result_bray_age
    statistic_result[["bray_adonis_BMI"]] <- permanova_result_bray_BMI
    statistic_result[["bray_adonis_depth"]] <- permanova_result_bray_depth
    statistic_result[["bray_adonis_sequencingplatform"]] <- permanova_result_bray_sequencingplatform
    statistic_result[["bray_adonis_study"]] <- permanova_result_bray_study

    statistic_result[["aitchison_adonis_continent"]] <- permanova_result_aitchison_continent
    statistic_result[["aitchison_adonis_country"]] <- permanova_result_aitchison_country
    statistic_result[["aitchison_adonis_gender"]] <- permanova_result_aitchison_gender
    statistic_result[["aitchison_adonis_age"]] <- permanova_result_aitchison_age
    statistic_result[["aitchison_adonis_BMI"]] <- permanova_result_aitchison_BMI
    statistic_result[["aitchison_adonis_depth"]] <- permanova_result_aitchison_depth
    statistic_result[["aitchison_adonis_sequencingplatform"]] <- permanova_result_aitchison_sequencingplatform
    statistic_result[["aitchison_adonis_study"]] <- permanova_result_aitchison_study


    return(statistic_result)

}

    
modified_ec_profile_distance_result <- distance_produce(beta_diversity = modified_ec_profile_beta_input)
saveRDS(modified_ec_profile_beta_input,str_c(output_dir, "modified_ec_profile_beta_input.rds"))



############################################################

# PCoA plot                                                #

############################################################

library(ggpubr)


PcoA_plot_continent_bray <- function(points_df = NA, eig = NA, R2 = NA, p = NA, disp_p = NA){
    color <-  c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")
    continent <- c("Africa", "America", "Asia", "Europe")

    annotate_label <- paste0("PERMANOVA, p = ", p, ", effect size =  ", round(R2, 3))

    pcoa_plot <- ggscatter(points_df, x = "x", y = "y",
               color = "continent", ellipse = TRUE) + #, ellipse = TRUE, ellipse.type = "convex") +
    scale_fill_manual(name = "Enzyme of Diff.Continents", values = color, limits = continent) + 
    scale_color_manual(name = "Enzyme of Diff.Continents", values = color, limits = continent) + 
    annotate("text", x = -0.22, y = 0.22, label = annotate_label, size = 7) +
    annotate("text", x = -0.32, y = 0.19, label = paste0("Beta dispersion p = ", disp_p), size = 7) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
    theme(text = element_text(size=20),
        axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


    return(pcoa_plot)
}

PcoA_plot_continent_ait <- function(points_df = NA, eig = NA, R2 = NA, p = NA, disp_p = NA){
    color <-  c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")
    continent <- c("Africa", "America", "Asia", "Europe")

    annotate_label <- paste0("PERMANOVA, p = ", p, ", effect size =  ", round(R2, 3))

    pcoa_plot <- ggscatter(points_df, x = "x", y = "y",
               color = "continent", ellipse = TRUE) + #, ellipse = TRUE, ellipse.type = "convex") +
    scale_fill_manual(name = "Enzyme of Diff.Continents", values = color, limits = continent) +
    scale_color_manual(name = "Enzyme of Diff.Continents", values = color, limits = continent) +
    annotate("text", x = -80, y = 65, label = annotate_label, size = 7) +
    annotate("text", x = -117, y = 58, label = paste0("Beta dispersion p = ", disp_p), size = 7) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    theme(text = element_text(size=20),
        axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

    return(pcoa_plot)
}

PCoA_plots_all_modify <- function(distance.bray = NA, permanova_result = NA, output_dir = NA, metadat_sel = NA, filetag = NA, R2 = NA, p = NA, disp_p = NA, bray = "yes"){

    pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
    points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
    colnames(points) <- c("x", "y", "z")
    eig <- pcoa$eig
    points <- cbind(points, metadat_sel)


    if (bray == "yes"){
            write_rds(pcoa, paste(output_dir, "pcoa_bray.rds", sep = "/"))
    continent_pcoa <- PcoA_plot_continent_bray(points = points, eig = eig, R2 = R2, p = p, disp_p = disp_p)
    }else {
            write_rds(pcoa, paste(output_dir, "pcoa_ait.rds", sep = "/"))
    continent_pcoa <- PcoA_plot_continent_ait(points = points, eig = eig, R2 = R2, p = p, disp_p = disp_p)
    }
    continent_pcoa <- ggExtra::ggMarginal(continent_pcoa, type="boxplot", groupColour = TRUE, groupFill = TRUE, size = 8)

    return(continent_pcoa)
}


# calculate the beta dispersion 
modified_metadata <- modified_ec_profile_beta_input[[1]]
bray_distance <- as.dist(modified_ec_profile_distance_result[["bray_distance"]])
aitchison <- as.dist(modified_ec_profile_distance_result[["aitchison"]])


modified_bray.bd <- betadisper(bray_distance, modified_metadata$continent)
modified_aitchison.bd <- betadisper(aitchison, modified_metadata$continent)

a <- permutest(modified_bray.bd)
b <- permutest(modified_aitchison.bd)
saveRDS(a, paste(output_dir, "dietaryEC_betadisp_bray.rds"))
saveRDS(b, paste(output_dir, "dietaryEC_betadisp_aitchison.rds"))


beta_disper_bray <- a
beta_disper_aitchison <- b 


distance_result <- modified_ec_profile_distance_result 
metadat_sel <- modified_ec_profile_beta_input[[1]]
pcoa_result_bray <- PCoA_plots_all_modify(distance.bray = distance_result[["bray_distance"]], 
                                   permanova_result = distance_result[["bray_adonis"]], 
                                   output_dir = output_dir, 
                                   metadat_sel = metadat_sel, 
                                   filetag = paste(category, "bray_curtis", sep = ""),
                                   R2 = distance_result[["bray_adonis_continent"]]$R2[1],
                                   p = distance_result[["bray_adonis_continent"]]$`Pr(>F)`[1],
                                   disp_p = beta_disper_bray$tab$`Pr(>F)`[1])




ggsave(paste(output_dir, "modified_pcoa_plot_with_boxplot.pdf", sep = "/"), pcoa_result_bray, width = 10, height = 7)

pcoa_result_aitchison <- list()
pcoa_result_aitchison <- PCoA_plots_all_modify(distance.bray = distance_result[["aitchison"]], 
                                   permanova_result = distance_result[["aitchison_adonis"]], 
                                   output_dir = output_dir, 
                                   metadat_sel = metadat_sel,
                                   filetag = paste(category, "aitchison", sep = ""),
                                   R2 = distance_result[["aitchison_adonis_continent"]]$R2[1],
                                   p = distance_result[["aitchison_adonis_continent"]]$`Pr(>F)`[1],
                                   disp_p = beta_disper_aitchison$tab$`Pr(>F)`[1],
                                   bray = "no")

ggsave(paste(output_dir, "modified_pcoa_plot_with_boxplot_aitchison.pdf", sep = "/"), pcoa_result_aitchison, height = 7, width = 10)

print("finished drawing pcoa plot with boxplot")

############################################################

# PCoA plot by adding env from compound class              #

############################################################


ec_profile_sel <- modified_ec_profile_after_np %>%  
	dplyr::select(., c("ec_simplied", metadata_sel_revised$internal_sample_id, "Compounds", "Tag"))

##   calculate the compound class for each sample
# =================================================
compound_class_table <- ec_profile_sel %>% 
        select(Compounds) %>% 
        unique() %>% 
        left_join(., nutrichem_md %>% select(CompoundID, CompoundClass, class3, class1, class2) %>% unique(), by = c("Compounds" = "CompoundID"))

## here we calculated how many ecs are related to one compound 
for (i in colnames(ec_profile_sel)[!colnames(ec_profile_sel) %in% c("ec_simplied", "Compounds", "Tag")]){

        print(i)

        target_sample <- ec_profile_sel %>% select(i, Compounds, ec_simplied) 
        colnames(target_sample)[1] <- "sample"
        target_sample_unique <- target_sample %>%
                filter(sample > 0) %>% 
                select(ec_simplied, Compounds) %>% 
                unique() %>% 
                group_by(Compounds) %>% 
                summarise(exist = n()) %>% 
                ungroup()
        colnames(target_sample_unique)[2] <- i
        
        compound_class_table <- compound_class_table %>% 
                left_join(., target_sample_unique)

}

## for class1 produce number of compouds 
# =================================================
compound_class_table[is.na(compound_class_table)] <- 0

compound_class_table_class1 <- compound_class_table %>%
        #select(-CompoundClass, -Compounds, -class3, -class2) %>% 
        pivot_longer(A01_02_1FE:HV9) %>% 
        group_by(class1, name) %>% 
        summarise(sum_value = sum(value)) %>% 
        ungroup() %>% 
        pivot_wider(names_from = name, values_from = sum_value) %>% 
        column_to_rownames(var = "class1")

compound_class_table_class1_t <- compound_class_table_class1 %>% 
        t(.) %>% as.data.frame()

bray_curtis_distance <- modified_ec_profile_distance_result[["bray_distance"]]
aitchison_distance <- modified_ec_profile_distance_result[["aitchison"]]


ec_pcoa <- list()
ec_env <- list()

library(vegan)
set.seed(1234)
#ec_pcoa[["modified_bray"]] <- cmdscale(bray_curtis_distance, k = (nrow(bray_curtis_distance) - 1), eig = TRUE)
ec_pcoa[["modified_bray"]] <- cmdscale(bray_curtis_distance, k = 3, eig = TRUE)
set.seed(1234)
ec_env[["modified_bray"]] <- envfit(ec_pcoa[["modified_bray"]]~., data = compound_class_table_class1_t, perm = 999, choices = c(1,2), display = 'sites')



set.seed(1234)
#ec_pcoa[["modified_ait"]] <- cmdscale(aitchison_distance, k = (nrow(aitchison_distance) - 1), eig = TRUE)
ec_pcoa[["modified_ait"]] <- cmdscale(aitchison_distance, k = 3, eig = TRUE)
set.seed(1234)
ec_env[["modified_ait"]] <- envfit(ec_pcoa[["modified_ait"]]~., data = compound_class_table_class1_t,, perm = 999, choices = c(1,2), display = 'sites')


## select top10 for visualization  
# =================================================

# how many class1 

# add label 

annotate_label_bray <- paste0("PERMANOVA, p = ", distance_result[["bray_adonis_continent"]]$`Pr(>F)`[1], 
                         ", effect size =  ", round(distance_result[["bray_adonis_continent"]]$R2[1], 3))
annotate_label_ait <- paste0("PERMANOVA, p = ", distance_result[["aitchison_adonis_continent"]]$`Pr(>F)`[1],
                        ", effect size =  ", round(distance_result[["aitchison_adonis_continent"]]$R2[1], 3))



# draw plot 
draw_pcoa_plot_with_species <- function(pcoa_site_merge = NA, pcoa_env_sign_order_top10 = NA, pc1 = NA, pc2 = NA, size_axis = NA, distance = NA){

  if (distance == "bray"){
        x_label = -0.22
        y_label = 0.22
        annotate_label <- annotate_label_bray
  }else{
        x_label = -80
        y_label = 65
        annotate_label <- annotate_label_ait

  }

  pcoa_site_merge$continent[pcoa_site_merge$continent == "Africa"] <- "AF"
  pcoa_site_merge$continent[pcoa_site_merge$continent == "America"] <- "AM"
  pcoa_site_merge$continent[pcoa_site_merge$continent == "Asia"] <- "AS"
  pcoa_site_merge$continent[pcoa_site_merge$continent == "Europe"] <- "EU"
  color <-  c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF") # , "#F39B7FFF"
  continent <- c("AF", "AM", "AS", "EU") # , "Oceania"
  p_pcoa_env <- ggscatter(pcoa_site_merge, x = "V1", y = "V2", color = "continent", ellipse = TRUE) +
    scale_fill_manual(name = "Enzyme of Diff.Continents", values = color, limits = continent) +
    scale_color_manual(name = "Enzyme of Diff.Continents", values = color, limits = continent) +
    #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) +
    labs(x = paste('PCoA1 (', pc1, '%)'), y = paste('PCoA2 (', pc2, '%)')) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
    geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    geom_segment(data = pcoa_env_sign_order_top10, aes(x = 0,y = 0, xend = PC1/size_axis, yend = PC2/size_axis), arrow = arrow(length=unit(0.1, 'cm')), size=0.4, color='black') +
    geom_text_repel(data = pcoa_env_sign_order_top10, aes(PC1/size_axis * 1.0, PC2/size_axis * 1.0, label = group), color = 'black', size = 4, max.overlaps = 50) +
    annotate("text", x = x_label, y = y_label, label = annotate_label, size = 7) +
    theme(text = element_text(size=20),
        axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

   p_pcoa_env <- ggExtra::ggMarginal(p_pcoa_env, type="boxplot", groupColour = TRUE, groupFill = TRUE, size = 8)

  return(p_pcoa_env)

}

produce_env <- function(input_pcoa_bray = NA, input_pcoa_ef = NA, combine_metadata = NA, output_dir = NA, filetag = NA, size_axis = NA, dis = NA){
  
  # first part: adjust env pvalue to p.adusted fdr value 
  pcoa_ef_adj <- input_pcoa_ef
  pcoa_ef_adj$vectors$pvals <- p.adjust(pcoa_ef_adj$vectors$pvals, method = 'fdr')
  
  
  # second part: save it to csv file 
  pcoa_env <- data.frame(cbind(pcoa_ef_adj$vectors$arrows, pcoa_ef_adj$vectors$r, pcoa_ef_adj$vectors$pvals))
  names(pcoa_env) <- c('PC1', 'PC2', 'r2', 'p.adj')
  write.csv(pcoa_env, paste0(output_dir, filetag, '_pcoa_env.csv'), quote = FALSE)
  
  
  # third part: axis explained variance 
  
  eig <- input_pcoa_bray$eig
  pc1 <- format(100 * eig[1] / sum(eig), digits=4)
  pc2 <- format(100 * eig[2] / sum(eig), digits=4)
  
  
  # forth part: the metadata prepare 
  
  sample_continent <- tibble(name = rownames(input_pcoa_bray$points)) %>% left_join(., combine_metadata, by = c("name" = "internal_sample_id")) %>% select(name, continent)
  
 
  
  # fifth part: env (species pick)
  
  pcoa_env_sign <- pcoa_env[which(pcoa_env$p.adj < 0.05),] # choose significant 
  pcoa_env_sign_order <- pcoa_env_sign %>% arrange(desc(r2)) 

  pcoa_env_sign_order <- pcoa_env_sign_order[!rownames(pcoa_env_sign_order) %in% c("Others", "Unclassified", "Others;Others", "Unclassified;", "Others;Others;Others"),]

  if (nrow(pcoa_env_sign_order) > 20){

        
        pcoa_env_sign_order_top20 <- pcoa_env_sign_order[1:20,] # choose the biggest top20 r2 
        pcoa_env_sign_order_top20$group <- rownames(pcoa_env_sign_order_top20) #%>% str_split(., pattern = "\\|", simplify = T) %>% .[,7]


  }

 if (nrow(pcoa_env_sign_order) > 10){

        pcoa_env_sign_order_top10 <- pcoa_env_sign_order[1:10,]
        pcoa_env_sign_order_top10$group <- rownames(pcoa_env_sign_order[1:10,])

  }else{

        pcoa_env_sign_order_top10 <- pcoa_env_sign_order
        pcoa_env_sign_order_top10$group <- rownames(pcoa_env_sign_order)
  }
  
   
 
  # sixth part: combine pcoa axis with continent info 
  pcoa_site <- input_pcoa_bray$points[,1:2] %>% as.data.frame(.)
  pcoa_site$samples <- rownames(input_pcoa_bray$points)
  pcoa_site_merge <- pcoa_site %>% left_join(., sample_continent, by = c("samples" = "name")) %>% select(V1, V2, samples, continent)
  
  
  # seventh part: draw plot 
  pcoa_plot <- draw_pcoa_plot_with_species(pcoa_site_merge = pcoa_site_merge, pcoa_env_sign_order_top10 = pcoa_env_sign_order_top10, pc1 = pc1, pc2 = pc2, size_axis = size_axis,
distance = dis)
  ggsave(file = paste0(output_dir, filetag, "_pcoa_env_plot.pdf"), pcoa_plot, width = 14, height = 8)
  saveRDS(pcoa_site_merge, paste0(output_dir, filetag, "_pcoa_site.rds"))
  
  return(pcoa_plot)
  
}

plot_pcoa_list <- list()
library(ggrepel)
plot_pcoa_list[["bray"]] <- produce_env(input_pcoa_bray = ec_pcoa[["modified_bray"]],
                                           input_pcoa_ef = ec_env[["modified_bray"]],
                                           combine_metadata = combine_metadata,
                                           output_dir = output_dir,
                                           filetag = "modified_bray_compoundclass_",
                                           size_axis = 5,
                                           dis = "bray")

plot_pcoa_list[["ait"]] <- produce_env(input_pcoa_bray = ec_pcoa[["modified_ait"]],
                                           input_pcoa_ef = ec_env[["modified_ait"]],
                                           combine_metadata = combine_metadata,
                                           output_dir = output_dir,
                                           filetag = "modified_ait_compoundclass_",
                                           size_axis = 0.05,
                                           dis = "ait")




########################################################

#    boxplot significance comparison                   #

########################################################

bray_pcaxis <- ec_pcoa[["modified_bray"]]$points[,1:2] %>% as.data.frame()
colnames(bray_pcaxis)[1:2] <- c("axis1", "axis2")

bray_pcaxis <- bray_pcaxis %>% 
        rownames_to_column(., var = "internal_sample_id") %>% 
        left_join(., combine_metadata %>% dplyr::select(internal_sample_id, continent))

library(rstatix)
bray_pcaxis_statistic_axis1 <- bray_pcaxis %>%
        pivot_longer(axis1:axis2) %>% 
        group_by(name) %>%  
        rstatix::wilcox_test(value ~ continent) %>% 
        adjust_pvalue(method = "fdr")

write_tsv(bray_pcaxis_statistic_axis1, str_c(output_dir, "/axis_bray_significance.tsv"))

aitchison_pcaxis <- ec_pcoa[["modified_ait"]]$points[,1:2] %>% as.data.frame()
colnames(aitchison_pcaxis)[1:2] <- c("axis1", "axis2")

aitchison_pcaxis <- aitchison_pcaxis %>%
        rownames_to_column(., var = "internal_sample_id") %>%
        left_join(., combine_metadata %>% dplyr::select(internal_sample_id, continent))

library(rstatix)
aitchison_pcaxis_statistic_axis1 <- aitchison_pcaxis %>%
        pivot_longer(axis1:axis2) %>%
        group_by(name) %>%
        rstatix::wilcox_test(value ~ continent) %>%
        adjust_pvalue(method = "fdr")

write_tsv(aitchison_pcaxis_statistic_axis1, str_c(output_dir, "/axis_aitchison_significance.tsv"))



save.image(str_c(output_dir, "step16_manual_check_adonis_env_compoundclass.RData"))



