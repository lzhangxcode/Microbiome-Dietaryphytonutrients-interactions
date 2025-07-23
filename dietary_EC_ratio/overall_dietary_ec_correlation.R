# the aim of this script is to produce alpha diversity and do correlation with percentage 
# author: lu.zhang
# date: 12rd,Dec,2023 

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
	  stop("Rscript script.R metaphlan_file combine_metadata dietary_percentage_table depth_file output_dir", call.=FALSE)
}


############################################################

# loading package                                          #

############################################################

library(tidyverse)
library(fossil)
library(ppcor)  
library(ggpubr)
library(ggplot2)

############################################################

# input file                                    

############################################################


print("reading input file ... ... ")
# metaphlan file 
print("reading metaphlan file ... ...")
metaphlan_file <- args[1]
metaphlan_file <- read_tsv(metaphlan_file)
colnames(metaphlan_file) <- str_remove_all(colnames(metaphlan_file), "_metaphlan_bowtie2") %>% str_remove_all(., "_T1")

# combine metadata 
print("reading metadata ... ...")
metadata <- read_tsv(args[2], guess_max = 3000)
metaphlan_file_sel <- metaphlan_file %>%
            dplyr::select(., c("clade_name", metadata$internal_sample_id))


# Percentage table 
print("reading percentage table ... ...")
percentage_table <- readRDS(args[3])

# sequencing depth file 
print("reading sequence depth file ... ...")
sequencing_depth <- readRDS(args[4])

# I need to change all manos dataset > 20million as 20million since he did a subsampling for those samples 

sequencing_depth_revised <- sequencing_depth 
sequencing_depth_revised$nonhost[sequencing_depth_revised$study_name == "Manos_Diabesity_Project" & sequencing_depth_revised$nonhost > 20000000] <- 20000000
sequencing_depth_revised$nonhost[sequencing_depth_revised$study_name == "Manos_Prospective_NAFLD_Project" & sequencing_depth_revised$nonhost > 20000000] <- 20000000
sequencing_depth_revised$nonhost[sequencing_depth_revised$study_name == "Manos_Atherosclerosis_Project" & sequencing_depth_revised$nonhost > 20000000] <- 20000000
sequencing_depth_revised$nonhost[sequencing_depth_revised$study_name == "Manos_Hypertension_Project" & sequencing_depth_revised$nonhost > 20000000] <- 20000000

print("finished revised depth file......")

# output dir 
print("reading output dir")
output_dir <- args[5]
print(str_c("the file will be written into this dir:", output_dir))
print("...................................................................................")

saveRDS(sequencing_depth_revised, str_c(output_dir, "/sequencing_depth_revised.rds"))
write_tsv(sequencing_depth_revised, str_c(output_dir, "/sequencing_depth_revised.tsv"))


############################################################

# calculate the alpha diversity                                     

############################################################
print("calculating alpha diversity ... ... ")
# function 
calculate_alpha_diversity <- function(ec_profile_input = NA, continent = NA){

    ec_profile_input <- ec_profile_input %>%
        column_to_rownames(var = "clade_name") 
    # richness 
    input_richness <- apply(ec_profile_input>0, 2, sum) %>% 
        data.frame(richness = .) %>% 
        rownames_to_column(var = "samples") %>%
        as_tibble()

    # simpson : https://github.com/vegandevs/vegan/issues/419 
    # here is gini simpson index 1-simpson index 
    input_simpson <- vegan::diversity(t(ec_profile_input), "simpson") 
    # shannon 
    input_shannon_H <- vegan::diversity(t(ec_profile_input), "shannon") 
    # Matrices are by default treated such that each row is a different taxon and each column is a sample or locality
    input_ACE <- tibble(samples = colnames(ec_profile_input), ACE = 0)
    for (col in 1:ncol(ec_profile_input)){
        ACE_index <- fossil::ACE(ec_profile_input[,col], taxa.row = TRUE) 
        input_ACE[col,"ACE"] <- ACE_index
    }
    
    # summary 
    summary_table <- tibble(samples = input_richness$samples,
                            Simpson = input_simpson,
                            Shannon = input_shannon_H,
                            continent = continent)

    return(summary_table)

}


Species_alpha_diversity <- calculate_alpha_diversity(ec_profile_input = metaphlan_file_sel, continent = "no")

write_delim(Species_alpha_diversity, paste(output_dir, "alpha_diversity_species_no_rarefy.tsv", sep = "/"))
write_rds(Species_alpha_diversity, paste(output_dir, "alpha_diversity_species_no_rarefy.rds", sep = "/"))


############################################################

# calculate the correlation                                     

############################################################

# spearman correlation 
spearman_correlation <- function(input_summary_table, percentage_table){

    statistics <- list()
    join_table <- left_join(input_summary_table, percentage_table, by = c("samples" = "sample_name"))
    spearman_correlation_of_percentage_Simpson <- cor.test(join_table$Simpson, join_table$percentage, method = "spearman")
    spearman_correlation_of_percentage_Shannon <- cor.test(join_table$Shannon, join_table$percentage, method = "spearman")

    statistics[["join_table"]] <- join_table
    statistics[["Simpson"]] <- spearman_correlation_of_percentage_Simpson
    statistics[["Shannon"]] <- spearman_correlation_of_percentage_Shannon

    return(statistics)

}




species_spearman_statics <- spearman_correlation(input_summary_table = Species_alpha_diversity, percentage_table = percentage_table)

saveRDS(species_spearman_statics, str_c(output_dir, "/spearman_correlation.rds"))


############################################################

# calculate the correlation (partial spearman)                                    

############################################################


# partial spearman correlation 


alpha_table_with_metadata <- Species_alpha_diversity %>% left_join(., metadata, by = c("samples" = "internal_sample_id")) %>% 
  left_join(., sequencing_depth, by = c("samples" = "internal_sample_id")) %>% 
  left_join(., percentage_table, by = c("samples" = "sample_name"))

colnames(alpha_table_with_metadata)[8] <- "subject_id"
colnames(alpha_table_with_metadata)[6] <- "study_name"

alpha_table_with_metadata$nonhost_log <- log2(alpha_table_with_metadata$nonhost)

alpha_table_with_metadata$sample_id.x[is.na(alpha_table_with_metadata$subject_id)]

shannon_partial <- pcor.test(alpha_table_with_metadata$percentage, alpha_table_with_metadata$Shannon, alpha_table_with_metadata$nonhost_log, method="spearman")
simpson_partial <- pcor.test(alpha_table_with_metadata$percentage, alpha_table_with_metadata$Simpson, alpha_table_with_metadata$nonhost_log, method="spearman")

summary_species_partial_correlation <- tibble(shannon_partial_r = shannon_partial$estimate,
                                              simpson_partial_r = simpson_partial$estimate,
                                              shannon_partial_p = shannon_partial$p.value,
                                              simpson_partial_p = simpson_partial$p.value)


############################################################

# draw the partial spearman correlation plot                                 

############################################################

r <- c(summary_species_partial_correlation$shannon_partial_r %>% round(., digits = 4),
       summary_species_partial_correlation$simpson_partial_r  %>% round(., digits = 4),
       summary_species_partial_correlation$inverse_simpson_partial_r %>% round(., digits = 4))

p_values <- c(summary_species_partial_correlation$shannon_partial_p,
	      summary_species_partial_correlation$simpson_partial_p,
	      summary_species_partial_correlation$inverse_simpson_partial_p)

p_values <- scales::scientific(p_values, digits = 3)

index <- c("Shannon", "Simpson")


p <- list()
for (x in c("Shannon", "Simpson")){

        r_target <- r[which(index == x)]
        p_target <- p_values[which(index == x)]
        label <- paste0("Parital Spearman, ","r", " = ", r_target, ", ", "p = ", p_target)
        label <- paste0("R", " = ", r_target, ", ", "p = ", p_target)

        p[[x]] <- ggscatter(alpha_table_with_metadata, x = x, y = "percentage",
          add = "reg.line", conf.int = TRUE, add.params = list(color = "red"),
          cor.coef = FALSE, cor.method = "spearman", cor.coef.size = 8,
          xlab = paste0(x, "(", "species", ")"), ylab = "N of plant-diet associated ECs/\nN of total microbiome ECs") +
          annotate("text", x = 0.55*max(alpha_table_with_metadata[,x]), y = 0.98*max(alpha_table_with_metadata[,"percentage"]),
                       label = label,
                       color = "black", size = 6) +
          theme(text = element_text(size=30),
                axis.line = element_line(colour = "black"),
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
		plot.margin = margin(1,1,1.5,1.2, "cm"))
}

print("ploting ...")
pdf(paste(output_dir, "partial_spearman_non_rarefy_percentage_alphadiversity_species_shannon.pdf", sep = "/"), width = 10*1.2, height = 6*1.2)
p[["Shannon"]]
dev.off()


pdf(paste(output_dir, "partial_spearman_non_rarefy_percentage_alphadiversity_species_simpson.pdf", sep = "/"), width = 10*1.2, height = 6*1.2)
p[["Simpson"]]
dev.off()

print("saving RData ...")
save.image(paste(output_dir, "non_rarefy_percentage_step14_subscript4_manual_check_partial.RData", sep = "/"))



