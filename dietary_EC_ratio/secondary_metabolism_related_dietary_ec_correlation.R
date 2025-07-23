# the aim of the script is to do parital correlation between the secondary metabolism percentage and species alpha diversity 
# the secondary metabolism percentage is based on total secondary metabolism related EC 
# 2024.09.19 
# lu

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
          stop("Rscript step14_subscript5_nonrarefy_ec_percentage_secondary_metabolism_manual_check.R alpha_diversity_file secondary_metabolism_ec_percentage_table_RData depth_file output_dir", call.=FALSE)
}


############################################################

# loading package                                          #

############################################################
library(tidyverse)
library(ggplot2)



############################################################

# input file

############################################################


# output
output_dir <- args[4]

# load the previous calculated secondary-metabolism related percentage 
load(str_c(output_dir, "secondary_enzyme_distribution_plot_total_secondary_ec.RData"))

# save the secondary percentage object into rds 
saveRDS(secondary_percentage_per_sample_enzyme_simple, str_c(output_dir, "secondary_percentage_per_sample_enzyme_simple.rds"))

# clean all other variables 
rm(list=ls()[! ls() %in% c("args")])

# ec percentage file 
print(args[2])
ec_percentage <- readRDS(args[2])
print("read ec percentage file ...")

# alpha diversity
Species_alpha_diversity <- readRDS(args[1])

# depth file 
depth_file <- readRDS(args[3])

# output dir 
output_dir <- args[4]


############################################################

# partial spearman correlation 

############################################################

library(ppcor) # version 1.1 

alpha_table_with_metadata <- Species_alpha_diversity %>% 
  left_join(., depth_file, by = c("samples" = "internal_sample_id")) %>% 
  left_join(., ec_percentage, by = c("samples" = "sample_name"))

alpha_table_with_metadata$nonhost_log <- log2(alpha_table_with_metadata$nonhost)
alpha_table_with_metadata$sample_id[is.na(alpha_table_with_metadata$subject_id)]
alpha_table_with_metadata[alpha_table_with_metadata$sample_id == "SID31232",]$subject_id <- "SID31232"
alpha_table_with_metadata[alpha_table_with_metadata$sample_id == "N104A",]$subject_id <- "N104A"
alpha_table_with_metadata[alpha_table_with_metadata$sample_id == "N075A",]$subject_id <- "N075A"
# depth file updated, metatable has already updated in very beginning, even though i chagne it here, it has no effect 

shannon_partial <- pcor.test(alpha_table_with_metadata$percentage, alpha_table_with_metadata$Shannon, alpha_table_with_metadata$nonhost_log, method="spearman")
simpson_partial <- pcor.test(alpha_table_with_metadata$percentage, alpha_table_with_metadata$Simpson, alpha_table_with_metadata$nonhost_log, method="spearman")
inverse_simpson_partial <- pcor.test(alpha_table_with_metadata$percentage, alpha_table_with_metadata$inverse_Simpson, alpha_table_with_metadata$nonhost_log, method="spearman")

summary_species_partial_correlation <- tibble(shannon_partial_r = shannon_partial$estimate,
                                              simpson_partial_r = simpson_partial$estimate,
                                              inverse_simpson_partial_r = inverse_simpson_partial$estimate,
                                              shannon_partial_p = shannon_partial$p.value,
                                              simpson_partial_p = simpson_partial$p.value,
                                              inverse_simpson_partial_p = inverse_simpson_partial$p.value)


############################################################

# draw the partial spearman correlation plot                                 

############################################################

library(ggpubr)

# draw plot 

r <- c(summary_species_partial_correlation$shannon_partial_r %>% round(., digits = 4),
       summary_species_partial_correlation$simpson_partial_r  %>% round(., digits = 4),
       summary_species_partial_correlation$inverse_simpson_partial_r %>% round(., digits = 4))

p_values <- c(summary_species_partial_correlation$shannon_partial_p,
              summary_species_partial_correlation$simpson_partial_p,
              summary_species_partial_correlation$inverse_simpson_partial_p)

p_values <- scales::scientific(p_values, digits = 3)

index <- c("Shannon", "Simpson", "inverse_Simpson")


join_table <- alpha_table_with_metadata <- Species_alpha_diversity %>%
  left_join(., ec_percentage, by = c("samples" = "sample_name"))

p <- list()

for (x in c("Simpson", "Shannon")){

        r_target <- r[which(index == x)]
        p_target <- p_values[which(index == x)]
        #label <- paste0("Partial Spearman, ","r", " = ", r_target, ", ", "p = ", p_target)
        label <- paste0("R", " = ", r_target, ", ", "p = ", p_target)

        p[[x]] <- ggscatter(join_table, x = x, y = "percentage",
          add = "reg.line", conf.int = TRUE, add.params = list(color = "red"),
          cor.coef = FALSE, cor.method = "spearman", cor.coef.size = 8,
          xlab = paste0(x, "(species)"), ylab = "N of plant-diet asscociated ECs/\nN of total microbiome ECs") +
          annotate("text", x = 0.55*max(join_table[,x]), y = 0.98*max(join_table[,"percentage"]),
                       label = label,
                       color = "black", size = 6) +
          ggtitle("Secondary Metabolism") + 
          theme(text = element_text(size=30),
                axis.line = element_line(colour = "black"),
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
                plot.margin = margin(1,1,1.5,1.2, "cm")
}



pdf(paste(output_dir, "partial_secondary_spearman_non_rarefy_percentage_alphadiversity_species_shannon_total_secondary_enzyme.pdf", sep = "/"), width = 10*1.1, height = 6*1.2)
p[["Shannon"]] 
dev.off()

pdf(paste(output_dir, "partial_secondary_spearman_non_rarefy_percentage_alphadiversity_species_simpson_total_secondary_enzyme.pdf", sep = "/"), width = 10*1.1, height = 6*1.2)
p[["Simpson"]]
dev.off()

save.image(paste0(output_dir, "/partial_secondary_compared_to_total_secondary.RData"))




