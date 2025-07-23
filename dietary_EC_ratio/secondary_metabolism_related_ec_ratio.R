# the aim of this script is to calculate the secondary-microbiome-dietary-ECs/ total secondary ECs
# author: lu.zhang 
# date: 2024.juli.15 

# the code is run in R

########################################################################

# load the data 

########################################################################
library(tidyverse)

# the RData I use for the EC list is : the variable is : links_table$ec
load("secondary_metabolites_extraction.RData")


# the needed profile 
modified_ec_profile <- readRDS(str_c("./extract_secondary_metablites/", "modified_ec_profile_after_np.rds"))

modified_ec_profile_simple <- modified_ec_profile %>%
        select(ec_simplied, Compounds) %>%
        unique()

# the all ec profile 
all_ec_profile_clean <- readRDS("./othercategory/all_ec_profile_clean.rds")


# the output dir 
output_dir <- "./secondary_ec_percentage/"

########################################################################

# calculate the total secondary EC number 

########################################################################
# head(links_table$ec_simplied) 

colnames(all_ec_profile_clean)[1] <- "ec_simplied"

calculate_per_ec <- function(input = all_ec_profile_clean, input_modified = modified_ec_profile, col1 = NA, col2 = "ec_simplied", compound_list = links_table$ec_simplied){ # links_table$ec_simplied

  input_sel <- input[, c(col1, col2)]
  colnames(input_sel)[1] <- "value"

  input_sel_total <- input_sel %>%
    filter(value > 0) %>% 
    filter(ec_simplied %in% compound_list)

  total_compounds <- input_sel_total$ec_simplied %>% unique() %>% length()

  input_sel_secondary <- input_sel %>%
	  filter(value > 0) %>% 
	  filter(ec_simplied %in% input_modified$ec_simplied) %>% 
          filter(ec_simplied %in% compound_list)

  secondary_compounds <- input_sel_secondary$ec_simplied %>% unique() %>% length()

  percentage <- secondary_compounds/total_compounds

  re <- list(total_compounds, secondary_compounds, percentage)

  return(re)

}

secondary_percentage_per_sample_enzyme <- tibble(sample = colnames(modified_ec_profile)) %>%
  filter(!sample %in% c("ec_simplied", "Compounds", "Tag")) %>%
  rowwise() %>%
  mutate(secondary = list(calculate_per_ec(col1 = sample)))

secondary_percentage_per_sample_enzyme_simple <- secondary_percentage_per_sample_enzyme %>%
  pull(secondary) %>%
  do.call(rbind, .) %>%
  as.data.frame()


secondary_percentage_per_sample_enzyme_simple$total_ec <- secondary_percentage_per_sample_enzyme_simple$V1 %>% as.numeric()
secondary_percentage_per_sample_enzyme_simple$secondary_pathway_ec <- secondary_percentage_per_sample_enzyme_simple$V2 %>% as.numeric()
secondary_percentage_per_sample_enzyme_simple$percentage <- secondary_percentage_per_sample_enzyme_simple$V3 %>% as.numeric()
secondary_percentage_per_sample_enzyme_simple$sample_name <- secondary_percentage_per_sample_enzyme$sample

mean(secondary_percentage_per_sample_enzyme_simple$percentage)
sd(secondary_percentage_per_sample_enzyme_simple$percentage)

# [1] 0.8988744
# [1] 0.01495499
 
################################################################

# plot 

################################################################

ec_percentage <- secondary_percentage_per_sample_enzyme_simple
ec_percentage$samples <- "Microbiome Samples"

plot_violin_enzyme <- ggplot(ec_percentage, aes(x = samples, y = percentage)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .15,
    outlier.shape = NA
  ) + # draw jitter points
  geom_point(
    size = 1.3,
    alpha = .3,
    col = "#b5dae6",
    position = position_jitter(
      seed = 1, width = .1
    )
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  scale_y_continuous(name = "Ratio of ECs involved\nin secondary metabolism") +
  xlab("") +
  theme_bw() +
  annotate("label", x = 0.7, y = round(mean(secondary_percentage_per_sample_enzyme_simple$percentage), 3),
           label = paste("Mean(sd):", "\n",round(mean(secondary_percentage_per_sample_enzyme_simple$percentage),3)," \u00B1 ", round(sd(secondary_percentage_per_sample_enzyme_simple$percentage),3),sep = ""),
           color = "black", size = 6) +
  theme(text = element_text(size=30),
        axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
scale_y_continuous(name = "N of plant-diet associated ECs/\nN of total secondary ECs") + ggtitle("Secondary Metabolism")

saveRDS(plot_violin_enzyme, paste(output_dir, "plot_violin_enzyme_total_secondary_ec.rds", sep = "/"))

pdf(paste(output_dir, "secondary_enzyme_distribution_plot_total_secondary_ec.pdf", sep = "/"))
plot_violin_enzyme
dev.off()

save.image(str_c(output_dir, "secondary_enzyme_distribution_plot_total_secondary_ec.RData"))


