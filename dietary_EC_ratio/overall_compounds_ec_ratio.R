# the aim of this script is to check dietary-ec ratio for all phytonutrients
# author: lu.zhang
# date: 16,Oct,2023 


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Rscript script.R modified_ec_profile ec_profile metatable output_dir", call.=FALSE)
}


############################################################

# loading package                                          #

############################################################

library(tidyverse)


############################################################

# input file 					

############################################################


# modified ec profile long format 
modified_ec_profile <- args[1] %>% readRDS(.)

# ec profile 
ec_profile_dir <- args[2]

# metatable 
metadata <- read_tsv(args[3], guess_max = 3000)

# output_dir 
output_dir <- args[4]

############################################################

# calculate the total number of EC 

############################################################

# read files ec unstratified file from the directory,those ec profile has been filtered the ones withrowsum==0 file 
ec_profile <- list()

files <- list.files(ec_profile_dir, pattern = "rds")

for (file in files){
	    continent <- str_remove(file, "all_ec_profile") %>% str_remove(".rds")
    print(continent)
        ec_profile[[continent]] <- readRDS(paste(ec_profile_dir, file, sep = "/"))
}


sample_number <- 0
ec_profile_sel <- list()
for (i in names(ec_profile)){
	    print(i)
    ec_profile_sel[[i]] <- ec_profile[[i]]
    sample_number <- sample_number + ncol(ec_profile[[i]])-1
}


summary_all_ec_percentage <- list()
for (i in names(ec_profile)){
	summary_all_ec_percentage[[i]] <- tibble(sample_name = rep("NA", ncol(ec_profile[[i]]) - 1),
						 all_ec_number = rep(0, ncol(ec_profile[[i]]) - 1),
						 continent = "NA")
        sample_names <- colnames(ec_profile_sel[[i]])[2:(ncol(ec_profile_sel[[i]]))]
        for (col in 2:(ncol(ec_profile_sel[[i]]))){
		sample_all_ecs <- ec_profile_sel[[i]][,col][ec_profile_sel[[i]][,col] > 0]
		summary_all_ec_percentage[[i]][(col-1),] <- tibble(sample_name = sample_names[col-1], 
		                                                   all_ec_number = length(sample_all_ecs), 																	 continent = i)
	}
}

summary_all_ec_percentage_tibble <- do.call(rbind, summary_all_ec_percentage)

summary_modified_ec_tibble <- tibble(sample_name = colnames(modified_ec_profile)[2:3069], modified_ec_number = rep(0, 3068))
for (k in colnames(modified_ec_profile)[2:3069]){
	print(k)
	sample_modified_ecs <- modified_ec_profile %>% select(ec_simplied, all_of(k)) 
	colnames(sample_modified_ecs)[2] <- "value"
	sample_modified_ecs_sel <- sample_modified_ecs %>% filter(value > 0)
	modified_ec_number <- sample_modified_ecs_sel$ec_simplied %>% unique() %>% length()
	summary_modified_ec_tibble$modified_ec_number[summary_modified_ec_tibble$sample_name == k] <- modified_ec_number

}

summary_all_ec_percentage_modified_ec_tibble <- summary_all_ec_percentage_tibble %>% 
	left_join(., summary_modified_ec_tibble) %>% 
	mutate(percentage = modified_ec_number/all_ec_number, group = "Microbiome Samples")

library(ggplot2)
plot_violin_v3 <- ggplot(summary_all_ec_percentage_modified_ec_tibble, aes(x = group, y = percentage)) +
	ggdist::stat_halfeye(adjust = .5,
		             width = .6,
		             .width = 0,
			     justification = -.2,
			     point_colour = "#b5dae6") +
	 	geom_boxplot(
		             width = .15,
			     outlier.shape = NA,
			     col = "#b5dae6") + # draw jitter points 
		geom_point(size = 1.3,
		     	   alpha = .3,
		           col = "#b5dae6",
		           position = position_jitter(seed = 1, width = .1)) +
	        coord_cartesian(xlim = c(1.2, NA), clip = "off") +
	        scale_y_continuous(name = "Ratio", breaks = c(0.66,0.68,0.70,0.72,0.74,0.76,0.78), limits = c(0.64,0.80)) +
	        xlab("") +
	        theme_bw() +
	        annotate("label", x = 0.7, y = 0.66,
			 label = paste("Mean(sd):", "\n",round(mean(summary_all_ec_percentage_modified_ec_tibble$percentage),3)," \u00B1 ", 
				       round(sd(summary_all_ec_percentage_modified_ec_tibble$percentage),3),sep = ""),
			 color = "black", size = 8) +
  		theme(text = element_text(size=30),
	        axis.line = element_line(colour = "black"),
		    panel.grid.minor = element_blank(),
		    panel.border = element_blank(),
		    panel.background = element_blank())




pdf(paste(output_dir, "distribution_plot.pdf", sep = "/"))
print(plot_violin_v3)
dev.off()

save.image(str_c(output_dir, "all_ec_distribution_plot.RData", sep = "/"))



