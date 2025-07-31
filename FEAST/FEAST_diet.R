library(FEAST)
library(tidyverse)
library(readxl)


setwd("/media/amarfil/data/Documents/SmallMolecules/Thai_American/FEAST/")

data <- read.delim("formatted_foods.tsv") %>% column_to_rownames(., "Sample")
metadata <- read_excel("metadata.xlsx") %>% column_to_rownames(., "SampleID") %>% varhandle::unfactor(.)

data[] <- apply(data, 2, function(x) as.integer(unlist(x))) %>% varhandle::unfactor(.) 

data <- as.matrix(data)

# feast -------------------------------------------------------------------

set.seed(1)
FEAST_output <- FEAST(C = data, 
                      metadata = metadata, 
                      different_sources_flag = 0, 
                      dir_path = "/media/amarfil/data/Documents/SmallMolecules/Thai_American/FEAST/", 
                      outfile = "foods")


# plot --------------------------------------------------------------------

Thai <- 0
First <- 0
LTR <- 0

for(i in 1:ncol(FEAST_output)) {
  if (grepl("Thai", colnames(FEAST_output[i]))) {
    Thai <- Thai + FEAST_output[, i]
  }
  if (grepl("First", colnames(FEAST_output[i]))) {
    First <- First + FEAST_output[, i]
  }
  if (grepl("LTR", colnames(FEAST_output[i]))) {
    LTR <- LTR + FEAST_output[, i]
  }
}

sources <- data.frame(Sample = rep(rownames(FEAST_output), 4), 
                      Source = factor(c(rep("Thai", 15), rep("First", 15), rep("LTR", 15), rep("Unknown", 15)), levels =c("Thai", "First", "LTR", "Unknown") ), 
                      Value = c(Thai, First, LTR, FEAST_output$Unknown)*100)

ggplot(sources, aes(x = Source, y = Value)) +
  geom_jitter(position = position_jitter(width = .2, height = 0), aes(colour = Source)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  stat_compare_means(comparisons = list(c("Thai", "First"), c("First", "LTR"), c("Thai", "LTR")), tip.length = 0, label.y = c(30, 40, 50)) +
  ggtitle("European-American as sink") +
  scale_x_discrete(labels = c("Thai" = "Thailand", "First" = "New\narrivals", "LTR" = "Long-term\nresidents")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  theme_test() +
  ylab("% Contribution") +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none")