library(FEAST)
library(tidyverse)
library(readxl)
library(ggpubr)
library(varhandle)


setwd("E:/backup/SmallMolecules/Thai_American/FEAST/")

data <- read.delim("../Thai_ECs_unstratified.tsv") 
metadata <- read_excel("diet/metadata.xlsx") %>% column_to_rownames(., "SampleID") %>% varhandle::unfactor(.)
sample_names <- read.delim("../filereport_read_run_PRJEB28687_tsv.txt") %>% 
  select("run_accession", "sample_alias")
dietary_ECs <- read.delim("../../new_tree/new_NutriChem.tsv")$ECs %>% 
  unfactor() %>% 
  strsplit(", ") %>% 
  unlist() %>% 
  unique()

data <- column_to_rownames(data, "X..Gene.Family") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("run_accession") %>% 
  mutate(run_accession = gsub("_Abundance.RPKs", "", run_accession)) %>% 
  merge(., sample_names, by = "run_accession") %>% 
  select(-run_accession) %>% 
  mutate(sample_alias = gsub("wgs.", "", sample_alias)) %>% 
  mutate(sample_alias = gsub("_", ".", sample_alias)) %>% 
  column_to_rownames("sample_alias") %>% 
  t() %>% 
  as.data.frame() %>% 
  subset(rownames(.) %in% dietary_ECs) %>% 
  select(rownames(metadata)) %>% 
  t() %>% 
  as.data.frame()

data[] <- apply(data, 2, function(x) as.integer(unlist(x))) %>% 
  varhandle::unfactor(.) 

data <- as.matrix(data)

# feast -------------------------------------------------------------------

set.seed(1)
FEAST_output <- FEAST(C = data, 
                      metadata = metadata, 
                      different_sources_flag = 0, 
                      dir_path = "dietaryECs/", 
                      outfile = "foods")


# plot --------------------------------------------------------------------

Thai <- 0
First <- 0
LTR <- 0

FEAST_output <- read.delim("foods_source_contributions_matrix.txt")

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
                      Source = factor(c(rep("Thai", 15), rep("First", 15), rep("LTR", 15), rep("Unknown", 15)), levels =c("Thai", "First", "LTR", "Unknown")), 
                      Value = c(Thai, First, LTR, FEAST_output$Unknown)*100)

ggplot(sources, aes(x = Source, y = Value)) +
  geom_jitter(position = position_jitter(width = .2, height = 0), aes(colour = Source)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  stat_compare_means(comparisons = list(c("Thai", "First"), c("First", "LTR"), c("Thai", "LTR")), tip.length = 0) +
  ggtitle("European-American as sink") +
  scale_x_discrete(labels = c("Thai" = "Thailand", "First" = "New\narrivals", "LTR" = "Long-term\nresidents")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  theme_test() +
  ylab("% Contribution") +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none") +
  ylim(0, 100)
