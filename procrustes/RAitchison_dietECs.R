setwd("E:/backup/SmallMolecules/Daily_dataset/ECs/imputed_0s/")

library(tidyverse)
library(robCompositions)
library(metagMisc)
library(ape)
library(vegan)

# load data ---------------------------------------------------------------

read_multiple_files <- function(path) {
  files <- list.files(path)
  my_list <- list()
  for (i in 1:length(files)) {
    file <- paste(files[i])
    file_name <- file
    file <- read.delim(files[i], sep = " ") %>% 
      column_to_rownames(., "Sample") %>%
      t(.) %>% 
      as.data.frame(.)
    rownames(file) <- gsub(pattern = "X", replacement = "", rownames(file))
    my_list[[file_name]] <- file
  }
  return(my_list)
}

files <- read_multiple_files(path = ".")


# filter dietary ECs ------------------------------------------------------

compounds_records <- read.delim("../../food_records/decay/merged_decay_foods_nutrichem.csv")$TAXID

dietary_ECs <- read.delim("../../../new_tree/new_NutriChem.tsv") %>% 
  subset(PlantID %in% compounds_records)

dietary_ECs <- dietary_ECs$ECs %>% 
  varhandle::unfactor() %>% 
  strsplit(", ") %>% 
  unlist() %>% 
  unique()

for (i in 1:length(names(files))) {
  files[[i]] <- subset(files[[i]], rownames(files[[i]]) %in% dietary_ECs)
}

merged_table <- files[[1]] %>% 
  rownames_to_column("ECs")

for(i in 2:length(names(files))) {
  sub_files <- files[[i]] %>% 
    rownames_to_column("ECs")
  
  merged_table <- merge(merged_table, sub_files, by = "ECs")
}

merged_table <- column_to_rownames(merged_table, "ECs")
raitchison <- as.matrix(vegdist(t(merged_table), "robust.aitchison"))

raitchison <- as.data.frame(raitchison) %>% 
  rownames_to_column("Sample")

# write.table(raitchison, "../../procrustes/all_subjects/all_samples/RAitchison_dietECs_updated.tsv", quote = F, sep = "\t", row.names = F)
