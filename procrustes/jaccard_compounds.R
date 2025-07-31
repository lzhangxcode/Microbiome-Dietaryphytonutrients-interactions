setwd("E:/backup/SmallMolecules/Daily_dataset/")

library(tidyverse)
library(vegan)


# load data ---------------------------------------------------------------

read_multiple_files <- function(path) {
  files <- list.files(path, recursive = F, include.dirs = F)
  files <- files %>% 
    subset(grepl("^compounds_nutrichem_decay", .)) %>% 
    subset(!grepl("RData", .))
  my_list <- list()
  for (i in 1:length(files)) {
    file <- paste(files[i])
    file_name <- file
    file <- read.delim(paste(path, files[i], sep = ""))
    my_list[[file_name]] <- file
  }
  return(my_list)
}

compound_tables <- read_multiple_files(path = "food_records/decay/nutrichem/")
load("Mantel/samples_metagenomics.RData")

# remove patients 11 and 12 -----------------------------------------------

compound_tables$compounds_nutrichem_decay_MCTs11.tsv <- NULL
compound_tables$compounds_nutrichem_decay_MCTs12.tsv <- NULL

# merge tables ------------------------------------------------------------

table_merged <- compound_tables[[1]]

for(i in 2:length(names(compound_tables))) {
  table_merged <- merge(table_merged, compound_tables[[i]], by = "Compounds")
}

table_merged <- column_to_rownames(table_merged, "Compounds")


# remove compounds and samples with more than 90% zeros -------------------

# only keep samples with metagenomics and diet records
table_merged <- select_if(table_merged, names(table_merged) %in% samples_metagenomics)


# compute bray curtis distance --------------------------------------------

jaccard <- as.matrix(vegdist(t(table_merged), "jaccard", binary = T))

jaccard <- as.data.frame(jaccard) %>% 
  rownames_to_column("Sample")

# write.table(jaccard, "procrustes/all_subjects/all_samples/Jaccard_Compounds.tsv", quote = F, sep = "\t", row.names = F)
