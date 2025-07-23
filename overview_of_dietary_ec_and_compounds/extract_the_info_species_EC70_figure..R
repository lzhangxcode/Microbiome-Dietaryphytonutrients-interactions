# draw figure for the heatmap for compounds-sp 
# read the files 
setwd("./minisp_70ECs")
sp_abundance_prevalence_sel_0.05_wide <- readRDS("sp_abundance_prevalence_sel_0.05_wide.rds")
rowannotation <- readRDS("rowannotation.rds") %>% rename(Phylum = phylum)
columnannotateion <- readRDS("columnannotateion.rds")

library(tidyverse)
compounds_order <- sp_abundance_prevalence_sel_0.05_wide %>% as.data.frame(.) %>% 
  rownames_to_column(., var = "species") %>% 
  pivot_longer(!species) %>% group_by(name) %>% 
  summarise(num_of_sp = sum(value)) %>% 
  left_join(., columnannotateion %>% rownames_to_column(., var = "name")) %>% 
  arrange(group, desc(num_of_sp)) 

sp_abundance_prevalence_sel_0.05_wide_order <- sp_abundance_prevalence_sel_0.05_wide %>% select(compounds_order$name)

columnannotateion2_order <- data.frame(Group = columnannotateion[colnames(sp_abundance_prevalence_sel_0.05_wide_order),])
rownames(columnannotateion2_order) <- colnames(sp_abundance_prevalence_sel_0.05_wide_order)

library(ComplexHeatmap)
columnannotateion2_order$Group[columnannotateion2_order$Group == "overlap"] <- "common"
Top_anno_order = HeatmapAnnotation(df = columnannotateion2_order, col = list(Group = c("common" = "#91d1c1", "unique" = "#dd0709")))

rowannotationv2 <- rowannotation %>% mutate(Phylum = str_remove(Phylum, "p__"))

ha <- rowAnnotation(df = rowannotationv2, 
                    col = list(Phylum = c("Actinobacteria" = "#619cff", 
                                          "Bacteroidetes" = "#f8756b", # 
                                          "Firmicutes" = "#4cb030", 
                                          "Lentisphaerae" = "#00bad7", 
                                          "Proteobacteria" = "#00bfae", 
                                          "Synergistetes" = "#e668f3", 
                                          "Verrucomicrobia" = "#ca78a6")))



# compound1: 
#Fustin
label2 <- sp_abundance_prevalence_sel_0.05_wide_order["Fustin"]
colnames(label2)[1] <- 'label2'
label2 %>% filter(label2 > 0)
# I will label : g__Lactococcus.s__Lactococcus_lactis
which(rownames(sp_abundance_prevalence_sel_0.05_wide_order) == "g__Lactococcus.s__Lactococcus_lactis")

# compound2 
#Caffeines
label3 <- sp_abundance_prevalence_sel_0.05_wide_order["Caffeines"]
colnames(label3)[1] <- 'label3'
label3 %>% filter(label3 > 0)

# choose g__Lachnoclostridium.s__Clostridium_citroniae
which(rownames(sp_abundance_prevalence_sel_0.05_wide_order) == "g__Lachnoclostridium.s__Clostridium_citroniae")


# compound3
#Butein 
label4 <- sp_abundance_prevalence_sel_0.05_wide_order["Butein"]
colnames(label4)[1] <- 'label4'
label4 %>% filter(label4 > 0)

# choose g__Eubacterium.s__Eubacterium_ramulus

# add extra compounds 
# compound4: 
#Tricetin 
#CDBNO:2571;CHEBI:507499;CHEBI:60045;CHEMLIST:4267852;CID:5281701

label5 <- sp_abundance_prevalence_sel_0.05_wide_order["CDBNO:2571;CHEBI:507499;CHEBI:60045;CHEMLIST:4267852;CID:5281701"]
colnames(label5)[1] <- 'label5'
label5_sp <- label5 %>% filter(label5 > 0)
# I will label : g__Bifidobacterium.s__Bifidobacterium_animalis
which(rownames(sp_abundance_prevalence_sel_0.05_wide_order) == "g__Bifidobacterium.s__Bifidobacterium_animalis") #13 

nutrichem <- read_tsv("./compound_associations.tsv")

metaphlan_file <- read_tsv("./merge_Cohort_metaphlan_species.tsv")


ha2_rename = rowAnnotation(foo = anno_mark(at = c(13,73,65,193), 
                                    labels = c("Bifidobacterium_animalis", "Lactococcus_lactis", "Clostridium_citroniae", "Eubacterium_ramulus")))
ha2_rename = rowAnnotation(foo = anno_mark(at = c(13,73,65,193), 
                                           labels = c("Bifidobacterium animalis", "Lactococcus lactis", "Clostridium citroniae", "Eubacterium ramulus")))

library(ComplexHeatmap)
pdf("./species_compounds_uncalssified_labelv3.pdf", width = 10)
Heatmap(sp_abundance_prevalence_sel_0.05_wide_order, 
        name = "A", col = c("1" = "Forest green", "0" = "#f5f5f5"),  
        heatmap_legend_param = list(at = c(1, 0), labels = c("Present",  "Absent")),  # Forest green
        show_row_names = F, show_column_names = F, top_annotation = Top_anno_order, 
        right_annotation = c(ha, ha2_rename), cluster_columns = F, border = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.1)) # border of the cells 
dev.off()



