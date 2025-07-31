setwd("E:/backup/SmallMolecules/Daily_dataset/")

library(vegan)
library(tidyverse)
library(ape)


# load data ---------------------------------------------------------------

compounds <- read.delim("procrustes/all_subjects/all_samples/Jaccard_Compounds.tsv")
dietECs <- read.delim("procrustes/all_subjects/all_samples/RAitchison_dietECs_updated.tsv")


# make PCoAs --------------------------------------------------------------

compounds <- column_to_rownames(compounds, "Sample")
dietECs <- column_to_rownames(dietECs, "Sample")

compounds[] <- apply(compounds, 2, function(x) ifelse(is.na(x), 0, x))
dietECs[] <- apply(dietECs, 2, function(x) ifelse(is.na(x), 0, x))

pcoa_compounds <- pcoa(compounds)
pcoa_ecs <- pcoa(dietECs)

# procrustes analyses -----------------------------------------------------

pcoa_f <- pcoa_compounds$vectors # compounds pcoa
pcoa_t <- pcoa_ecs$vectors # ECs pcoa

# make sure the samples are in the same order
pcoa_f <- subset(pcoa_f, rownames(pcoa_f) %in% rownames(pcoa_t))
pcoa_t <- subset(pcoa_t, rownames(pcoa_t) %in% rownames(pcoa_f))
pcoa_f <- pcoa_f[rownames(pcoa_t), ]

# make sure that they both have the same number of axis
col_axes <- intersect(colnames(pcoa_f), colnames(pcoa_t))
pcoa_f <- dplyr::select(as.data.frame(pcoa_f), col_axes)
pcoa_t <- dplyr::select(as.data.frame(pcoa_t), col_axes)

# procrustes
set.seed(1)
pro <- procrustes(pcoa_f, pcoa_t) # we transform the dependant variable (ECs)
pval <- protest(pcoa_f, pcoa_t, perm = 999)$signif # test for significance

beta_pro <- data.frame(pro$X) # vectors of compounds
trans_pro <- data.frame(pro$Yrot) # new vectors of ECs

# add subject ID and data type to tables (for plotting purposes)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Phytonutrients (Jaccard's)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Phytonutrient-associated enzymes (Aitchison's)"

# merge the 2 tables
colnames(trans_pro) <- colnames(beta_pro)
plot <- rbind(beta_pro, trans_pro)


# plot --------------------------------------------------------------------

ggplot(plot) +
  geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#1285c7", "#bf530f")) +
  theme_classic() +
  geom_line(aes(x = Axis.1, y = Axis.2, group = UserName), col = "darkgrey", alpha = 0.6) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  #ggtitle("Procrustes Analysis\nPlant-diet ECs - Compounds") +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  annotate(geom = "text", x = .2, y = -.3, label = "p = 0.001") 
