# this script is used for draw the piechart for shared compounds and unique compounds between probiotic and gut microbiota 
# lu.zhang
# 2023.11.25 
############################################################################

# prepare input 

############################################################################

library(tidyverse)

# gut microbiota total after manual check compounds 
modified_ec_profile_sep <- readRDS("modified_ec_profile_after_np.rds") 

# probiotic compounds 
probiotic_species_link <- "probiotic_compounds_ec_links.tsv" %>% read_tsv(.)

# gut microbiota secondary compounds 
nutrichem_mapped_compounds_with_ec <- read_tsv("nutrichem_mapped_compounds_with_ec.tsv")

modified_ec_profile_combine_simple_sel_compound <- modified_ec_profile_sep %>% 
	filter(Compounds %in% nutrichem_mapped_compounds_with_ec$CompoundID) 

# probiotic secondary compounds 
probiotic_compounds_ec_sep <- probiotic_species_link %>% 
        mutate(species = strain) 

probiotic_compounds_ec_sep_comp <- probiotic_compounds_ec_sep %>% filter(Compounds %in% nutrichem_mapped_compounds_with_ec$CompoundID) 


# output_dir 


output_dir <- "piechart_probiotic_gutmirobiota/"


#########################################################################################


#                       summary of the number of compounds                              # 


#########################################################################################

# all compounds number 

modified_ec_profile_sep$Compounds %>% unique() %>% length() #775 compounds 
probiotic_compounds_ec_sep$Compounds %>% unique() %>% length() 

all_compounds_gut <- unique(modified_ec_profile_sep$Compounds)
all_compounds_pro <- unique(probiotic_compounds_ec_sep$Compounds)

# secondary compounds 

secondary_gut <- modified_ec_profile_combine_simple_sel_compound$Compounds %>% unique() #186 compounds 
secondary_pro <- unique(probiotic_compounds_ec_sep_comp$Compounds)


#########################################################################################


#                       piechart of the compounds                                       # 


#########################################################################################

# reference tutorial:  https://r-charts.com/part-whole/pie3d/


# install.packages("plotrix")
library(plotrix)
library(ggsci)

a <- intersect(all_compounds_gut, all_compounds_pro) %>% length()
b <- setdiff(all_compounds_gut, all_compounds_pro) %>% length()

all_input <- c(a, b)
mypal <- pal_npg("nrc", alpha = 1)(9) %>% .[7:8]
pdf(paste(output_dir, "all_compounds_comparion_3d_pie.pdf", sep = "/"))
pie3D(all_input, col = mypal, border = "black", labels = all_input, labelcol = "black", labelcex = 1, explode = 0.1, main = "All Dietary Compounds")
par(xpd=TRUE)
legend(0.1, 0.8,legend=c("common chemical space", "unique chemical space"), cex=0.8, yjust=0.2, xjust = 0.1,
       fill = mypal) 
dev.off()


mypal2 <- pal_npg("nrc", alpha = 1)(9) %>% .[5:6]
c <- intersect(secondary_pro, secondary_gut) %>% length()
d <- setdiff(secondary_gut, secondary_pro) %>% length()
all_input_second <- c(c, d)
pdf(paste(output_dir, "secondary_compounds_comparion_3d_pie.pdf", sep = "/"))
pie3D(all_input_second, col = mypal, border = "black", labels = all_input_second, labelcol = "black", labelcex = 1, explode = 0.1, main = "Dietary Compounds Related to \n Secondary Metabolism")
par(xpd=TRUE)
legend(0.1, 0.8,legend=c("common chemical space", "unique chemical space"), cex=0.8, yjust=0.2, xjust = 0.1,
       fill = mypal)
dev.off()


save.image(str_c(output_dir, "/piechart.RData"))



