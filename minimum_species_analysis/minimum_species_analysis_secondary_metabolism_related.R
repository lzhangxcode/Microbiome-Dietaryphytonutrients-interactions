# this script is used for get the mini.sp list for secondary-metabolism pathways' related metabolites after prevalence filtering for the species 
# author:lu.zhang
# date: 2023.11.15-11.20 
################################################################################

#                       loading packages                                       #

################################################################################

library(ggplot2)
library(tidyverse) 
library(ggrepel)


########################################################################################################################

#                       input files                                                                                    #

########################################################################################################################

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
  stop("At least four argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) == 8) {
        print("Rscript minimum_species_analysis_secondary_metabolism_related.R dietary_ec_profile et al (total 8 parameters) ...")
}


# 1. the total gut microbiota compounds 
# modified ec profile after np check 
modified_ec_profile_sep <- readRDS(args[1]) 

# 2. secondary compounds list only consider compounds 
nutrichem_mapped_compounds_with_ec <- read_tsv(args[2])

modified_ec_profile_combine_simple_sel_compound <- modified_ec_profile_sep %>% 
        filter(Compounds %in% nutrichem_mapped_compounds_with_ec$CompoundID) #186 secondary compounds in total 

# 3. species and metabolites links (all compounds one) 
gutmicrobiota_species_link <- readRDS(args[3])


# 4. probiotic links file 
probiotic_species_link <- args[4] %>% read_tsv(.)

# 5. metaphlan file 

metaphlan_file <- args[5] %>% read_tsv(.)


# 6. metatable 

metadata <- args[6] %>% read_tsv(., guess_max = 3000)

# 7. database file
database_file <- readRDS(args[7])

# 8. output dir 

output_dir <- args[8]

########################################################################################################################

#                   probiotic secondary metabolites link file : only consider compounds                                #

########################################################################################################################

probiotic_compounds_ec_sep <- probiotic_species_link %>% 
        mutate(species = strain)  

probiotic_compounds_ec_sep_comp <- probiotic_compounds_ec_sep %>% filter(Compounds %in% nutrichem_mapped_compounds_with_ec$CompoundID) 


########################################################################################################################

#                   overlap compounds between probiotic and gut microbiota                                             #

########################################################################################################################

# overlap & secondary compounds 

secondary_compounds_overlap <- intersect(probiotic_compounds_ec_sep_comp$Compounds, modified_ec_profile_combine_simple_sel_compound$Compounds)
secondary_compounds_unique <- setdiff(modified_ec_profile_combine_simple_sel_compound$Compounds, probiotic_compounds_ec_sep_comp$Compounds)
# 116 overlapped compounds and 70 unique compounds 

########################################################################################################################

#                                   Function : # links the EC to compounds by using secondary database file            #

########################################################################################################################

secondary_metabolite_list_modified <- nutrichem_mapped_compounds_with_ec %>% mutate(Compounds = CompoundID, tag = "secondary_links")

# gut microbiota species-ec link file 

extract_species_compounds_links_secondary <- function(species_ec_link = NA, secondary_metabolite_list_modi = NA){

        # left_join them together and separate_rows by using , in compounds column 

        compounds_species_link_second <- species_ec_link %>% 
            left_join(., secondary_metabolite_list_modi, by = c("Compounds")) %>%
            dplyr::filter(tag == "secondary_links") 

        
        return(compounds_species_link_second)
}

secondary_compounds_species_link <- extract_species_compounds_links_secondary(species_ec_link = gutmicrobiota_species_link,
									      secondary_metabolite_list_modi = secondary_metabolite_list_modified)



########################################################################################################################

#                                   Function : # filter only the sp pass 5% prevalence                                 #

########################################################################################################################
colnames(metaphlan_file) <- colnames(metaphlan_file) <- str_remove_all(colnames(metaphlan_file), "_metaphlan_bowtie2") %>% str_remove_all(., "_T1")
metaphlan_file_sel <- metaphlan_file %>%
        select(clade_name, metadata$internal_sample_id)
# pass 5% prevalence bacteria 
metaphlan_file_sel_targetedsp <- metaphlan_file_sel%>%
        mutate(species = str_split(clade_name, pattern = "\\|g__", simplify = T) %>% .[,2]) %>% 
        mutate(species = str_c("g__", species) %>% str_replace("\\|", ".")) %>% 
        mutate(kingdom = str_split(clade_name, pattern = "\\||p__", simplify = T) %>% .[,1]) %>% 
        filter(kingdom == "k__Bacteria")


metaphlan_file_sel_targetedsp_prevalence <- metaphlan_file_sel_targetedsp %>%
        select(-species, -kingdom) %>%
        pivot_longer(!clade_name) %>%
        filter(value > 0) %>%
        group_by(clade_name) %>%
        summarise(n = n()) %>%
        mutate(prevalence = n/3068) %>%
        ungroup() %>%
        mutate(species = str_split(clade_name, pattern = "\\|g__", simplify = T) %>% .[,2]) %>%
        mutate(species = str_c("g__", species) %>% str_replace("\\|", "."))%>%
        arrange(species)
metaphlan_file_sel_targetedsp_prevalence_sel <- metaphlan_file_sel_targetedsp_prevalence %>% filter(prevalence >= 0.05)  

print(str_c("there are ", length(unique(metaphlan_file_sel_targetedsp_prevalence_sel$clade_name)), " species left after 5% prevalence filter in metaphlan file"))

secondary_compounds_species_link_prevalence <- secondary_compounds_species_link %>% 
        filter(species %in% metaphlan_file_sel_targetedsp_prevalence_sel$species) #121 species 

print(str_c("there are ", length(unique(secondary_compounds_species_link_prevalence$species)), " species after all filtering used for secondary metabolites analysis"))
print(str_c("there are ", length(unique(secondary_compounds_species_link_prevalence$Compounds)), " compounds after all filtering used for secondary metabolites analysis"))

########################################################################################################################

#                                   Function : # check how many compounds left                                         #

########################################################################################################################


overlap_compounds_after_contribution <- intersect(secondary_compounds_species_link_prevalence$Compounds, secondary_compounds_overlap) %>%
        intersect(., probiotic_compounds_ec_sep_comp$Compounds)

compounds_species_link_sel <- secondary_compounds_species_link_prevalence %>% filter(Compounds %in% overlap_compounds_after_contribution)
probiotic_compounds_ec_sep_sel <- probiotic_compounds_ec_sep_comp %>% filter(Compounds %in% overlap_compounds_after_contribution)


########################################################################################################################

#                                   Function : # greedy algorithem for mini.sp                                         #

########################################################################################################################

#  greedy to get mini.sp 
greedy_to_extract_minimum_species <- function(combine_sel = NA){


        # I need to produce a matrix for greedy algorithm, first column is set : species, the second column is : compound name
        combine_sel_greedy <- combine_sel %>%
                select(species, Compounds) %>%
                unique() # some duplicated after remove the enzyme information & Tag information

        colnames(combine_sel_greedy) <- c("set", "element")
        set.seed(12345678)
        res <- RcppGreedySetCover::greedySetCover(combine_sel_greedy, FALSE)
        # check idential
        check_greedy_result <- identical(sort(unique(res$element)),sort(unique(combine_sel_greedy$element)))
        if (check_greedy_result == TRUE){
                return(res)

        }else {
                break
        }

    # res$set

}

set.seed(2)
minimum_species_probiotic <- greedy_to_extract_minimum_species(combine_sel = probiotic_compounds_ec_sep_sel) 
# return 4 species 

set.seed(2)
minimum_species_gutmicrobiota <- greedy_to_extract_minimum_species(combine_sel = compounds_species_link_sel)
# return 8 species 

saveRDS(minimum_species_probiotic, paste0(output_dir, "/overlap_compounds_secondary_probiotic_only_comp.rds"))
saveRDS(minimum_species_gutmicrobiota, paste0(output_dir, "/overlap_compounds_secondary_gutmicrobiota_only_comp.rds"))


########################################################################################################################

#                                   Function : # unique compounds mini.sp analysis                                     #

########################################################################################################################


unique_compounds_after_contribution <- intersect(secondary_compounds_species_link_prevalence$Compounds, secondary_compounds_unique) #43 compounds 
 
compounds_species_link_unique_sel <- secondary_compounds_species_link_prevalence %>% filter(Compounds %in% unique_compounds_after_contribution)

set.seed(2)
minimum_species_gutmicrobiota_unique <- greedy_to_extract_minimum_species(combine_sel = compounds_species_link_unique_sel) 
# return 11 species 

saveRDS(minimum_species_gutmicrobiota_unique, paste0(output_dir, "/unique_compounds_secondary_gutmicrobiota_only_comp.rds"))
saveRDS(compounds_species_link_unique_sel, paste0(output_dir, "/compounds_species_link_unique_sel_secondary.rds"))

########################################################################################################################

#                                   Function : # stacked plot for mini.sp                                              #

########################################################################################################################

probiotic_input <- minimum_species_probiotic %>% 
        group_by(set) %>% 
        summarise(compounds = paste(element, collapse = ","), total_num = n()) %>% 
        ungroup() %>% 
        arrange(desc(total_num)) %>% 
        mutate(., compounds_num = cumsum(total_num)) %>% 
        mutate(., percentage = compounds_num/nrow(minimum_species_probiotic), order = seq(1:nrow(.))) %>% 
        mutate(group = "probiotic") %>%
        mutate(contribution = total_num/nrow(minimum_species_gutmicrobiota))



gutmicrobiota_input <- minimum_species_gutmicrobiota %>%           
        group_by(set) %>%
        summarise(compounds = paste(element, collapse = ","), total_num = n()) %>%
        arrange(desc(total_num)) %>%
        mutate(., compounds_num = cumsum(total_num)) %>%
        mutate(., percentage = compounds_num/nrow(minimum_species_gutmicrobiota), order = seq(1:nrow(.))) %>% 
        mutate(group = "gut microbiota") %>% 
        mutate(set = set %>% str_split("s__", simplify =T) %>% .[,2] %>% str_replace_all("_", " ")) %>% 
        mutate(contribution = total_num/nrow(minimum_species_gutmicrobiota))

my36colors <- c('#58A4C3', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C','#585658','#9FA3A8','#E0D4CA', '#E5D2DD', '#C5DEBA', '#5F3D69',   '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')

stack_input <- rbind(gutmicrobiota_input, probiotic_input) %>% mutate(species = ifelse(order > 10, "others", set), new_order =  ifelse(order > 10, 11, order))
saveRDS(stack_input, paste0(output_dir, "/overlap_compounds_secondary_stack_compounds_only_compounds.rds"))

stack_input_modify <- stack_input %>% 
        group_by(species, group, new_order) %>% 
        summarise(total_num = sum(total_num), contribution = sum(contribution)) %>% 
        ungroup() %>% 
        arrange(group,new_order)


stack_input_modify$group[stack_input_modify$group == "probiotic"] <- "probiotics"

mx <- ggplot(stack_input_modify, aes(x = group, fill = species, y = total_num)) +
    geom_bar(stat = "identity", colour = "black", linewidth = 0.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 14, colour = "black", hjust = 1, face= "bold"),
    axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 10, face = "bold.italic", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "", y = "N of Compounds", fill = "Min.Species") +       
    scale_fill_manual(values = my36colors) +
    guides(fill=guide_legend(ncol=1))

ggsave(paste0(output_dir, "/overlap_compounds_secondary_stacked_plot_only_comp.pdf"), mx)


########################################################################################################################

#                                   Function : # minimum curve for unique compounds                                    #

########################################################################################################################


curve_input <- minimum_species_gutmicrobiota_unique
curve_input$num <- 1

curve_input_combine <- curve_input %>%
    group_by(set) %>%
    summarise(sum = sum(num)) %>%
    arrange(desc(sum)) %>%
    mutate(., compounds_num = cumsum(sum)) %>%
    add_column(species_num = seq(1:nrow(.)))

curve_input_combine$compounds_percentage <- curve_input_combine$compounds_num/nrow(curve_input)

curve_input_combine$species_num <- as.factor(curve_input_combine$species_num)

b_update_annotation <- curve_input_combine

b_update <- ggplot(data=curve_input_combine, aes(x=species_num, y=compounds_percentage, group=1)) +
    geom_line(color = "#0072B2")+
    geom_point(color = "#0072B2")+
    geom_point(color = ifelse(curve_input_combine$compounds_percentage > 0.85, "#0072B2", "red"))+
    theme_minimal() +
    geom_hline(yintercept=0.8, linetype="dashed", color = "red") +
    geom_vline(xintercept=6, linetype="dashed", color = "red") +
    scale_y_continuous(

    # Features of the first axis
    name = "Ratio of covered compounds",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*43, name="N of covered compounds")
        ) +
    xlab("Minimum species rank") +
    theme(axis.text.x = element_text(angle = 45, size = 14, colour = "black", hjust = 1, face= "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 12, face = "bold"))

b_update_annotation$species <- str_split(b_update_annotation$set, pattern = "s__", simplify = TRUE)[,2] %>% str_replace_all("_", " ")

b_update_add_text <- b_update +
        geom_label_repel(data=b_update_annotation, label = b_update_annotation$species,
        fontface = 'bold.italic', color = 'black', box.padding = 0.5,
        #box.padding = unit(1.3, "lines"),
        segment.color = 'grey50', direction = "both", min.segment.length = unit(0, 'lines')) +
        #geom_text_repel(data=b_update_annotation, box.padding = 0.5, max.overlaps = Inf, label = b_update_annotation$species) + 
        xlim(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")) +
        theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)))

ggsave(file = paste(output_dir, "unique_compounds_secondary_curve_percent_add_text_only_comp.pdf", sep = "/"), b_update_add_text, width = 8)

ggsave(file = paste(output_dir, "unique_compounds_secondary_curve_percent_only_comp.pdf", sep = "/"), b_update, width = 12, height = 4)

save.image(paste0(output_dir, "/overlap_and_unique_compounds_secondary_only_comp.RData"))



