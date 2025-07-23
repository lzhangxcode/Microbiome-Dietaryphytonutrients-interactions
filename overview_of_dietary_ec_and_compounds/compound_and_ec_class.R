# script for compound class & EC classes 
# author:lu.zhang
# 2024.04.25 

library(tidyverse)
################################################################

# prepare input 

################################################################

# input needed for the full EC classes file 
np_database_overall_fullEC <- read_tsv("nutrichem_md_manual_with_non_microbiotaEC_fullEC.tsv")

# input needed for the microbiota EC ~ compounds file 
microbiota_np_ec <- readRDS("modified_ec_profile_after_np.rds")

# nutrichem original database 
nutrichem <- read_tsv("compound_associations.tsv")

# full ec - compound database 
full_db_ec <- read_tsv("database_with_full_enzyme_name_compound.tsv")

################################################################

# check compounds number + ec number

################################################################

# ec contained 
full_db_ec_np <- full_db_ec %>% 
  filter(Compounds %in% np_database_overall_fullEC$CompoundID) 

length(unique(microbiota_np_ec$ec_simplied)) #1908 


################################################################

# plot the barplots for the compounds 

################################################################

nutrichem_md <- nutrichem %>% 
  mutate(class3 = str_split(CompoundClass, pattern = "__", simplify = T) %>% .[,2]) %>% 
  mutate(class3 = ifelse(class3 == "", "Unclassified", class3)) %>% 
  mutate(class1 = str_split(class3, pattern = ";", simplify = T) %>% .[,1]) %>% 
  mutate(class2 = str_split(class3, pattern = ";", simplify = T) %>% .[,2]) %>% 
  mutate(class2 = paste(class1, class2, sep = ";"))

# 1388 compounds  
nutrichem_md_compoundswithec <- nutrichem_md %>% 
  filter(CompoundID %in% unique(np_database_overall_fullEC$CompoundID))

# 755 compounds
nutrichem_md_compoundsmicrobiota <- nutrichem_md %>% 
  filter(CompoundID %in% unique(microbiota_np_ec$Compounds))

extract_class_number <- function(input_df = NA){
  # extract the class1 
  class_df <- input_df %>% select(CompoundID, class1, class2, class3, CompoundClass) %>% unique()
  class1_df <- table(class_df$class1) %>% as.data.frame()
  # barplot 
  df_barplot <- class1_df
  df_barplot <- df_barplot %>% 
    filter(Var1 != "Others")
  df_barplot$Freq[df_barplot$Var1 == "Unclassified"] <- 0
  df_barplot$Var1[df_barplot$Var1 == "Unclassified"] <- "Others"
 
  df_barplot <- df_barplot %>% 
    arrange(Freq) %>% 
    mutate(Var1 = factor(Var1, levels = Var1))
  return(df_barplot)
  
}

nutrichem_md_compoundswithec_plotinput <- extract_class_number(input_df = nutrichem_md_compoundswithec)
nutrichem_md_compoundsmicrobiota_plotinput <- extract_class_number(input_df = nutrichem_md_compoundsmicrobiota)

# draw the plot 

draw_barplot <- function(plot_input = NA){
  
  df_barplot_p <- ggplot(plot_input, aes(x = Freq, y = Var1, fill = Var1))+
    geom_col(width = 0.7) +
    xlab("N of Compounds") + 
    ylab("Compound classes") + 
    ggsci::scale_fill_jco() + 
    theme_bw() + 
    theme(axis.text.y = element_text(size = 7), 
          axis.text.x = element_text(size = 4),
          axis.title=element_text(size=8, face = "bold"),
          legend.position="none") + 
    annotate("text", x = 1.5, y = 1.25, label = "....")+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  
  return(df_barplot_p)
  
}

output_dir <- "./compoundclass_ECclass/"

# combine these two category together 

draw_barplot_combine <- function(plot_input = NA){
  
  df_barplot_p <- ggplot(plot_input, aes(x = Freq, y = Var1, fill = Group))+ #alpha = Group
    #geom_col(width = 0.7) +
    geom_bar(position="dodge", stat="identity", width = 0.9) + 
    xlab("N of phytonutrients") + 
    ylab("Phytonutrient classes") + 
    ggsci::scale_fill_jco() + 
    theme_bw() + 
    theme(axis.text.y = element_text(size = 15), 
          axis.text.x = element_text(size = 12),
          axis.title=element_text(size=15, face = "bold"),
          legend.position = c(.75, .15),
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 12),
          legend.background = element_rect(fill="white",
                                           size=0.2, linetype="solid", 
                                           colour ="black")) + 
    #scale_alpha_manual(values = c(1, 0.55)) + 
    scale_fill_manual(values = c("#708aaa", "#00bfac")) + 
    annotate("text", x = 1.5, y = 1.05, label = "....") + 
    guides(shape = guide_legend(override.aes = list(size = 0.02)),
           fill = guide_legend(override.aes = list(size = 0.02))) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  
  return(df_barplot_p)
  
}


combine_plot_input <- bind_rows(nutrichem_md_compoundswithec_plotinput %>% mutate(Group = "Phytonutrients associated with EC"),
                                nutrichem_md_compoundsmicrobiota_plotinput %>% mutate(Group = "Phytonutrients biotransf. by microbiota"))

combine_plot_input$Group_old <- combine_plot_input$Group
combine_plot_input$Group[combine_plot_input$Group == "Phytonutrients associated with EC"] <- "Total"
combine_plot_input$Group[combine_plot_input$Group == "Phytonutrients biotransf. by microbiota"] <- "Microbiome-associated"
combine_plot_input$Group <- factor(combine_plot_input$Group, levels = c("Microbiome-associated","Total"))

combine_plot_compounds <- draw_barplot_combine(plot_input = combine_plot_input)

ggsave(str_c(output_dir, "nutrichem_class_compounds_combine_plot_compounds_barplot.pdf"), combine_plot_compounds, dpi = 600, width = 8)
ggsave(str_c(output_dir, "nutrichem_class_compounds_combine_plot_compounds_barplot.png"), combine_plot_compounds, dpi = 600, width = 8)

################################################################

# add the ec class 

################################################################
#1. -. -.-  Oxidoreductases.
#2. -. -.-  Transferases.
#3. -. -.-  Hydrolases.
#4. -. -.-  Lyases.
#5. -. -.-  Isomerases.
#6. -. -.-  Ligases.
#7. -. -.-  Translocases.
# https://enzyme.expasy.org/cgi-bin/enzyme/enzyme-search-cl?1

ec_classes <- tibble(ec_classes = c("1", "2", "3", "4", "5", "6", "7"), 
                     ec_classes_name = c("Oxidoreductases",
                                         "Transferases",
                                         "Hydrolases",
                                         "Lyases",
                                         "Isomerases",
                                         "Ligases",
                                         "Translocases"))

microbiota_ec_list <- tibble(ec = microbiota_np_ec$ec_simplied %>% unique()) %>% 
  mutate(ec_classes = ec %>% str_split(., pattern = "\\.", simplify = T) %>% .[,1]) 
full_ec_list <- tibble(ec = full_db_ec_np$ec_simplied %>% unique()) %>% 
  mutate(ec_classes = ec %>% str_split(., pattern = "\\.", simplify = T) %>% .[,1]) %>% 
  left_join(., ec_classes) %>% 
  filter(!ec == "0.0.0.0") %>% mutate(Group = "All ECs")


microbiota_ec_list_freq <- table(microbiota_ec_list$ec_classes_name) %>% 
  as.data.frame() %>% 
  mutate(Group = "Microbiota ECs")
full_ec_list_freq <- table(full_ec_list$ec_classes_name) %>% 
  as.data.frame() %>% 
  mutate(Group = "All ECs") %>% 
  arrange(Freq) 

join_ec_clases <- bind_rows(microbiota_ec_list_freq, full_ec_list_freq) %>% 
  arrange(Freq) %>% 
  mutate(Var1 = factor(Var1, levels = full_ec_list_freq$Var1))
  


# combine these two category together 

draw_barplot_combine_ec <- function(plot_input = NA){
  
  df_barplot_p <- ggplot(plot_input, aes(x = Freq, y = Var1, fill = Group))+ #alpha = Group
    geom_bar(position="dodge", stat="identity", width = 0.8) + 
    xlab("N of ECs") + 
    ylab("EC classes") + 
    ggsci::scale_fill_jco() + 
    theme_bw() + 
    theme(axis.text.y = element_text(size = 15), 
          axis.text.x = element_text(size = 12),
          axis.title=element_text(size=15, face = "bold"),
          legend.position = c(.75, .15),
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 12),
          legend.background = element_rect(fill="white",
                                           size=0.2, linetype="solid", 
                                           colour ="black")) + 
    scale_fill_manual(values = c("#708aaa", "#00bfac"))+
    guides(shape = guide_legend(override.aes = list(size = 0.01)),
           fill = guide_legend(override.aes = list(size = 0.01))) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  
  return(df_barplot_p)
  
}


join_ec_clases$Group_old <- join_ec_clases$Group
join_ec_clases$Group[join_ec_clases$Group == "All ECs"] <- "Total"
join_ec_clases$Group[join_ec_clases$Group == "Microbiota ECs"] <- "Microbiome-associated"
join_ec_clases$Group <- factor(join_ec_clases$Group, levels = c("Microbiome-associated","Total"))

combine_plot_ecs <- draw_barplot_combine_ec(plot_input = join_ec_clases)

ggsave(str_c(output_dir, "ec_class_combine_barplot.pdf"), combine_plot_ecs, dpi = 600, width = 8)


combine_plot_ecs_compounds <- ggpubr::ggarrange(combine_plot_compounds, combine_plot_ecs, widths = c(1, 0.9))
ggsave(str_c(output_dir, "ec_and_compound_class_combine_barplot.pdf"), combine_plot_ecs_compounds, dpi = 600, width = 16)

################################################################

# saving RData 

################################################################

save.image(str_c(output_dir, "extract_compound_ec_class.RData"))
load(str_c(output_dir, "extract_compound_ec_class.RData"))

combine_plot_ecs_compounds <- ggpubr::ggarrange(combine_plot_compounds, combine_plot_ecs+ xlab("N of enzymes"), widths = c(1, 0.9))
ggsave(str_c(output_dir, "ec_and_compound_class_combine_barplotv2.pdf"), combine_plot_ecs_compounds, dpi = 600, width = 16)

# change again the figure 
combine_plot_compoundsv3 <- combine_plot_compounds + 
    theme(axis.line = element_line(colour = "black"),
	  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 

combine_plot_ecsv3 <- combine_plot_ecs + xlab("N of enzymes") + 
    theme(axis.line = element_line(colour = "black"),
	  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())


combine_plot_ecs_compoundsv3 <- ggpubr::ggarrange(combine_plot_compoundsv3, combine_plot_ecsv3, widths = c(1, 0.9))

ggsave(str_c(output_dir, "ec_and_compound_class_combine_barplotv3.pdf"), combine_plot_ecs_compoundsv3, dpi = 600, width = 16)
