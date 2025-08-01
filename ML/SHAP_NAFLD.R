setwd("E:/backup/SmallMolecules/revision_NatMicrobiology/new_NASH_groups/")

library(doParallel)
library(caret)
library(MLeval)
library(varhandle)
library(tidyverse)
library(fastshap)
library(ggtext)

# load data
load("ML_noNASH_HealthySteatosis.RData")
model <- ranger_mod

norm_data <- model$trainingData %>%
  rownames_to_column("Sample") %>% 
  rename("Group" = ".outcome") %>% 
  mutate(Group = as.factor(Group))

raw_data <- training_data %>% 
  rownames_to_column("Sample") %>%  
  add_column(Group = training_grouping)

# train model (or load model)
set.seed(2)
model <-  train(x = norm_data[2:(ncol(norm_data)-1)],
                y = norm_data$Group,
                trControl = trainControl(method = "none", classProbs = TRUE),
                metric = "ROC",
                method = "ranger",
                importance = "impurity",
                num.trees = 500,
                max.depth = 3,
                tuneGrid = expand.grid(mtry =  1,
                                       splitrule = "gini",
                                       min.node.size = 1),
                num.threads = 1)

# pred_wrapper: To get the class probabilities for positive class
pfun <- function(object, newdata){
  predict(object, data = newdata)$predictions[, 1]
} 

# run fastshap (important to set seed!!)
set.seed(10)
final_model_shap <- explain(model$finalModel,
                            X = model$trainingData[, -ncol(model$trainingData)], 
                            pred_wrapper = pfun, 
                            nsim = 50) %>% 
  as.data.frame

# format tables for plots
gathered_final_model_shap <-
  tidyr::gather(final_model_shap, key = "key", value = "shap_value") %>%
  mutate(shap_value = shap_value * 100) %>% 
  group_by(key) %>%
  mutate(mean_shap_value = mean(shap_value)) %>%
  mutate(mean_abs_shap_value = mean(abs(shap_value))) %>%
  mutate(sum_shap_value = sum(shap_value)) %>%
  mutate(max_abs_shap_value = max(abs(shap_value))) %>%
  mutate(Group = raw_data$Group %>% factor(levels = c("Control", "NAFLD"))) %>%
  ungroup() %>% 
  mutate(normalised_abundance = norm_data[, 2:(ncol(norm_data)-1)] %>% tidyr::gather(key, value) %>% .$value) %>%
  mutate(abundance = raw_data[, 2:(ncol(raw_data)-1)] %>% tidyr::gather(key, value) %>% .$value) 

gathered_final_model_shap_contribution <-
  gathered_final_model_shap[, c("key", "shap_value", "abundance")] %>% 
  group_by(key) %>% 
  mutate(glm_slope = coef(glm(shap_value ~ abundance))[[2]]) %>% 
  mutate(mean_abs_shap_value = mean(abs(shap_value))) %>% 
  ungroup()

gathered_final_model_shap_contribution <-
  unique(gathered_final_model_shap_contribution[, c("key", "glm_slope", "mean_abs_shap_value")]) %>%
  mutate(sum_mean_abs_shap_value = sum(mean_abs_shap_value)) %>%
  mutate(contribution_pct = mean_abs_shap_value/sum_mean_abs_shap_value * 100)

# arrange features by mean shap value
top_features <- gathered_final_model_shap %>% 
  arrange(-mean_abs_shap_value) %>% 
  .$key %>% 
  unique() %>% 
  as.vector() 

# plots for feature contribution
final_model_shap_plot_left <-
  ggplot(data = gathered_final_model_shap_contribution %>% 
           mutate(Sign = ifelse(glm_slope > 0, "Control", "NAFLD")) %>% 
           mutate(key = gsub("g__\\w+.s__", "", key)) %>% 
           mutate(key = gsub("_", " ", key)) %>% 
           mutate(key = gsub("\\|", " \\| ", key))) +
  coord_flip() +
  geom_bar(aes(x = reorder(key, contribution_pct),
               y = contribution_pct,
               fill = Sign),
           stat = "identity",
           width = 0.65,
           alpha = 0.9) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 10, 2)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.2),
        legend.position= "top",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.title.x= element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c("Control" = "#ebdd86", "NAFLD" = "#b79f00")) +
  ylab("Feature contribution (%)")

final_model_shap_plot_right <- 
  ggplot(data = gathered_final_model_shap) +
  coord_flip(ylim = c(-4, 6)) +
  geom_hline(yintercept = 0, colour = "grey60") + # the y-axis beneath

  ggforce::geom_sina(aes(x = reorder(key, abs(mean_abs_shap_value)), 
                         y = shap_value, 
                         color = normalised_abundance),
                     method = "counts", 
                     maxwidth = 0.8, 
                     alpha = 0.5, 
                     size = 0.4) +
  scale_colour_gradientn(
    colors = c("#004688", "#ED0000"),
    breaks=c(0,.25),
    labels=c(" Low","High "),
    limits = c(0, 0.25),
    oob = scales::squish,
    guide = guide_colorbar(barwidth = 10, barheight = 0.5)) +
  
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size = 0.2),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="top",
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank()) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-6, 6, 2)) +
  labs(y = "SHAP value (impact on model output)", x = "", color = "")


final_model_shap_plot <- ggpubr::ggarrange(plotlist = list(final_model_shap_plot_left, final_model_shap_plot_right), 
                                           widths = c(3, 3),
                                           ncol = 2,
                                           align = "h")


# individual plots for each feature (x=feature abundance, y=shap value)
list_of_dependence_plot <- 
  lapply(top_features, function(unique_key){
    
    unique_key <- gsub("g__\\w+.s__", "", unique_key) %>% 
      gsub("\\|", " \\| ", .) %>% 
      gsub("_", " ", .)
    
    current_data <- gathered_final_model_shap %>% 
      mutate(key = gsub("g__\\w+.s__", "", key)) %>% 
      mutate(key = gsub("\\|", " \\| ", key)) %>% 
      mutate(key = gsub("_", " ", key)) %>% 
      filter(key == unique_key) 
  
    unique_abundance <- current_data$abundance %>% sort %>% unique
    
    best_cutoff <- lapply(unique_abundance, function(x){
      upper_left  <- current_data %>% filter(abundance < x , shap_value >= 0) %>% nrow
      upper_right <- current_data %>% filter(abundance >= x, shap_value >= 0) %>% nrow
      lower_left  <- current_data %>% filter(abundance < x , shap_value < 0) %>% nrow
      lower_right <- current_data %>% filter(abundance >= x, shap_value < 0) %>% nrow
      
      max_number <- max((upper_left + lower_right), (upper_right + lower_left))
      return(data.frame(x, max_number))
    }) %>% 
      Reduce(rbind, .) %>% 
      filter(max_number == max(max_number)) %>% .$x %>% .[1]
    
    current_plot <- 
      ggplot(data = current_data %>% 
               mutate(quan98 = quantile(abundance, na.rm = T, probs = 0.98)) %>% # remove extremely large outliers
               filter(abundance <= quan98),
             aes(x = abundance, y = shap_value, colour = Group)) +
      geom_point(shape = 20, size = 3,
                 alpha = 0.7) +
      geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", size = 0.2) + 
      geom_vline(xintercept = best_cutoff, color = "grey50", linetype = "dashed", size = 0.2) + 
      ylab("SHAP") +
      xlab(unique_key %>% 
             gsub(" \\| ", "\n", .)) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(size = 0.2), 
            axis.text = element_text(size = 8),
            aspect.ratio = 1) +
      scale_color_manual(values=c("Control" = "#ebdd86", "NAFLD" = "#b79f00"))
    
    return(current_plot)
  })

ggpubr::ggarrange(plotlist = list_of_dependence_plot, ncol = 5, nrow = 5, common.legend = T, align = "hv")