setwd("E:/backup/SmallMolecules/healthy_data/diseases_analysis/IBD/ML/IBD/")

library(tidyverse)
library(caTools)
library(caret)
library(Boruta)
library(varhandle)
library(cowplot)
library(ggpubr)
library(reshape2)


# load data ---------------------------------------------------------------

load("IBD_model_naturalCompounds.RData")
data <- read.delim("../../validation_cohort-ECs.tsv", check.names = F)
metadata <- read.delim("../../validation_cohort_metadata.tsv", check.names = F) %>% 
  select(c("sample_id", "study_condition"))

# process data ------------------------------------------------------------

metadata <- metadata %>% 
  rename("Group" = "study_condition") %>% 
  rename("Sample" = "sample_id")

data <- data %>% 
  rename("EC" = "# Gene Family") %>% 
  subset(EC %in% features$.)

# merge column
data <- remove_rownames(data) %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  merge(., metadata, by = "Sample") %>% 
  column_to_rownames("Sample") %>% 
  mutate(Group = gsub("control", "Control", Group))

test_data <- data[, 1:(ncol(data)-1)] %>% 
  unfactor() %>% 
  as.data.frame()

test_grouping <- as.factor(data[, ncol(data)])

test_data[] <- apply(test_data, 2, as.numeric) 


# test --------------------------------------------------------------------

set.seed(1)

confusion_matrix_ranger <- predict(ranger_mod, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(ranger_mod, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("Control", "IBD"),
                         predictor = probabilities_ranger[, "IBD"])


# plot ROC ----------------------------------------------------------------

test_grouping <- ifelse(test_grouping == "Control", 0, 1)

roc_colour <- PRROC::roc.curve(scores.class0 = probabilities_ranger$IBD, weights.class0 = test_grouping, curve = T)

ggplot(as.data.frame(roc_colour$curve), aes(x = V1, y = V2, colour = V3)) +
  geom_line(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey") +
  theme_test() +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(text = element_text(size = 15),
        legend.title = element_blank(),
        legend.key.height = unit(2, "cm"),
        aspect.ratio = 1) +
  ggtitle("IBD vs Control", subtitle = paste("AUC = ", sprintf("%.3f", pROC_ranger$auc), sep = "")) +
  scale_colour_gradientn(colours = rainbow(5),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         limits = c(0, 1))

# plot boxplot ------------------------------------------------------------

predictor <- pROC_ranger$predictor %>% 
  as.data.frame() %>% 
  add_column(Real = test_grouping)

ggplot(predictor, aes(x = Real, y = ., fill = Real)) +
  geom_boxplot(lwd = 1) +
  geom_abline(slope = 0, intercept = 0.5, linetype = "dashed", color = "grey") +
  theme_test() +
  ylab("Model prediction") +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("Control" = "#a9deff", 
                               "IBD" = "#558adf"))

# plot confusion matrix ---------------------------------------------------

confusion_matrix_ranger$table <- as.data.frame(confusion_matrix_ranger$table) %>% 
  add_column(., "Label" = ifelse(.$Prediction == .$Reference, "Hit", "Miss")) 

ggplot(data = confusion_matrix_ranger$table, mapping = aes(x = Prediction, y = Reference, fill = Label)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 0.5, fontface = "bold", alpha = 1, size = 7) +
  coord_fixed() +
  scale_fill_manual(values = c("#558adf", "white")) +
  scale_x_discrete(position = "top") +
  labs(x = "Prediction", y = "Reference") +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "none")