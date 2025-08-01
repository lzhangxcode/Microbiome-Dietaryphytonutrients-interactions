setwd("E:/backup/SmallMolecules/healthy_data/diseases_analysis/NAFLD/")

library(tidyverse)
library(caTools)
library(caret)
library(Boruta)
library(varhandle)
library(cowplot)
library(ggpubr)
library(reshape2)


# load data ---------------------------------------------------------------

foods <- read.delim("../../nutrichem/plant_disease_associations.tsv") %>% 
  subset(DiseaseID == "DOID:9452") %>% 
  .$PlantID %>% 
  unique()

database <- read.delim("../../nutrichem/full_database_filtered.tsv") %>% 
  subset(PlantID %in% foods) %>% 
  .$ECs %>% 
  strsplit(., ", ") %>% 
  unlist() %>% 
  unique()

NASH_samples <- read.csv("../../../revision_NatMicrobiology/NASH_cohort_groups.csv", sep = ";")

data_nafld_strat <- read.delim("merged_NAFLD_ECs.tsv", check.names = F) %>% 
  rename_all(~gsub("_Abundance-RPKs", "", .x)) 
  
data_control_strat <- read.delim("merged_Control_ECs.tsv", check.names = F) %>% 
  rename_all(~gsub("_Abundance-RPKs", "", .x)) 


# process data ------------------------------------------------------------

data_strat <- merge(data_nafld_strat, data_control_strat, by = "# Gene Family") %>% 
  rename("EC" = "# Gene Family") %>% 
  subset(grepl("s__", EC))

data_strat <- cbind(colsplit(data_strat$EC, "\\|", c("EC", "Specie")), data_strat[, 2:ncol(data_strat)]) %>% 
  subset(EC %in% database) %>% 
  unite(., EC, c(EC, Specie), sep = "|") %>% 
  remove_rownames() %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>%
  add_column(Group = c(rep("NAFLD", ncol(data_nafld_strat) -1), rep("Control", ncol(data_control_strat) -1))) 

# split -------------------------------------------------------------------

seeds <- sample.int(10000, 100)

training_set <- subset(data_strat, !(Sample %in% NASH_samples$Sample))

training_data <- training_set[, 1:(ncol(training_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("Sample")

training_grouping <- as.factor(training_set[, ncol(training_set)])

# run boruta 100 times
features <- data.frame(Iteration = 1:100,
                       Features = NA)

pb <- txtProgressBar(min = 0, max = 100, style = 3)

for (i in 1:100) {
  
  set.seed(seeds[i])
  boruta_res <- Boruta(x = training_data, 
                       y = training_grouping, 
                       doTrace = 0)  
  
  features_i <- getSelectedAttributes(boruta_res)
  
  stats <- attStats(boruta_res) %>% 
    subset(., decision == "Confirmed") %>% 
    .[order(.$meanImp, decreasing = T), ] %>% 
    .[1:20, ]
  
  features_i <- subset(features_i, features_i %in% rownames(stats))
  
  features$Features[i] <- paste(features_i, collapse = ", ")
  
  setTxtProgressBar(pb, i)
  
}

close(pb)

# get feature list
features <- features$Features %>% 
  strsplit(", ") %>% 
  unlist() %>% 
  table() %>% 
  as.data.frame()


# train -------------------------------------------------------------------

training_data <- select_if(training_data, names(training_data) %in% features$.)

trControl <- trainControl(method = "cv",
                          number = 5,
                          search = "random",
                          sampling = "down",
                          verboseIter = F,
                          returnData = T,
                          savePredictions = T,
                          allowParallel = T,
                          classProbs = T,
                          seeds = NULL)

set.seed(1)

ranger_mod <- train(x = training_data,
                    y = training_grouping,
                    trControl = trControl,
                    method = "rf")

# test --------------------------------------------------------------------

NASH_samples <- subset(NASH_samples, Group %in% c("Healthy", "Steatosis"))

test_set <- subset(data_strat, Sample %in% NASH_samples$Sample)

test_data <- test_set[, 1:(ncol(test_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("Sample")

test_grouping <- as.factor(test_set[, ncol(test_set)])

test_data <- select_if(test_data, names(test_data) %in% features$.) 

set.seed(1)

confusion_matrix_ranger <- predict(ranger_mod, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(ranger_mod, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("Control", "NAFLD"),
                         predictor = probabilities_ranger[, "NAFLD"])


# plot ROC ----------------------------------------------------------------

predictor <- pROC_ranger$predictor %>% 
  as.data.frame() %>% 
  add_column(Real = test_grouping)

test_grouping <- ifelse(test_grouping == "Control", 0, 1)

roc_colour <- PRROC::roc.curve(scores.class0 = probabilities_ranger$NAFLD, weights.class0 = test_grouping, curve = T)

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
  ggtitle("NAFLD vs Control", subtitle = paste("AUC = ", sprintf("%.3f", pROC_ranger$auc), sep = "")) +
  scale_colour_gradientn(colours = rainbow(5),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         limits = c(0, 1))

# plot boxplot ------------------------------------------------------------

ggplot(predictor, aes(x = Real, y = ., fill = Real)) +
  geom_boxplot(lwd = 1) +
  geom_abline(slope = 0, intercept = 0.5, linetype = "dashed", color = "grey") +
  theme_test() +
  ylab("Model prediction") +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12)) +
  scale_fill_manual(values = c("Control" = "#b8f685", 
                               "NAFLD" = "#107c00"))

# plot confusion matrix ---------------------------------------------------

confusion_matrix_ranger$table <- as.data.frame(confusion_matrix_ranger$table) %>% 
  add_column(., "Label" = ifelse(.$Prediction == .$Reference, "Hit", "Miss")) 

ggplot(data = confusion_matrix_ranger$table, mapping = aes(x = Prediction, y = Reference, fill = Label)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 0.5, fontface = "bold", alpha = 1, size = 7) +
  coord_fixed() +
  scale_fill_manual(values = c("#107c00", "white")) +
  scale_x_discrete(position = "top") +
  labs(x = "Prediction", y = "Reference") +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "none")