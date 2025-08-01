setwd("E:/backup/SmallMolecules/healthy_data/diseases_analysis/IBD/")

library(tidyverse)
library(readxl)
library(caTools)
library(Boruta)
library(caret)
library(varhandle)
library(cowplot)
library(ggpubr)
library(reshape2)


# load and format data ----------------------------------------------------

disease_plants <- read_excel("../../../IBD_data/IBD_UC_foods.xlsx")$PlantID
nutrichem <- read_tsv("../../../healthy_data/nutrichem/full_database_filtered.tsv") %>% 
  subset(., PlantID %in% disease_plants) %>% 
  .$ECs %>% 
  strsplit(", ") %>% 
  unlist() %>% 
  unique()

data <- read.csv("../../../IBD_data/vs_controls/data/clean_humann3_full.csv") %>% 
  subset(grepl("s__", ECs)) %>% 
  cbind(colsplit(.$ECs, "\\|", c("EC", "Specie")), .[, 2:ncol(.)]) %>% 
  subset(EC %in% nutrichem) %>% 
  unite(., EC, c(EC, Specie), sep = "|") %>% 
  remove_rownames() %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample")

metadata <- read.csv("../../../IBD_data/vs_controls/metadata.csv") %>% 
  mutate(diagnosis = ifelse(diagnosis == "nonIBD", "Control", "IBD")) %>% 
  subset(sample_alias %in% data$Sample) %>% 
  unique() %>% 
  select(c("sample_alias", "diagnosis")) %>% 
  rename("Sample" = "sample_alias")

data <- merge(data, metadata, by = "Sample")


# split -------------------------------------------------------------------

set.seed(1)
split <- sample.split(data$diagnosis, SplitRatio = .8)

training_set <- data[split == T, ]

training_data <- training_set[, 1:(ncol(training_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("Sample")

training_grouping <- as.factor(training_set[, ncol(training_set)])

test_set <- data[split == F, ]

test_data <- test_set[, 1:(ncol(test_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("Sample")

test_grouping <- as.factor(test_set[, ncol(test_set)])


# boruta 100 times --------------------------------------------------------

features <- data.frame(Iteration = 1:100,
                       Features = NA)

seeds <- sample.int(10000, 100)

# run boruta 100 times
for (i in 1:100) {
  
  set.seed(seeds[i])
  boruta_res <- Boruta(x = training_data, 
                       y = training_grouping, 
                       doTrace = 0)  
  
  features_i <- getSelectedAttributes(boruta_res)
  
  # select top 20 features from each iteration
  stats <- attStats(boruta_res) %>% 
    subset(., decision == "Confirmed") %>% 
    .[order(.$meanImp, decreasing = T), ] %>% 
    .[1:20, ]
  
  features_i <- subset(features_i, features_i %in% rownames(stats))
  
  features$Features[i] <- paste(features_i, collapse = ", ")
  
  setTxtProgressBar(pb, i)
  
}

# get final feature list
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

test_data <- select_if(test_data, names(test_data) %in% features$.)

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