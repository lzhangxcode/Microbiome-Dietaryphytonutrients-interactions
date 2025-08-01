setwd("E:/backup/SmallMolecules/healthy_data/diseases_analysis/CRC/")

library(tidyverse)
library(caTools)
library(caret)
library(Boruta)
library(varhandle)
library(cowplot)
library(ggpubr)
library(reshape2)


# load data ---------------------------------------------------------------

data <- read.delim("data/merged_ECs_CRC_stratified.tsv", check.names = F)
metadata <- read.delim("metadata/merged_metadata_CRC.tsv")
nutrichem <- read.delim("data/CRC_foods_compounds.tsv")$ECs %>% 
  strsplit(", ") %>% 
  unlist() %>% 
  unique()


# process data ------------------------------------------------------------

test_samples <- subset(metadata, grepl("YuJ", study_name))

metadata <- metadata %>% 
  select(c("sample_id", "study_condition", "age", "gender", "BMI")) %>% 
  rename("Group" = "study_condition") %>% 
  rename("Sample" = "sample_id")

data <- data %>% 
  subset(!grepl("UNGROUPED", `# Gene Family`)) %>% 
  subset(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  subset(!grepl("\\|unclassified", `# Gene Family`)) %>% 
  subset(grepl("s__", `# Gene Family`)) %>% 
  rename("EC" = "# Gene Family") %>% 
  setNames(gsub("_Abundance-RPKs", "", names(.)))

# split EC column into EC and Species
data <- cbind(colsplit(data$EC, "\\|", c("EC", "Specie")), data[, 2:ncol(data)]) %>% 
  subset(EC %in% nutrichem) 

# merge column
data <- unite(data, EC, c(EC, Specie), sep = "|") %>% 
  remove_rownames() %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  merge(., metadata, by = "Sample") %>% 
  .[, -(10605:10607)] %>% # remove metadata columns
  column_to_rownames("Sample") %>% 
  subset(Group != "adenoma")


# split -------------------------------------------------------------------

training_set <- subset(data, !(rownames(data)) %in% test_samples$sample_id)

training_data <- training_set[, 1:(ncol(training_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() 

training_grouping <- as.factor(training_set[, ncol(training_set)])

test_set <- subset(data, rownames(data) %in% test_samples$sample_id)

test_data <- test_set[, 1:(ncol(test_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame()

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
  
  # get top 20 features from this iteration
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

confusion_matrix_ranger <- predict(ranger_mod, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(ranger_mod, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("control", "CRC"),
                         predictor = probabilities_ranger[, "CRC"])


# plot ROC ----------------------------------------------------------------

test_grouping <- ifelse(test_grouping == "control", 0, 1)

roc_colour <- PRROC::roc.curve(scores.class0 = probabilities_ranger$CRC, weights.class0 = test_grouping, curve = T)

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
  ggtitle("CRC vs Control", subtitle = paste("AUC = ", sprintf("%.3f", pROC_ranger$auc), sep = "")) +
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
  scale_fill_manual(values = c("control" = "#6182a8", 
                               "CRC" = "#d86b6f"))


# plot confusion matrix ---------------------------------------------------

confusion_matrix_ranger$table <- as.data.frame(confusion_matrix_ranger$table) %>% 
  add_column(., "Label" = ifelse(.$Prediction == .$Reference, "Hit", "Miss")) 

ggplot(data = confusion_matrix_ranger$table, mapping = aes(x = Prediction, y = Reference, fill = Label)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 0.5, fontface = "bold", alpha = 1, size = 7) +
  coord_fixed() +
  scale_fill_manual(values = c("#3b9dda", "white")) +
  scale_x_discrete(position = "top") +
  labs(x = "Prediction", y = "Reference") +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "none")