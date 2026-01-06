#STEP 1: Data preprocessing for ML#

BiocManager::install(c("caret","RANN","randomForest","Boruta","kernlab","pRDC"))
BiocManager::install("pROC")
# 1. Load Required Packages
# 'caret' is the main package used for the workflow.
# It often requires other packages for specific functions (e.g., 'RANN' for knnimpute).
library(caret)
library(RANN) # Might be needed for knnimpute
library(randomForest)
library(Boruta)
library(kernlab)
library(pROC)
library(ggplot2)
library(readr)
library(dplyr)
# 2. Load the Gene Expression Data
# The lecture does not provide the exact file path or loading function.
# Placeholder for loading your data (genes as rows, samples as columns):
# exp_data <- read.csv("your_gene_expression_matrix.csv", row.names = 1)
# Note: Ensure your data is numeric.
data <- read.table("data.csv",header = TRUE,sep = ",",row.names = 1)
phenotype_data<- phenotype_data%>%
                 mutate(characteristics_ch1=
                          recode(characteristics_ch1,
                                 "subject status: normal, healthy subject"="normal",
                                 "subject status: Chinese imported falciparum malaria patient"="malaria"
                                 ))

boxplot(data,
        outline = FALSE,
        col = "skyblue",
        main = "Before log Transformation")
par(mfrow = c(1,2))
str(data)

# 3. Log Transformation [06:46]
# Apply log10 transformation with a pseudo-count of +1 to normalize the right-skewed data.
hist(data[,12], main ="Before Log Transformation", col="skyblue")
hist(log10(data + 1)[,12], main ="After Log Transformation", col="coral")
data_log <- log10(data + 1)

# 4. Transpose Data [07:12]
# Transpose the matrix to meet ML requirements: Samples as rows, Genes (features) as columns.
data_t <- as.data.frame(t(data_log))

#  Pre-processing
# filtering near zero variables

nzv<- preProcess(data_log, method = "nzv", uniqueCut = 15)
data_t<- predict(nzv,data_t)

# Center and scale the data
processdata <- preProcess(data_t, method = c("center", "scale"))
data_t<- predict(processdata,data_t)

# Handling Missing Values if any
#knnimpute <- preProcess(data_t, method = "knnimpute")
#data_t<- predict(knnimpute,data_t)#

# Check for missing values after imputation (should be 0)
# sum(is.na(data_t))
anyNA(data_t)
sum(is.na(data_t))

# Add Class Labels
groups <- factor(phenotype_data$characteristics_ch1,
                 levels = c("malaria", "normal"),
                 label = c(1, 0))
class(groups)
levels(groups)
table(groups)

# 7. Bind Labels to Data 
# Create the final data frame with features and the response variable (Group).
data_t <- cbind(groups,data_t)

# 8. Data Partition (70/30 Split) 
# Split data into training and test sets using stratified sampling to preserve class proportions.

# Set seed for reproducible results
set.seed(123) # You can use any integer

# Create the stratified partition index (70% for training)
trainIndex <- createDataPartition(data_t$groups, # Split based on the target variable
                                  p = 0.70,             # Percentage of data for the training set
                                  list = FALSE          # Return indices as a simple vector
)

# Create the Training and Test sets [14:38]
train_data <- data_t[trainIndex, ]
test_data <- data_t[-trainIndex, ]

# 9. Separate Features (X) and Labels (Y) [15:13]
# Separate predictor variables (X) from the response variable (Y) to avoid data leakage.

# Training data
X_train <- train_data[, -1] # All columns except the last (Group)
Y_train <- train_data$groups                   # The last column (Group)

# Test data
X_test <- test_data[, -1]
Y_test <- test_data$groups

# X_train, Y_train, X_test, and Y_test are now ready for model training.


# 1. Load Required Packages
# 'Boruta' for the Boruta feature selection algorithm.
# 'caret' for the RFE (Recursive Feature Elimination) framework.
# 'randomForest' is implicitly used by Boruta and explicitly by the RFE function in caret.
library(Boruta)
library(caret)
library(randomForest) # Required for the rfFuncs and the model within RFE

#  Boruta Feature Selection (Random Forest Based) ---

# Set seed for reproducibility before any random operation
set.seed(42)

# 2. Run the Boruta algorithm 
# X_train: Predictor variables (genes)
# Y_train: Response variable (class labels: Normal/Malaria)
# doTrace=1 prints the progress of the algorithm to the console
boruta_results <- Boruta(x = X_train, 
                         y = Y_train, 
                         doTrace = 1)

# 3. Print the summary of results 
# Shows counts of Confirmed, Tentative, and Rejected features
print(boruta_results)

# 4. Fix Tentative features
# Finalizes the status (Confirmed or Rejected) of any features that were borderline (Tentative)
boruta_final <- TentativeRoughFix(boruta_results)
print(boruta_final) # View the final counts

# 5. Extract the confirmed important features (genes) [05:40]
Boruta_features <- getSelectedAttributes(boruta_final)
write.csv(Boruta_features,"boruta_selected_genes.csv")

# Check the number of selected features (e.g., 88 in the lecture's data)
# length(Boruta_features) 

# --- Section 2: Recursive Feature Elimination (RFE) ---

# 6. Define RFE Control Parameters [06:41]
# Specifies the training model, cross-validation method, and number of folds.
# functions = rfFuncs tells RFE to use Random Forest to rank features
# method = "cv" specifies cross-validation
# number = 3 sets up 3-fold cross-validation
rfe_contol <- rfeControl(functions = rfFuncs, 
                          method = "cv", 
                          number = 3,
                          verbose = TRUE)

# 7. Define the feature subset sizes to test
# The model will be trained and tested with 1, 2, 3, 4, 5, 10, and 50 features.
sizes <- c(1:5, 10, 50) 

# 8. Run RFE on the training data 
# X_train: Predictor variables (genes)
# Y_train: Response variable (class labels)
# sizes: The number of features RFE should test (optional, defaults to testing all)
rfe_results <- rfe(x = X_train, 
                   y = Y_train, 
                   sizes = sizes, 
                   rfeControl = rfe_contol)

# 9. Print the RFE results
# Shows the best subset size and the accuracy achieved
print(rfe_results)

# 10. Visualize RFE performance (Accuracy vs. Number of Features) 
 plot(rfe_results,type = c("g", "o"))

# 11. Extract the optimal features selected by RFE
rfe_selected_features <- predictors(rfe_results)
write.csv(rfe_selected_features,"rfe_selected_features.csv")


# Check the number of selected features 
# length(RFE_selected_features)

# --- Section 3: Check Overlap ---

# 12. Find the common features between Boruta and RFE
common_features <- intersect(Boruta_features, 
                             rfe_selected_features)

common_features2 <- intersect(boruta_final, 
                             rfe_selected_features)

# Check the number of common features
# length(common_features) # 
# 13. Save the results for the next session
# Save the selected features for use in model training 
# save(Boruta_selected_features, RFE_selected_features, common_features, 
#      file = "selected_features.RData")



#Model training and evaluation
# Subset data with Boruta features [01:08]

# X_train_boruta: Only Boruta selected genes in training set
train_boruta <- X_train[, Boruta_features]

# X_test_boruta: Only Boruta selected genes in test set
test_boruta <- X_test[, Boruta_features]

ctrl <- trainControl(method = "cv",
                     number = 10,
                     verboseIter = TRUE)
#Random forest training
RF_boruta <- train(x = train_boruta,
                   y = Y_train,
                   method = "rf",
                   importance = TRUE,
                   trControl = ctrl)
plot(varImp(RF_boruta), top = 10)

# Combine features and labels for the train function
training_data_boruta <- data.frame(X_train_boruta, Class = Y_train)


# Train Support Vector Machine (method = "svmRadial" for RBF kernel) [04:34]
# Note: SVM tuning grid is typically defined here, but using default for simplicity
svm_boruta <- train( x= train_boruta,
                     y= Y_train,
                     method = "svmRadial", 
                     trControl = ctrl,
                     tuneGrid = data.frame(C=c(0.25,0.5,1),
                                           sigma = 0.5),
                     prob.model= TRUE
)

# Train Artificial Neural Network (method = "nnet") [06:11]
# Note: A simple tuning grid for size (neurons) is defined here for demonstration
nnet_grid <- expand.grid(size = c(1, 2), # Number of hidden units
                         decay = c(0))   # Weight decay (regularization) [07:11]

ann_boruta <- train( x= train_boruta, 
                     y = Y_train, 
                     method = "nnet", 
                     trControl = ctrl,
                     tuneGrid = data.frame(size = 1:2,decay =0),
                     MaxNWts = 2000 # Adjust max iterations if necessary
)


# --- Evaluate Models ---

# Make Predictions (Class Labels) 
for (name in names(model_list_boruta)) 
  # Generate class predictions on the test set
  rf_pred_boruta <- predict(RF_boruta, newdata = test_boruta)
  svm_pred_boruta <- predict(svm_boruta, newdata = test_boruta)
  ann_pred_boruta <- predict(ann_boruta, newdata = test_boruta)
  
  # Calculate Confusion Matrix 
  rf_conf_boruta<- confusionMatrix(rf_pred_boruta,Y_test)
  svm_conf_boruta<- confusionMatrix(svm_pred_boruta,Y_test)
  ann_conf_boruta<- confusionMatrix(ann_pred_boruta,Y_test)
  

# Print Accuracy (example for RF) 
# confusion_matrices_boruta[["RF"]]$overall['Accuracy']

  results_boruta <- data.frame(
    model = c("Random Forest", "SVM","ANN"),
    Accuracy = c(confusionMatrix(rf_pred_boruta,Y_test)$overall["Accuracy"],
                 confusionMatrix(svm_pred_boruta,Y_test)$overall["Accuracy"],
                 confusionMatrix(ann_pred_boruta,Y_test)$overall["Accuracy"]
      
    )
  )
# Predict probabilities (type = "prob")
rf_prob_boruta <- predict(RF_boruta, newdata = test_boruta, type = "prob")
svm_prob_boruta <- predict(svm_boruta, newdata = test_boruta, type = "prob")
ann_prob_boruta <- predict(ann_boruta, newdata = test_boruta, type = "prob")

# Calculate and Plot ROC Curves/AUC (Probabilities) 
rf_roc_boruta <- roc(Y_test, as.numeric(rf_prob_boruta[,1]))
svm_roc_boruta <- roc(Y_test, as.numeric(svm_prob_boruta[,1]))
ann_roc_boruta <- roc(Y_test, as.numeric(ann_prob_boruta[,1]))

rf_auc_boruta <- auc(rf_roc_boruta)
svm_auc_boruta <- auc(svm_roc_boruta)
ann_auc_boruta <- auc(ann_roc_boruta)

dev.off()

# Combine results for ROC curve plotting (optional but common)
par(mfrow= c(1,3))
#RF
plot(rf_roc_boruta,
     col ="black",
     lwd = 2,
     lty = 1,
     main = paste("Random Forest\nAUC =", round(rf_auc_boruta,3)))
#SVM
plot(svm_roc_boruta,
     col ="green",
     lwd = 2,
     lty = 2,
     main = paste("SVM (Radial)\nAUC =", round(svm_auc_boruta,3)))
#ANN

plot(ann_roc_boruta,
     col ="blue",
     lwd = 2,
     lty = 3,
     main = paste("ANN\nAUC =", round(ann_auc_boruta,3)))

# Need the 'pROC' package for ROC/AUC
# library(pROC)
# roc_rf <- roc(Y_test, prob_rf$X1)
# auc_rf <- auc(roc_rf)
# plot(roc_rf, main="ROC Curves (Boruta Features)")
# ... (add other ROC curves)










