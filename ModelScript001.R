# Model001.R
# features - gene and expression and cell viability via PCA, categorical
# features remain untouched
# outcomes - reduced set of outcomes via first few PCs

# Housekeeping -----------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("data.table")) install.packages("data.table") else library("data.table")

if (!("parallel" %in% installed.packages())) {install.packages("parallel")}
if (!("doParallel" %in% installed.packages())) {install.packages("doParallel")}
if (!("beepr" %in% installed.packages())) {install.packages("beepr")}

options("datatable.print.topn" = 20)

# Load data --------------------------------------------------------------------

train <- 
  list(features = read.csv(file.path("files","train_features.csv")),
       targets = read.csv(file.path("files","train_targets_scored.csv")),
       nonscored = read.csv(file.path("files","train_targets_nonscored.csv")))

test <- 
  list(features = read.csv(file.path("files","test_features.csv")))

# Data wrangling  --------------------------------------------------------------

train$features <-
  train$features %>% 
  mutate(sig_id = as.factor(sig_id),
         cp_type = as.factor(cp_type),
         cp_time = as.factor(cp_time),
         cp_dose = as.factor(cp_dose))

train$targets <-
  train$targets %>% 
  mutate(sig_id = as.factor(sig_id))

test$features <-
  test$features %>% 
  mutate(sig_id = as.factor(sig_id),
         cp_type = as.factor(cp_type),
         cp_time = as.factor(cp_time),
         cp_dose = as.factor(cp_dose))

train$pca.features <- 
  train$features %>% 
  select(- c("sig_id", "cp_type", "cp_time", "cp_dose")) %>% 
  as.matrix() %>% 
  prcomp

train$pca.outcomes <-
  train$targets %>% 
  select(-sig_id) %>% 
  prcomp()

# Create model with PCAs applied to features and outcomes ----------------------

grid <- expand.grid(predFixed = c(40,60,80,120,160,240,480,600,872),
                    minNode = 2)

numPrincipalComponents <- 872

# create and register processor cluster
clust <- parallel::makePSOCKcluster(8)
doParallel::registerDoParallel(clust)

# train model
fit_Rborist_872 <-
  train(method = "Rborist",
        x = cbind(select(train$features, cp_type, cp_dose, cp_time),
                  train$pca.features$x[,1:numPrincipalComponents]),
        y = train$pca.outcomes$x[,1],
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9, number = 1,
                                 verboseIter = TRUE, allowParallel = TRUE))

# stop cluster
parallel::stopCluster(clust)
doParallel::stopImplicitCluster()

predictions <-                                                                                                          # OVERTRAINING! PREDICTIONS CALCULATED ON TRAINING SET!
  predict(fit_Rborist_872,
          cbind(select(train$features, cp_type, cp_dose, cp_time),
                train$pca.features$x[,1:numPrincipalComponents])) %>% 
  matrix(nrow = 1) %>% 
  t() %*% 
  t(train$pca.outcomes$rotation[,1])

RMSE(predictions,
     as.matrix(select(train$targets, - sig_id)))

p <- pmax(pmin(predictions, 1 − (10^(−15))), (10^(−15)))
y <- train$targets %>% select(-sig_id) %>% as.matrix()

score <- - sum(y * log(p) + (1-y) * log(1 - p)) / (206 * 23814)

score

# confusion matrix
confusionMatrix(data = predictions %>% round(0) %>% as.factor(),
                reference = y %>% as.factor(),
                positive = "1")

# beep!
beepr::beep(5)
