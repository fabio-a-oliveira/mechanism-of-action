# Model003.R
# features - gene and expression and cell viability raw data, categorical
# features remain untouched
# outcomes - first few principal components, full set reconstituted after modeling;
# choice of amount of PCs based on accuracy measurement
# individual models trained for each of the PCs
# 10% hold-out samples to evaluate results at the end
# no tuning - best outcome from model was with 240 random variables per tree, 
# but 120 and 80 give reasonable results in less time

# Housekeeping ---------------------------------------------------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("data.table")) install.packages("data.table") else library("data.table")

if (!("parallel" %in% installed.packages())) {install.packages("parallel")}
if (!("doParallel" %in% installed.packages())) {install.packages("doParallel")}
if (!("beepr" %in% installed.packages())) {install.packages("beepr")}

options("datatable.print.topn" = 20)

# Load data ------------------------------------------------------------------------------------------------------------

train <- 
  list(features = read.csv(file.path("files","train_features.csv")),
       targets = read.csv(file.path("files","train_targets_scored.csv")),
       nonscored = read.csv(file.path("files","train_targets_nonscored.csv")))

test <- 
  list(features = read.csv(file.path("files","test_features.csv")))

# Data wrangling  ------------------------------------------------------------------------------------------------------

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

train$pca.outcomes <-
  train$targets %>% 
  select(-sig_id) %>% 
  prcomp()

# Partition data - create hold-out samples -----------------------------------------------------------------------------

train$heldOut <- createDataPartition(train$pca.outcomes$x[,1], # considering 1st PC as output
                             p = .1, times = 1, list = FALSE)

# Select number of principal components to represent the outcomes ------------------------------------------------------

numPrincipalComponents <- 30

y <- 
  train$targets %>% 
  select(-sig_id) %>% 
  as.matrix()

y_app <-
  train$pca.outcomes$x[,1:numPrincipalComponents, drop = FALSE] %*%
  t(train$pca.outcomes$rotation[,1:numPrincipalComponents, drop = FALSE])
y_app <- pmax((pmin(y_app,1-10^(-15))),10^(-15))

# confusion matrix
confusionMatrix(reference = y %>% round(0) %>% as.factor(),
                data = y_app %>% round(0) %>% as.factor(),
                positive = "1")

# score
- sum(y * log(y_app) + (1-y) * log(1 - y_app)) / (206 * nrow(train$targets))

# RMSE
RMSE(y, y_app)

# Create model for 1st PCA of outcomes ---------------------------------------------------------------------------------

grid <- expand.grid(predFixed = c(40,60,80,120,160,240,480,600,872),
                    minNode = 2)

# create and register processor cluster
clust <- parallel::makePSOCKcluster(8)
doParallel::registerDoParallel(clust)

# train model
fit002 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,1],
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9, number = 1,
                                 verboseIter = TRUE, allowParallel = TRUE))

# stop cluster
parallel::stopCluster(clust)
doParallel::stopImplicitCluster()

# Measure performance and save model -----------------------------------------------------------------------------------

# predictions on held out samples
predictions <-
  predict(fit002,
          train$features[train$heldOut,] %>% select(-sig_id)) %>% 
  matrix(nrow = 1) %>% 
  t() %*% 
  t(train$pca.outcomes$rotation[,1])

RMSE(predictions,
     as.matrix(select(train$targets[train$heldOut,], - sig_id)))

p <- pmax(pmin(predictions, 1 − (10^(−15))), (10^(−15)))
y <- train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix()

score <- - sum(y * log(p) + (1-y) * log(1 - p)) / (206 * length(train$heldOut))

score

# confusion matrix
confusionMatrix(data = predictions %>% round(0) %>% as.factor(),
                reference = y %>% as.factor(),
                positive = "1")

# save output
save(list = "fit002", file = "ModelFit002.Rdata")

# beep!
beepr::beep(5)

