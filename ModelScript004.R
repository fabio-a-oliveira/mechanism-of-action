# Model004.R
# outcomes simplified via a few principal components
# one dedicated model to predict each principal component via random forest
# simplified random forest - predFixed = 80 (not optimal, but close enough), ntree = 50
# models configured to use as little memory as practicable

# Housekeeping ---------------------------------------------------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("data.table")) install.packages("data.table") else library("data.table")

if (!("parallel" %in% installed.packages())) {install.packages("parallel")}
if (!("doParallel" %in% installed.packages())) {install.packages("doParallel")}
if (!("foreach" %in% installed.packages())) {install.packages("foreach")}
if (!("beepr" %in% installed.packages())) {install.packages("beepr")}

options("datatable.print.topn" = 20)

# Declaration of some useful functions ---------------------------------------------------------------------------------

limitRange <- function(x){
  delta <- 10^(-15)
  pmin(pmax(x,delta),1-delta)}

score <- 
  function(reference,prediction){
    N <- nrow(reference)
    M <- ncol(reference)
    score <- - sum(reference*log(prediction) + (1-reference)*log(1-prediction)) / (N*M)
  }

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

proportionHeldOut <- .1

train$heldOut <- createDataPartition(train$pca.outcomes$x[,1], # considering 1st PC as output
                             p = proportionHeldOut, times = 1, list = FALSE)

# Create models for each of the principal components of the outcomes ---------------------------------------------------

numPrincipalComponents <- 10

grid <- expand.grid(predFixed = 80,
                    minNode = 2)

ctrl <- trainControl(method = "none",
                     verboseIter = TRUE, 
                     allowParallel = TRUE,
                     returnData = FALSE,
                     returnResamp = 'none',
                     savePredictions = FALSE,
                     classProbs = FALSE,
                     trim = TRUE)

# # create and register processor cluster
# clust <- parallel::makePSOCKcluster(8)
# doParallel::registerDoParallel(clust)

# train models
fit004_001 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,1],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_002 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,2],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_003 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,3],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_004 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,4],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_005 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,5],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_006 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,6],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_007 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,7],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_008 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,8],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_009 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,9],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

fit004_010 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,10],
        tuneGrid = grid,
        ntree = 50,
        trControl = ctrl)

# # stop cluster
# parallel::stopCluster(clust)
# doParallel::stopImplicitCluster()

# Measure performance --------------------------------------------------------------------------------------------------

# target, approximate and predicted outcomes
y <- 
  train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix()

y_app <- 
  train$pca.outcomes$x[train$heldOut,1:numPrincipalComponents] %*% 
  t(train$pca.outcomes$rotation[,1:numPrincipalComponents])

y_hat <- 
  cbind(predict(fit004_001, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_002, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_003, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_004, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_005, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_006, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_007, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_008, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_009, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
        predict(fit004_010, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1)) %>%  
  magrittr::multiply_by_matrix(t(train$pca.outcomes$rotation[,1:numPrincipalComponents]))

# calculate scores
score_target_app <- score(reference = y,
                          prediction = limitRange(y_app))
score_app_prediction <- score(reference = limitRange(y_app),
                              prediction = limitRange(y_hat))
score_target_prediction <- score(reference = y,
                                 prediction = limitRange(y_hat))
# print scores
score_target_app
score_app_prediction
score_target_prediction

# confusion matrix
# confusionMatrix(data = predictions %>% round(0) %>% as.factor(),
#                 reference = y %>% as.factor(),
#                 positive = "1")

# beep!
beepr::beep(5)

# Save models ------------------------------------------------------------------

save(list = paste("fit004",1:numPrincipalComponents %>% as.character() %>% str_pad(width=3,pad="0"), sep = "_"),
     file = "ModelFit004.Rdata")



