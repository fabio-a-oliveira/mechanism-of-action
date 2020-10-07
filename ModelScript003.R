# Model003.R
# same as model 002, but trying to increase speed and reduce model size by
# reducing number of trees

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

# Create model for 1st PCA of outcomes with 50 trees ---------------------------

grid <- expand.grid(predFixed = c(40,60,80,120,160,240),
                    minNode = 2)

# create and register processor cluster
clust <- parallel::makePSOCKcluster(8)
doParallel::registerDoParallel(clust)

# train model with 50 trees
fit003_50 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,1],
        tuneGrid = grid,
        ntree = 50,
        trControl = trainControl(method = "LGOCV", p = .9, number = 1,
                                 verboseIter = TRUE, allowParallel = TRUE))

# stop cluster
parallel::stopCluster(clust)
doParallel::stopImplicitCluster()

# Measure performance ----------------------------------------------------------

# predictions on held out samples
predictions <-
  predict(fit003_50,
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

# beep!
beepr::beep(5)

# Create model for 1st PCA of outcomes with 100 trees --------------------------

grid <- expand.grid(predFixed = c(40,60,80,120,160,240),
                    minNode = 2)

# create and register processor cluster
clust <- parallel::makePSOCKcluster(8)
doParallel::registerDoParallel(clust)

# train model with 100 trees
fit003_100 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,1],
        tuneGrid = grid,
        ntree = 100,
        trControl = trainControl(method = "LGOCV", p = .9, number = 1,
                                 verboseIter = TRUE, allowParallel = TRUE))

# stop cluster
parallel::stopCluster(clust)
doParallel::stopImplicitCluster()

# Measure performance ----------------------------------------------------------

# predictions on held out samples
predictions <-
  predict(fit003_100,
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

# beep!
beepr::beep(5)

# Create model for 1st PCA of outcomes with 200 trees --------------------------

grid <- expand.grid(predFixed = c(40,60,80,120,160,240),
                    minNode = 2)

# create and register processor cluster
clust <- parallel::makePSOCKcluster(8)
doParallel::registerDoParallel(clust)

# train model with 50 trees
fit003_200 <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$pca.outcomes$x[-train$heldOut,1],
        tuneGrid = grid,
        ntree = 200,
        trControl = trainControl(method = "LGOCV", p = .9, number = 1,
                                 verboseIter = TRUE, allowParallel = TRUE))

# stop cluster
parallel::stopCluster(clust)
doParallel::stopImplicitCluster()

# Measure performance ----------------------------------------------------------

# predictions on held out samples
predictions <-
  predict(fit003_200,
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

# beep!
beepr::beep(5)

# Save models ------------------------------------------------------------------

save(list = c("fit003_50", "fit003_100", "fit003_200"), 
     file = "ModelFit003.Rdata")



