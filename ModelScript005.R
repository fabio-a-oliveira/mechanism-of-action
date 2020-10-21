# Model005.R

# quick experiments around configuration from model 4

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
    score <- - sum(reference*log(limitRange(prediction)) + (1-reference)*log(1-limitRange(prediction))) / (N*M)
    score
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

# Simplify problem by restricting to single most common combination of categorical values ------------------------------

train$features <-
  train$features %>% 
  # filter(cp_time == "48", cp_type == "trt_cp", cp_dose == "D1") %>% 
  select(- c(cp_time, cp_type, cp_dose))

train$targets <-
  train$targets %>% 
  filter(sig_id %in% train$features$sig_id)

# Run Principal Component Analysis -------------------------------------------------------------------------------------

train$pca.outcomes <-
  train$targets %>% 
  select(-sig_id) %>% 
  prcomp()

train$pca.features <-
  train$features %>% 
  select(-sig_id) %>% 
  prcomp()

# Partition data - create hold-out samples -----------------------------------------------------------------------------

# proportionHeldOut <- .1 # for actual training
# proportionHeldOut <- .98 # for quick experiments
proportionHeldOut <- .5 # for quick lm (can't have deficient rank..)

# train$heldOut <- createDataPartition(train$pca.outcomes$x[,1], # considering 1st PC as output
#                              p = proportionHeldOut, times = 1, list = FALSE)
train$heldOut <- sample(x = 1:nrow(train$targets),
                        size = round(proportionHeldOut * nrow(train$targets),0),
                        replace = FALSE)


# Create models for each of the principal components of the outcomes ---------------------------------------------------

numPrincipalComponents <- 10

# to use Rborist
algorithm <- "Rborist"
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


# for evaluating influence of minNode

# grid <- expand.grid(predFixed = 80,
#                     minNode = c(2,4,8))
# 
# ctrl <- trainControl(method = "LGOCV", p = .9, number = 1,
#                      verboseIter = TRUE, 
#                      allowParallel = TRUE,
#                      returnData = FALSE,
#                      returnResamp = 'none',
#                      savePredictions = FALSE,
#                      classProbs = FALSE,
#                      trim = TRUE)

# to use lm
# algorithm <- "lm"
# ctrl <- trainControl(method = "none",
#                      verboseIter = TRUE,
#                      allowParallel = TRUE,
#                      returnData = FALSE,
#                      returnResamp = 'none',
#                      savePredictions = FALSE,
#                      classProbs = FALSE,
#                      trim = TRUE)

# # create and register processor cluster
# clust <- parallel::makePSOCKcluster(8)
# doParallel::registerDoParallel(clust)

# # train models and get predictions
# predictions <- sapply(1:numPrincipalComponents, 
#        function(PC){
#          fit <- 
#            train(method = algorithm,
#                  x = train$features[-train$heldOut,] %>% select(-sig_id),
#                  y = train$pca.outcomes$x[-train$heldOut,PC],
#                  trControl = ctrl)
#          predictions <- predict(fit, train$features[train$heldOut,] %>% select(-sig_id))
# })

# to use knn

numDimensions <- 3
algorithm <- "knn"
grid = data.frame(k = c(3,5,10,15,20))
ctrl <- trainControl(method = "LGOCV",
                     p = .8, number = 1,
                     verboseIter = TRUE,
                     allowParallel = TRUE,
                     returnData = FALSE,
                     returnResamp = 'none',
                     savePredictions = FALSE,
                     classProbs = FALSE,
                     trim = TRUE)


# train models for the numMechanisms most common mechanisms
numMechanisms <- 10
common_mechanisms <-
  train$targets %>% 
  pivot_longer(cols = - "sig_id",
               names_to = "mechanism",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  select(-value) %>% 
  group_by(mechanism) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(count)) %>% 
  top_n(n = numMechanisms, wt = count) %>% 
  pull(mechanism)

predictions <- sapply(1:numMechanisms, 
                      function(N){
                        fit <- 
                          train(method = algorithm,
                                # x = train$features[-train$heldOut,] %>% select(-sig_id),
                                x = train$pca.features$x[-train$heldOut, 1:numDimensions],                              # reduced dimension set for KNN
                                y = train$targets[-train$heldOut,common_mechanisms[N]],
                                tuneGrid = grid,
                                trControl = ctrl)
                        # predictions <- predict(fit, train$features[train$heldOut,] %>% select(-sig_id))
                        predictions <- predict(fit, train$pca.features$x[train$heldOut, 1:numDimensions])               # predictions for KNN
                        
                      })

# # stop cluster
# parallel::stopCluster(clust)
# doParallel::stopImplicitCluster()

# beep!
beepr::beep(5)

# Measure performance --------------------------------------------------------------------------------------------------

# target, approximate and predicted outcomes
y <-
  train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix()

y_app <- 
  train$pca.outcomes$x[train$heldOut,1:numPrincipalComponents] %*% 
  t(train$pca.outcomes$rotation[,1:numPrincipalComponents])

y_hat <- predictions

# y_hat <- 
#   cbind(predict(fit005_001, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_002, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_003, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_004, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_005, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_006, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_007, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_008, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_009, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1),
#         predict(fit005_010, train$features[train$heldOut,] %>% select(-sig_id)) %>% matrix(ncol = 1)) %>%  
#   magrittr::multiply_by_matrix(t(train$pca.outcomes$rotation[,1:numPrincipalComponents]))

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

# Measure performance for models created to predict most common mechanisms ---------------------------------------------

sapply(1:10, function(N){
  RMSE(predictions[,N],
     train$targets[train$heldOut,common_mechanisms[N]])
})


df <- data.frame(sample = train$heldOut[,1],
                 predictions = predictions[,1], 
                 target = train$targets[train$heldOut,common_mechanisms[1]])

df %>% 
  filter(target == 1)

df %>% 
  # filter(target == 1) %>% 
  summarise(score = score(reference = data.frame(target),
                          prediction = data.frame(predictions)))

score(reference = df %>% select(target),
      prediction = df %>% select(predictions) %>% mutate(predictions = limitRange(predictions)))

df %>% 
  ggplot(aes(x = target %>% as.factor(), y = predictions)) +
  geom_boxplot()


# Save models ----------------------------------------------------------------------------------------------------------

# save(list = paste("fit005",1:numPrincipalComponents %>% as.character() %>% str_pad(width=3,pad="0"), sep = "_"),
#      file = "ModelFit005.Rdata")


# Analysis -------------------------------------------------------------------------------------------------------------

# confusion matrices

cm_target_app <- confusionMatrix(data = y_app %>% limitRange %>% round(0) %>% as.factor,
                                 reference = y %>% as.factor,
                                 positive = "1")

cm_app_prediction <- confusionMatrix(data = y_hat %>% limitRange %>% round(0) %>% as.factor,
                                     reference = y_app %>% limitRange %>% round(0) %>% as.factor,
                                     positive = "1")

cm_target_prediction <- confusionMatrix(data = y_hat %>% limitRange %>% round(0) %>% as.factor,
                                        reference = y %>% as.factor,
                                        positive = "1")


cm_target_app
cm_app_prediction
cm_target_prediction
# 
# # how does each model fit its target (vector of principal components) on the held out samples
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,1],
#            prediction = predict(fit005_001, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE001 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,2],
#            prediction = predict(fit005_002, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE002 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,3],
#            prediction = predict(fit005_003, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE003 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,4],
#            prediction = predict(fit005_004, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE004 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,5],
#            prediction = predict(fit005_005, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE005 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,6],
#            prediction = predict(fit005_006, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE006 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,7],
#            prediction = predict(fit005_007, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE007 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,8],
#            prediction = predict(fit005_008, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE008 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,9],
#            prediction = predict(fit005_009, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE009 = RMSE(reference,prediction))
# 
# data.frame(reference = train$pca.outcomes$x[train$heldOut,10],
#            prediction = predict(fit005_010, train$features[train$heldOut,] %>%
#                                   select(-sig_id)) %>% matrix(ncol = 1)) %>%
#   summarise(RMSE010 = RMSE(reference,prediction))
# 
# # is the prediction improving as we add more PCs?
# 
# sapply(1:numPrincipalComponents, function(PCs){
# 
#   Y <- train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix()
# 
#   Y_app <-
#     train$pca.outcomes$x[train$heldOut,1:PCs] %*%
#     t(train$pca.outcomes$rotation[,1:PCs])
# 
#   Y_hat <-
#     y_hat %*% train$pca.outcomes$rotation[,1:PCs] %*% t(train$pca.outcomes$rotation[,1:PCs])
# 
#   RMSE <- RMSE(Y,Y_hat)
# })
# 
# # histogram of PCs on target vectors
# train$pca.outcomes$x[-train$heldOut,1] %>% histogram
# train$pca.outcomes$x[-train$heldOut,2] %>% histogram
# train$pca.outcomes$x[-train$heldOut,3] %>% histogram
# train$pca.outcomes$x[-train$heldOut,4] %>% histogram
# train$pca.outcomes$x[-train$heldOut,5] %>% histogram
# train$pca.outcomes$x[-train$heldOut,6] %>% histogram
# train$pca.outcomes$x[-train$heldOut,7] %>% histogram
# train$pca.outcomes$x[-train$heldOut,8] %>% histogram
# train$pca.outcomes$x[-train$heldOut,9] %>% histogram
# train$pca.outcomes$x[-train$heldOut,10] %>% histogram
# 
# # variable importance for each model
# fit005_001 %>% varImp
# fit005_002 %>% varImp
# fit005_003 %>% varImp
# fit005_004 %>% varImp
# fit005_005 %>% varImp
# fit005_006 %>% varImp
# fit005_007 %>% varImp
# fit005_008 %>% varImp
# fit005_009 %>% varImp
# fit005_010 %>% varImp
# 
# # principal components rotation matrix
# # train$pca.outcomes$rotation[,1, drop = FALSE] %>%
# #   as.data.frame() %>%
# #   rownames_to_column("mechanism") %>%
# #   arrange(PC1)
