# Housekeeping -----------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("rpart")) install.packages("rpart") else library("rpart")
# if (!require("MASS")) install.packages("MASS") else library("MASS")
if (!require("data.table")) install.packages("data.table") else library("data.table")

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

# Basic EDA --------------------------------------------------------------------

# features and dimensions of the features object
train$features %>% names
train$features %>% glimpse
train$features %>% nrow

# glimpse at the categorical variables
train$features$sig_id %>% n_distinct()
train$features$cp_time %>% table()
train$features$cp_type %>% table()
train$features$cp_dose %>% table()

# is the order of observations compatible between features and targets?
identical(train$features$sig_id, train$targets$sig_id) %>% mean()

# combinations of the categorical variables
train$features %>% 
  group_by(cp_time, cp_type, cp_dose) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(count))

# NAs in the feature space?
train$features %>% 
  as.matrix() %>% 
  is.na() %>% 
  sum()

test$features %>% 
  as.matrix() %>% 
  is.na() %>% 
  sum()

# statistics on the features - gene expression
train$features %>%
  # rbind(train$features, test$features) %>% 
  summarise(across(contains("g."), 
                   list(mean = mean, sd = sd, min = min, max = max))) %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>% 
  separate(col = "variable",
           into = c("variable","statistic"),
           sep = "_") %>% 
  pivot_wider(names_from = "statistic",
              values_from = "value") %>% 
  arrange(variable) %>% 
  mutate(variable = as.factor(variable)) %>% 
  ggplot(aes(x = variable, y = mean)) +
  geom_point(alpha = .25) +
  geom_point(aes(y = min), color = 'red', size = 1, alpha = .25) +
  geom_point(aes(y = max), color = 'red', size = 1, alpha = .25) +
  labs(title = "Statistics for each of the gene expression features",
       x = "genes", y = "statistics",
       caption = "Mean, min. and max. gene expression values for each gene") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# statistics on the features - cell viability
train$features %>%
  # rbind(train$features, test$features) %>% 
  summarise(across(contains("c."), 
                   list(mean = mean, sd = sd, min = min, max = max))) %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>% 
  separate(col = "variable",
           into = c("variable","statistic"),
           sep = "_") %>% 
  pivot_wider(names_from = "statistic",
              values_from = "value") %>% 
  arrange(variable) %>% 
  mutate(variable = as.factor(variable)) %>% 
  ggplot(aes(x = variable, y = mean)) +
  geom_point(alpha = .5) +
  geom_point(aes(y = min), color = 'red', size = 1, alpha = .5) +
  geom_point(aes(y = max), color = 'red', size = 1, alpha = .5) +
  labs(title = "Statistics for each of the cell viability features",
       x = "genes", y = "statistics",
       caption = "Mean, min. and max. cell viability values for each gene") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# structure of the targets object
train$targets %>% glimpse

# mechanisms by frequency
train$targets %>% 
  pivot_longer(cols = - "sig_id",
               names_to = "mechanism",
               values_to = "value") %>% 
  filter(value == 1) %>% 
  select(-value) %>% 
  group_by(mechanism) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(count))

# how many mechanisms per sample?
train$targets %>% 
  column_to_rownames(var = "sig_id") %>% 
  as.matrix() %>% 
  rowSums() %>% 
  table()

# how many different outcomes?
outcome_combinations <-
  train$targets %>% 
  pivot_longer(cols = -"sig_id",
               names_to = "mechanism",
               values_to = "target") %>% 
  filter(target == 1) %>% 
  select(-target) %>% 
  arrange(sig_id,mechanism) %>%
  group_by(sig_id) %>% 
  mutate(number = row_number(),
         mechanisms = n()) %>% 
  pivot_wider(names_from = "number",
              values_from = "mechanism") %>% 
  unite(-"sig_id",-"mechanisms", col = "output", na.rm = TRUE, sep = " | ") %>% 
  group_by(output, mechanisms) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  arrange(desc(count))

outcome_combinations

# how many combinations have more than 1 mechanism?
outcome_combinations %>% 
  filter(mechanisms > 1) %>% 
  pull(output) %>% 
  unique() %>% 
  length()

rm(outcome_combinations)

# are mechanisms correlated?
train$targets %>% 
  column_to_rownames(var = "sig_id") %>% 
  as.matrix() %>% 
  cor() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mechanism1") %>% 
  pivot_longer(cols = -"mechanism1",
               names_to = "mechanism2",
               values_to = "correlation") %>% 
  filter(mechanism1 < mechanism2) %>% 
  arrange(desc(correlation))

# are features correlated?
correlations <-
  train$features %>% 
  column_to_rownames("sig_id") %>% 
  select(-c(cp_type, cp_time, cp_dose)) %>% 
  as.matrix() %>% 
  cor() %>% 
  as.data.frame() %>% 
  rownames_to_column("feature1") %>% 
  pivot_longer(cols = -"feature1",
               names_to = "feature2",
               values_to = "correlation") %>% 
  filter(feature1 < feature2) %>% 
  arrange(desc(correlation)) %>% 
  as.data.frame()

correlations %>% 
  filter(str_detect(feature1, "c.") & str_detect(feature2, "c.")) %>% 
  arrange(correlation) %>% 
  head(20)

correlations %>% 
  filter(str_detect(feature1, "g.") & str_detect(feature2, "g.")) %>% 
  arrange(desc(correlation)) %>% 
  head(20)

correlations %>% 
  filter(str_detect(feature1, "c.") & str_detect(feature2, "g.")) %>% 
  arrange(desc(correlation)) %>% 
  head(20)

rm(correlations)

# correlations between features and outcomes
correlations <-
  cor(train$features %>% select(contains("g.") | contains("c.")) %>% as.matrix,
      train$targets %>% select(-sig_id) %>% as.matrix) %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  pivot_longer(cols = -"predictor",
               names_to = "outcome",
               values_to = "correlation")

correlations %>% 
  as.data.frame() %>% 
  arrange(desc(correlation)) %>% 
  head(30)

correlations %>% 
  as.data.frame() %>% 
  arrange(correlation) %>% 
  head(30)

correlations %>% 
  ggplot(aes(x = correlation)) +
  geom_density() +
  scale_y_continuous(trans = 'sqrt')

rm(correlations)

# PCA on non categorical data
pc <- 
  train$features %>% 
  select(- c("sig_id", "cp_type", "cp_time", "cp_dose")) %>% 
  as.matrix() %>% 
  prcomp

var <- 
  data.frame(pc = 1:length(pc$sdev),
             sd = pc$sdev,
             prop = cumsum(pc$sdev)/sum(pc$sdev))

var %>% 
  ggplot(aes(x = pc, y = prop)) +
  geom_line() + 
  geom_point(data = var %>% filter(prop >= .7) %>% slice_min(prop, 1), 
             color = 'red', size = 3) +
  geom_point(data = var %>% filter(prop >= .9) %>% slice_min(prop, 1), 
             color = 'red', size = 3) +
  labs(title = "Cumulative standard proportion of standard deviation per principal component",
       x = "number of principal components", y = "proportion of total standard deviation",
       caption = "cumulative proportion of total standard deviation for each principal component, with red points at 70% (441) and 90% (704)") +
  scale_y_continuous(labels = scales::percent)

pc.df <- 
  pc$x[,1:2] %>% 
  as.data.frame() %>% 
  cbind(train$targets$sig_id) %>% 
  select(sig_id = `train$targets$sig_id`, PC1, PC2)

train$targets %>% 
  pivot_longer(cols = -sig_id, values_to = "score", names_to = "mechanism") %>% 
  filter(score == 1) %>% 
  select(-score) %>% 
  left_join(pc.df, by = "sig_id") %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

rm(var,pc.df)

# how do the principal components correlate with each mechanism

correlations <-
  cor(pc$x %>% as.matrix,
      train$targets %>% select(-sig_id) %>% as.matrix) %>% 
  as.data.frame() %>% 
  rownames_to_column("PC") %>% 
  pivot_longer(cols = -"PC", values_to = "correlation", names_to = "mechanism")

correlations %>% 
  arrange(desc(correlation)) %>% 
  data.table()

correlations %>% 
  arrange(correlation)

rm(correlations)


# Proof of concept - predict most common mechanism -----------------------------

grid <- expand.grid(predFixed = c(10,50,100,200),
                    minNode = 2)

fit_Rborist <- 
  train(method = "Rborist",
        x = select(train$features, -sig_id),
        y = train$targets$nfkb_inhibitor,
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9,
                                 number = 1,
                                 verboseIter = TRUE))

# variable importance
fit_Rborist %>% varImp
fit_Rborist %>% varImp %>% plot

# OK vs NOK
data.frame(target = train$targets$nfkb_inhibitor,
           prediction = predict(fit_Rborist, select(train$features, -sig_id))) %>% 
  summarise(OK = sum(target == round(prediction,0)),
            NOK = sum(target != round(prediction,0)))
           
# how does it miss?
data.frame(target = train$targets$nfkb_inhibitor,
           prediction = predict(fit_Rborist, select(train$features, -sig_id))) %>% 
  filter(target != round(prediction,0))

# Predict most common mechanism with PCA ---------------------------------------

grid <- expand.grid(predFixed = c(10,20,30,40,60,80),
                    minNode = 2)

fit_Rborist_PCA <- 
  train(method = "Rborist",
        x = pc$x[,1:100],
        y = train$targets$nfkb_inhibitor,
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9,
                                 number = 1,
                                 verboseIter = TRUE))

# variable importance
fit_Rborist_PCA %>% varImp
fit_Rborist_PCA %>% varImp %>% plot

# OK vs NOK
data.frame(target = train$targets$nfkb_inhibitor,
           prediction = predict(fit_Rborist, select(train$features, -sig_id))) %>% 
  summarise(OK = sum(target == round(prediction,0)),
            NOK = sum(target != round(prediction,0)))

# how does it miss?
data.frame(target = train$targets$nfkb_inhibitor,
           prediction = predict(fit_Rborist, select(train$features, -sig_id))) %>% 
  filter(target != round(prediction,0))

histogram(predict(fit_Rborist, select(train$features, -sig_id)))

# confusion matrix
confusionMatrix(data = 
                  predict(fit_Rborist, select(train$features, -sig_id)) %>% 
                  round(0) %>% as.factor(),
                reference = train$targets$nfkb_inhibitor %>% as.factor(),
                positive = "1")

# Experiment with number of PCAs -----------------------------------------------

grid <- expand.grid(predFixed = c(40,60,80,120,140,160),
                    minNode = 2)

fit_Rborist_PCA_200 <- 
  train(method = "Rborist",
        x = pc$x[,1:200],
        y = train$targets$nfkb_inhibitor,
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9,
                                 number = 1,
                                 verboseIter = TRUE))

grid <- expand.grid(predFixed = c(40,80,120,160,200,240),
                    minNode = 2)

fit_Rborist_PCA_400 <- 
  train(method = "Rborist",
        x = pc$x[,1:400],
        y = train$targets$nfkb_inhibitor,
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9,
                                 number = 1,
                                 verboseIter = TRUE))

# Create model with PCAs applied to features and outcomes ----------------------

grid <- expand.grid(predFixed = c(40,60,80,120,140,160,200),
                    minNode = 2)

numPrincipalComponents <- 872

fit_Rborist_all <-
  train(method = "Rborist",
        x = cbind(select(train$features, cp_type, cp_dose, cp_time),
                  train$pca.features$x[,1:100]),
        y = train$pca.outcomes$x[,1],
        tuneGrid = grid,
        trControl = trainControl(method = "LGOCV", p = .9, number = 1,
                                 verboseIter = TRUE))

predictions <- 
  predict(fit_Rborist_all,
          cbind(select(train$features, cp_type, cp_dose, cp_time),
                train$pca.features$x[,1:100])) %>% 
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

# Examine results from model 003 ---------------------------------------------------------------------------------------

# evaluate results on held out samples

# targets
y <- train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix()

# PCA approximation
y_app <- 
  train$pca.outcomes$x[train$heldOut,1, drop = FALSE] %*% t(train$pca.outcomes$rotation[,1])

# predictions
y_hat <- 
  predict(fit003_200,
          train$features[train$heldOut,] %>% select(-sig_id)) %>% 
  matrix(nrow = 1) %>% 
  t() %*% 
  t(train$pca.outcomes$rotation[,1])

# scores

limitRange <- function(x,delta){pmin(pmax(x,delta),1-delta)}

score <- 
  function(reference,prediction){
    N <- nrow(reference)
    M <- ncol(reference)
    score <- - sum(reference*log(prediction) + (1-reference)*log(1-prediction)) / (N*M)
}

delta <- 10^(-15)

score_target_app <-
  score(reference = y,
        prediction = limitRange(y_app,delta))

score_app_prediction <-
  score(reference = limitRange(y_app,delta),
        prediction = limitRange(y_hat,delta))

score_target_prediction <-
  score(reference = y,
        prediction = limitRange(y_hat,delta))

score_target_app
score_app_prediction
score_target_prediction

cm_target_app <- confusionMatrix(reference = y %>% as.factor(),
                                 data = y_app %>% limitRange(0) %>% round(0) %>% as.factor(),
                                 positive = "1")

cm_app_prediction <- confusionMatrix(reference = y_app %>% limitRange(0) %>% round(0) %>% as.factor(),
                                     data = y_hat %>% limitRange(0) %>% round(0) %>% as.factor(),
                                     positive = "1")

cm_target_prediction <- confusionMatrix(reference = y %>% as.factor(),
                                        data = y_app %>% limitRange(0) %>% round(0) %>% as.factor(),
                                        positive = "1")
cm_target_app
cm_app_prediction
cm_target_prediction

# Experiment with options in the caret package -------------------------------------------------------------------------

N = 500

x <- matrix(runif(N * 3,min=0,max=1),
            ncol=3)

colnames(x) <- c("x1","x2","x3")

epsilon <- matrix(rnorm(N,0,.1),ncol=1) 

y <- 
  x[,1,drop=FALSE] + 
  2 * x[,2,drop=FALSE] + 
  x[,3,drop=FALSE]^2 + 
  epsilon

colnames(y) <- "outcome"

fit <- train(x = x,
             y = y[,1],
             method = "Rborist",
             ntree = 50,
             tuneGrid = data.frame(predFixed = 1, minNode = 2),
             trControl = trainControl(method = 'none',
                                      returnData = FALSE,
                                      returnResamp = 'none',
                                      savePredictions = FALSE,
                                      classProbs = FALSE,
                                      trim = TRUE,
                                      allowParallel = TRUE))

fit
object.size(fit) / 2^20
object.size(fit$finalModel) / 2^20

# Naive predictions ----------------------------------------------------------------------------------------------------

# predict global mean
prediction <- train$targets %>% select(-sig_id) %>%  as.matrix() %>% mean()
target <- train$targets %>% select(-sig_id) %>% as.matrix()
score(reference = target, prediction = prediction)

# predict mean for each mechanism
prediction <-
  target %>% 
  colMeans() %>% 
  matrix(nrow = nrow(target), ncol = ncol(target), byrow = TRUE)
colnames(prediction) <- colnames(target)
score(reference = target, prediction = prediction)

# how good can the score get with perfect predictions for one mechanism and the mean for the others?
prediction[,"nfkb_inhibitor"] <- target[,"nfkb_inhibitor"]
score(reference = target, prediction = prediction)


# Examine most popular mechanism in details ----------------------------------------------------------------------------

numMechanisms <- 1
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

x <- train$features
y <- train$targets[,'nfkb_inhibitor']

df <-
  data.frame(target = y,
             prediction = mean(y))
  
df %>% 
  summarise(score = score(reference = target %>% as.data.frame(), 
                          prediction = prediction %>% as.data.frame()))

# data frame in long format
df <-
  train$targets %>% 
  select(sig_id, nfkb_inhibitor) %>% 
  mutate(nfkb_inhibitor = as.factor(nfkb_inhibitor)) %>% 
  left_join(train$features, by = "sig_id") %>% 
  pivot_longer(cols = -c(sig_id,nfkb_inhibitor), names_to = "feature", values_to = "activation")

# mean activation for each of the features, grouped by output
df %>% 
  group_by(nfkb_inhibitor, feature) %>% 
  summarise(mean_activation = mean(activation)) %>% 
  ggplot(aes(x = feature, y = mean_activation, color = nfkb_inhibitor)) +
  geom_point()

# 

df %>% 
  filter(feature %in% c("g.0","g.1")) %>% 
  pivot_wider(names_from = "feature", values_from = "activation") %>% 
  ggplot(aes(x = g.0, y = g.1, color = nfkb_inhibitor)) +
  geom_point()

# Fit different models and compare results -----------------------------------------------------------------------------

# classification tree

fitTree <-
  train(method = "rpart",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$targets$nfkb_inhibitor[-train$heldOut] %>% as.factor(),
        trControl = trainControl(method = 'LGOCV', p = .9, number = 1,
                                 verboseIter = TRUE))

# random forest

fitForest <-
  train(method = "Rborist",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$targets$nfkb_inhibitor[-train$heldOut] %>% as.factor(),
        verbose = TRUE,
        tuneGrid = data.frame(predFixed = 50, minNode = 2),
        trControl = trainControl(method = 'LGOCV', p = .9, number = 1,
                                 verboseIter = TRUE))
  
# LDA

fitLDA <-
  train(method = "lda",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$targets$nfkb_inhibitor[-train$heldOut] %>% as.factor(),
        trControl = trainControl(method = 'LGOCV', p = .9, number = 1,
                                 verboseIter = TRUE))
# generalized linear model

fitGLM <-
  train(method = "glm",
        x = train$features[-train$heldOut,] %>% select(-sig_id),
        y = train$targets$nfkb_inhibitor[-train$heldOut] %>% as.factor(),
        trControl = trainControl(method = 'LGOCV', p = .9, number = 1,
                                 verboseIter = TRUE))

# KNN

fitKNN <-
  train(method = "knn",
        x = train$pca.features$x[-train$heldOut,1:3],
        y = train$targets$nfkb_inhibitor[-train$heldOut] %>% as.factor(),
        tuneGrid = data.frame(k = c(3,5,7,10)),
        trControl = trainControl(method = 'LGOCV', p = .9, number = 1,
                                 verboseIter = TRUE))

# confusion matrices

confusionMatrix(fitTree, positive = "1")
confusionMatrix(fitForest, positive = "1")
confusionMatrix(fitLDA, positive = "1")
confusionMatrix(fitGLM, positive = "1")
confusionMatrix(fitKNN, positive = "1")

# examine models

predictions <- 
  data.frame(groundTruth = train$targets$nfkb_inhibitor[train$heldOut] %>% as.factor,
             predictionTree = predict(fitTree, train$features[train$heldOut,] %>% select(-sig_id)),
             predictionForest = predict(fitForest, train$features[train$heldOut,] %>% select(-sig_id)),
             predictionLDA = predict(fitLDA, train$features[train$heldOut,] %>% select(-sig_id)),
             predictionGLM = predict(fitGLM, train$features[train$heldOut,] %>% select(-sig_id)),
             predictionKNN = predict(fitKNN, train$pca.features$x[train$heldOut,1:3]))

predictions %>% 
  filter(groundTruth == 1)

predictions %>% 
  filter(groundTruth != predictionForest)

predictions %>% 
  mutate(groundTruth = groundTruth %>% as.character %>% as.integer(),
         predictionTree = predictionTree %>% as.character %>% as.integer(),
         predictionForest = predictionForest %>% as.character %>% as.integer(),
         predictionLDA = predictionLDA %>% as.character %>% as.integer(),
         predictionGLM = predictionGLM %>% as.character %>% as.integer(),
         predictionKNN = predictionKNN %>% as.character %>% as.integer()) %>% 
  mutate(mean = (predictionTree+predictionForest+predictionLDA+predictionGLM+predictionKNN)/5) %>% 
  summarise(RMSE_tree = RMSE(groundTruth,predictionTree),
            RMSE_forest = RMSE(groundTruth,predictionForest),
            RMSE_LDA = RMSE(groundTruth,predictionLDA),
            RMSE_GLM = RMSE(groundTruth,predictionGLM),
            RMSE_KNN = RMSE(groundTruth,predictionKNN),
            RMSE_mean = RMSE(groundTruth,mean))
  

# Correlation between predictors and outcome ---------------------------------------------------------------------------

x <- train$features[-train$heldOut,] %>% select(-sig_id) %>% as.matrix()
pc <- train$pca.features$x[-train$heldOut,] %>% as.matrix()
y <- train$targets$nfkb_inhibitor[-train$heldOut] %>% as.matrix()

features_correlations <-
  data.frame(feature = colnames(x),
             correlation = cor(x,y)) %>% 
  arrange(desc(correlation))

pc.correlations <-
  data.frame(PC = colnames(pc),
             correlation = cor(pc, y)) %>% 
  arrange(desc(abs(correlation)))

highest.correlations <-
  pc.correlations %>% 
  head(3) %>% 
  .$PC

fit_KNN_high_cor <- train(method = "knn",
                          x = train$pca.features$x[-train$heldOut,highest.correlations],
                          y = train$targets$nfkb_inhibitor[-train$heldOut] %>% as.factor(),
                          preProcess = preProcess(method = c("center","scale")),
                          tuneGrid = data.frame(k = c(3,5,7,10)),
                          trControl = trainControl(method = 'LGOCV', p = .9, number = 1,
                                                   verboseIter = TRUE))

confusionMatrix(fit_KNN_high_cor, positive = "1")
