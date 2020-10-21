# Header -----------------------------------------------------------------------------------------

# Model 012

# Application of PCA on the centroids for each class (inspired by reduced rank LDA, ESL 4.3.3)
# Ensemble of LDA and KNN, both fit on subspace of PCA on centroids

# Housekeeping ---------------------------------------------------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("rpart")) install.packages("rpart") else library("rpart")
# if (!require("MASS")) install.packages("MASS") else library("MASS")
if (!require("data.table")) install.packages("data.table") else library("data.table")
if (!require("foba")) install.packages("foba") else library("foba")
if (!require("pls")) install.packages("pls") else library("pls")
if (!require("matrixStats")) install.packages("matrixStats") else library("matrixStats")

options("datatable.print.topn" = 20)

# Declaration of some useful functions ---------------------------------------------------------------------------------

delta <- 10^(-15)

limitRange <- function(x,xmin,xmax){
  pmin(pmax(x,xmin),xmax)}

score <- 
  function(reference,prediction){
    N <- length(reference)
    M <- length(reference)
    delta <- 10^(-15)
    - sum(reference*log(limitRange(prediction,delta,1-delta)) + (1-reference)*log(1-limitRange(prediction,delta,1-delta))) / (N*M)
  }

# Load and prepare data ------------------------------------------------------------------------------------------------

# read files
files <-
  list(train_features = read_csv(file.path("files","train_features.csv")),
       train_targets_scored = read_csv(file.path("files","train_targets_scored.csv")),
       train_targets_nonscored = read_csv(file.path("files","train_targets_nonscored.csv")),
       test_features = read_csv(file.path("files","test_features.csv")))

# create train, test and validate sets
proportionHeldOut <- .2
heldOutSamples <- 
  sample_frac(files$train_features, 
              size = proportionHeldOut,
              replace = FALSE) %>% 
  pull(sig_id)

train <-
  list(features = 
         files$train_features %>%
         filter(!(sig_id %in% heldOutSamples)) %>% 
         filter(cp_type != "ctl_vehicle") %>% 
         select(-cp_type,-cp_time,-cp_dose) %>% 
         mutate(sig_id = as.factor(sig_id)),
       outcomes = 
         files$train_targets_scored %>% 
         filter(!(sig_id %in% heldOutSamples)) %>% 
         mutate(sig_id = as.factor(sig_id)))

train$outcomes <- 
  files$train_targets_scored %>% 
  filter(sig_id %in% train$features$sig_id) %>% 
  mutate(sig_id = as.factor(sig_id))
train$pca.features <- train$features %>% column_to_rownames("sig_id") %>% prcomp()
train$pca.outcomes <- train$outcomes %>% column_to_rownames("sig_id") %>% prcomp()

validate <-
  list(features = 
         files$train_features %>%
         filter(sig_id %in% heldOutSamples) %>% 
         filter(cp_type != "ctl_vehicle") %>% 
         select(-cp_type,-cp_time,-cp_dose) %>% 
         mutate(sig_id = as.factor(sig_id)))
validate$outcomes <- 
  files$train_targets_scored %>% 
  filter(sig_id %in% validate$features$sig_id) %>% 
  mutate(sig_id = as.factor(sig_id))
validate$pca.features <- validate$features %>% column_to_rownames("sig_id") %>% prcomp()

test <-
  list(features = 
         files$test_features %>% 
         filter(cp_type != "ctl_vehicle") %>% 
         select(-cp_type,-cp_time,-cp_dose) %>% 
         mutate(sig_id = as.factor(sig_id)))
test$pca.features <- test$features %>% select(-sig_id) %>% prcomp()

# Find centroids and define rotation matrix ---------------------------------------------------

centroids <- 
  train$outcomes %>% 
  pivot_longer(cols = -"sig_id",
               names_to = "mechanism",
               values_to = "target") %>% 
  filter(target == 1) %>% 
  select(-target) %>% 
  left_join(train$features, by = "sig_id") %>% 
  group_by(mechanism) %>% 
  summarise(across(-sig_id, mean))

pca.centroids <-
  centroids %>% 
  as.data.frame() %>% 
  column_to_rownames("mechanism") %>% 
  as.matrix %>% 
  prcomp()

pca.centroids %>% summary

pca.centroids$rotation

# plot centroids on PC1/PC2 plane
pca.centroids$x[,1:2] %>% 
  as.data.frame %>% 
  rownames_to_column("mechanism") %>% 
  ggplot(aes(x = PC1, y = PC2, color = mechanism)) +
  geom_point() +
  theme(legend.position = 'none')

# plot centroids on PC3/PC4 plane
pca.centroids$x[,3:4] %>% 
  as.data.frame %>% 
  rownames_to_column("mechanism") %>% 
  ggplot(aes(x = PC3, y = PC4, color = mechanism)) +
  geom_point() +
  theme(legend.position = 'none')

# Create LDA model using PCA decomposition -------------------------------------------------------

numPrincipalComponents <- 10

# x and y training vectors
y.train <- 
  train$outcomes %>% 
  pivot_longer(cols = -"sig_id",
               names_to = "mechanism",
               values_to = "target") %>% 
  filter(target == 1) %>% 
  select(-target) %>% 
  mutate(mechanism = as.factor(mechanism))

x.train <-
  left_join(y.train, train$features, by = "sig_id") %>% 
  select(-sig_id,-mechanism) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(pca.centroids$rotation[,1:numPrincipalComponents])

y.train <- y.train %>% select(-sig_id)

# model training
fit_LDA <- train(method = 'lda',
             x = x.train, 
             y = y.train$mechanism,
             trControl = trainControl(method = 'none',verboseIter = TRUE))

# x and y validation vectors
y.validation <- 
  validate$outcomes

x.validation <-
  left_join(y.validation %>% select(sig_id), 
            validate$features, by = "sig_id") %>% 
  select(-sig_id) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(pca.centroids$rotation[,1:numPrincipalComponents])

# make predictions on validation set
y.hat <- 
  left_join(y.validation %>% select(sig_id),
             predict(fit_LDA,x.validation,type='prob') %>% as.data.frame %>% mutate(sig_id = y.validation$sig_id),
            by = "sig_id")

# evaluate performance
y.validation.long <-
  y.validation %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target")

y.hat.long <-
  y.hat %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction")

results <-
  left_join(y.validation.long,
            y.hat.long,
            by = c("sig_id","mechanism"))

results %>% 
  replace_na(list(prediction = 0)) %>% 
  mutate(predictionLimited = limitRange(prediction,delta,1-delta)) %>% 
  summarise(logLoss = mean(- target * log(predictionLimited) - 
                             (1-target) * log(1-predictionLimited)),
            RMSE = RMSE(target,prediction))

# Create KNN model using PCA decomposition -------------------------------------------------------

numPrincipalComponents <- 50

# x and y training vectors
y.train <- 
  train$outcomes %>% 
  pivot_longer(cols = -"sig_id",
               names_to = "mechanism",
               values_to = "target") %>% 
  filter(target == 1) %>% 
  select(-target) %>% 
  mutate(mechanism = as.factor(mechanism))

x.train <-
  left_join(y.train, train$features, by = "sig_id") %>% 
  select(-sig_id,-mechanism) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(pca.centroids$rotation[,1:numPrincipalComponents])

y.train <- y.train %>% select(-sig_id)

# model training
fit_KNN <- train(method = 'knn',
                 tuneGrid = data.frame(k = 400),
                 x = x.train, 
                 y = y.train$mechanism,
                 trControl = trainControl(method = 'none',verboseIter = TRUE))

# x and y validation vectors
y.validation <- 
  validate$outcomes

x.validation <-
  left_join(y.validation %>% select(sig_id), 
            validate$features, by = "sig_id") %>% 
  select(-sig_id) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(pca.centroids$rotation[,1:numPrincipalComponents])

# make predictions on validation set
y.hat <- 
  left_join(y.validation %>% select(sig_id),
            predict(fit_KNN,x.validation,type='prob') %>% as.data.frame %>% mutate(sig_id = y.validation$sig_id),
            by = "sig_id")

# evaluate performance
y.validation.long <-
  y.validation %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target")

y.hat.long <-
  y.hat %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction")

results <-
  left_join(y.validation.long,
            y.hat.long,
            by = c("sig_id","mechanism"))

results %>% 
  replace_na(list(prediction = 0)) %>% 
  mutate(predictionLimited = limitRange(prediction,delta,1-delta)) %>% 
  summarise(logLoss = mean(- target * log(predictionLimited) - 
                             (1-target) * log(1-predictionLimited)),
            RMSE = RMSE(target,prediction))

# Make ensemble predictions ----------------------------------------------------------------


y.validation.long <-
  validate$outcomes %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target")

# y.hat LDA

numPrincipalComponents <- 10 # same as used above on model creation

x.validation.lda <-
  left_join(validate$outcomes %>% select(sig_id), 
            validate$features, by = "sig_id") %>% 
  select(-sig_id) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(pca.centroids$rotation[,1:numPrincipalComponents])

y.hat.lda.long <-
  left_join(y.validation %>% select(sig_id),
            predict(fit_LDA,x.validation.lda,type='prob') %>% as.data.frame %>% mutate(sig_id = y.validation$sig_id),
            by = "sig_id") %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction_lda")

# y.hat KNN

numPrincipalComponents <- 50 # same as used above on model creation

x.validation.knn <-
  left_join(validate$outcomes %>% select(sig_id), 
            validate$features, by = "sig_id") %>% 
  select(-sig_id) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(pca.centroids$rotation[,1:numPrincipalComponents])

y.hat.knn.long <-
  left_join(y.validation %>% select(sig_id),
            predict(fit_KNN,x.validation.knn,type='prob') %>% as.data.frame %>% mutate(sig_id = y.validation$sig_id),
            by = "sig_id") %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction_knn")

# combine results

results <- y.validation.long %>% 
  left_join(y.hat.lda.long, by = c("sig_id","mechanism")) %>% 
  left_join(y.hat.knn.long, by = c("sig_id", "mechanism")) %>% 
  mutate(prediction = (prediction_knn + prediction_lda)/2) %>% 
  replace_na(list(prediction = 0)) %>% 
  mutate(prediction = limitRange(prediction,delta,1-delta))

# calculate logLoss
results %>% 
  summarise(logLoss = mean(- target * log(prediction) - (1-target)*log(1-prediction)))

# plot results on predictions plane
results %>% 
  filter(mechanism == "nfkb_inhibitor") %>% 
  sample_n(size = 1000) %>% 
  ggplot(aes(x = prediction_lda, y = prediction_knn, color = as.factor(target))) +
  geom_point()

