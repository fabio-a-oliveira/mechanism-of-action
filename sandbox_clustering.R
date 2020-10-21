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

# k-means clustering - influence of number of centers on ss ratio -----------------------------------------

numMechanisms <- 206

numCenters <- seq(1,numMechanisms,10)

ssRatio <- sapply(numCenters, function(numCenters){
  print(paste("Number of centers:", numCenters))
  fit <- train$features %>% 
  column_to_rownames("sig_id") %>% 
  kmeans(centers = numCenters,iter.max = 100)

fit$betweenss/fit$totss
})

plot(numCenters, ssRatio)

# do clusters capture groups with different prevalence ?-------------------------------------------

numMechanisms <- 206

numCenters <- seq(1,numMechanisms,10)

std <- sapply(numCenters, function(numCenters){
  
  print(paste("Number of centers:", numCenters))
  
  fit <- train$features %>% 
    column_to_rownames("sig_id") %>% 
    kmeans(centers = numCenters,iter.max = 100)
  
  clusters <-
    fit$cluster %>% 
    as.data.frame %>% 
    rownames_to_column("sig_id") %>% 
    rename(cluster = ".") %>% 
    mutate(cluster = as.ordered(cluster)) %>% 
    as.tibble
  
  outcomes <-
    train$outcomes %>% 
    pivot_longer(cols = -"sig_id",
                 names_to = "mechanism",
                 values_to = "target") %>% 
    left_join(clusters, by = "sig_id")
  
  outcomes %>% 
    group_by(mechanism,cluster) %>% 
    summarise(prevalence = sum(target)/n(), .groups = 'drop') %>% 
    group_by(mechanism) %>% 
    summarise(std = sd(prevalence), .groups = 'drop') %>% 
    arrange(desc(std)) %>% 
    pull(std)
  
})

std <- t(std)
colnames(std) <- 
  colnames(train$outcomes %>% select(-sig_id))
rownames(std) <- numCenters
std <- as.data.frame(std)

std %>% 
  rownames_to_column("numCenters") %>% 
  pivot_longer(cols = -"numCenters",
               names_to = "mechanism",
               values_to = "SD") %>% 
  filter(!is.na(SD)) %>% 
  mutate(numCenters = as.integer(numCenters)) %>% 
  ggplot(aes(x = numCenters %>% as.factor, y = SD)) +
  geom_boxplot()

# Use cluster target means as predictions -------------------------------------------------------

numMechanisms <- 206
numCenters <- 
  c(1,2,3,5,seq(10,100,10),seq(120,200,20),206)

results <- sapply(numCenters, function(numCenters){
  
  # fit model - find cluster centers
  fit <- train$features %>% 
    column_to_rownames("sig_id") %>% 
    kmeans(centers = numCenters,iter.max = 100)
  
  # prevalence in the training set
  prevalence <-
    train$outcomes %>% 
    mutate(cluster = fit$cluster) %>% 
    pivot_longer(cols = c(-"sig_id",-"cluster"),
                 names_to = "mechanism",
                 values_to = "target") %>% 
    group_by(mechanism,cluster) %>% 
    summarise(prevalence = sum(target)/n(), .groups = 'drop')
  
  # clusters in the validation set
  clusters <- sapply(1:nrow(validate$features), function(row){
    squared.distances <- 
      rowSums((fit$centers - as.integer(validate$features[row,-1]))^2)
    cluster <- names(squared.distances)[which.min(squared.distances)]
  })
  
  # targets and predictions on outcome matrix
  outcomes <- 
    validate$outcomes %>% 
    mutate(cluster = clusters %>% as.integer) %>% 
    pivot_longer(cols = c(-"sig_id",-"cluster"),
                 names_to = "mechanism", values_to = "target") %>% 
    left_join(prevalence, by = c("mechanism","cluster"))
  
  # score on the validation set
  outcomes %>% 
    mutate(score = - target * log(limitRange(prevalence, delta, 1-delta)) 
           - (1-target) * log(1 - limitRange(prevalence, delta, 1-delta))) %>% 
    summarise(logLoss = mean(score),
              RMSE = RMSE(prevalence,target))

})

results.df <- data.frame(logLoss = results[1,],
                         RMSE = results[2,])

results
