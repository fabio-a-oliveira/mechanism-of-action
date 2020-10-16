# Housekeeping -----------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("data.table")) install.packages("data.table") else library("data.table")

options("datatable.print.topn" = 20)

# Load data --------------------------------------------------------------------

train <- 
  list(features = read.csv(file.path("files","train_features.csv")),
       targets = read.csv(file.path("files","train_targets_scored.csv")),
       nonscored = read.csv(file.path("files","train_targets_nonscored.csv")))

test <- 
  list(features = read.csv(file.path("files","test_features.csv")))

# Transform categorical data into factors --------------------------------------

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

# Perform PCA on matrix of outcomes --------------------------------------------

pc <- 
  train$targets %>% 
  select(-sig_id) %>% 
  prcomp()

N <- 206
N <- 6

approximation <-  pc$x[,1:N] %*% t(pc$rotation[,1:N])

approximation %>% as.vector %>% cut(breaks = seq(0,1,.1), include.lowest = TRUE) %>% table()

((approximation %>% round(0)) == as.matrix(select(train$targets,-sig_id))) %>% 
  mean()

# Performance function ---------------------------------------------------------

# p <- pmax(pmin(round(approximation,0),1−10^−15),10^−15)
p <- pmax(pmin(approximation,1−10^−15),10^−15)
y <- train$targets %>% select(-sig_id) %>% as.matrix()

score <- - sum(y * log(p) + (1-y) * log(1 - p)) / (206 * 23814)

score
