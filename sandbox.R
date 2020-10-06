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
  summarise(count = n(), .groups = 'drop')

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
