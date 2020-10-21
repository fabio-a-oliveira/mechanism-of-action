# Model009.R

# GLM model on selection of PCs, with range limiting
# competition submission 001 - 2020-10-14

# Housekeeping ---------------------------------------------------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("rpart")) install.packages("rpart") else library("rpart")
# if (!require("MASS")) install.packages("MASS") else library("MASS")
if (!require("data.table")) install.packages("data.table") else library("data.table")
if (!require("kknn")) install.packages("kknn") else library("kknn")

options("datatable.print.topn" = 20)
options("max.print" = 50)

# Declaration of some useful functions ---------------------------------------------------------------------------------

delta <- 10^(-15)

limitRange <- function(x,xmin,xmax){
  pmin(pmax(x,xmin),xmax)}

score <- 
  function(reference,prediction){
    N <- nrow(reference)
    M <- ncol(reference)
    delta <- 10^(-15)
    - sum(reference*log(limitRange(prediction,delta,1-delta)) + 
            (1-reference)*log(1-limitRange(prediction,delta,1-delta))) /  (N*M)
  }

# Load data ------------------------------------------------------------------------------------------------------------

train <- 
  list(features = read_csv(file.path("files","train_features.csv")),
       targets = read_csv(file.path("files","train_targets_scored.csv")),
       nonscored = read_csv(file.path("files","train_targets_nonscored.csv")))

test <- 
  list(features = read_csv(file.path("files","test_features.csv")))

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

# Create predictions object --------------------------------------------------------------------------------------------

proportionHeldOut <- .1
train$heldOut <- sample(x = 1:nrow(train$features),
                        size = round(nrow(train$targets) * proportionHeldOut,0),
                        replace = FALSE)

predictions <- matrix(0, 
                      nrow = nrow(train$features) + nrow(test$features),
                      ncol = ncol(train$targets) - 1)
rownames(predictions) <- c(train$features$sig_id %>% as.character, test$features$sig_id %>% as.character())
colnames(predictions) <- colnames(train$targets %>% select(-sig_id))

# Run Principal Component Analysis -------------------------------------------------------------------------------------

train$pca.outcomes <-
  train$targets[-train$heldOut,] %>% 
  select(-sig_id) %>% 
  prcomp()

train$pca.features <-
  train$features[-train$heldOut,] %>% 
  select(-c(sig_id,cp_type,cp_time,cp_dose)) %>% 
  prcomp()

# Make predictions for each of the mechanisms --------------------------------------------------------------------------

# parameters
numMechanisms <- ncol(train$targets) - 1
method <- "glm"
numPrincipalComponents <- 50

# full set of data for which to calculate model performance
new.data <- 
  rbind(train$features %>% select(-cp_type,-cp_time,-cp_dose),
        test$features %>% select(-cp_type,-cp_time,-cp_dose)) %>% 
  column_to_rownames("sig_id") %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(train$pca.features$rotation)

# observations that are automatically 0
ctl_vehicle <- c(train$features %>% filter(cp_type == "ctl_vehicle") %>% pull(sig_id) %>% as.character,
                 test$features %>% filter(cp_type == "ctl_vehicle") %>% pull(sig_id) %>% as.character)
  
# for each of the mechanisms  
for (mec in 1:numMechanisms) {
  
  mechanism <- names(train$targets)[mec + 1]
  
  print(paste("mechanism #: ",
              mec,
              "| mechanism title: ",
              mechanism))

  top_pcs <- paste("PC", 1:numPrincipalComponents, sep = "")
  
  x <- train$pca.features$x[,top_pcs]
  y <- train$targets[-train$heldOut,] %>% 
    select(all_of(mechanism)) %>% as.matrix() %>% as.factor
  
  fit <- train(method = method,
               x = x,
               y = y,
               trControl = trainControl(method = 'none'))
  
  prediction <- predict(fit, new.data, type = 'prob') %>% select("probability" = "1")
  
  predictions[,mechanism] <- prediction$probability %>% limitRange(xmin = mean(as.integer(as.character(y))), xmax = .98)
  predictions[ctl_vehicle,mechanism] <- 0

}


# Evaluate performance -------------------------------------------------------------------------------------------------

# evaluate performance on training samples
confusionMatrix(data = predictions[train$features$sig_id[-train$heldOut],] %>% round(0) %>% as.factor,
                reference = train$targets[-train$heldOut,] %>% select(-sig_id) %>% as.matrix %>% as.factor,
                positive = "1")
RMSE(pred = predictions[train$features$sig_id[-train$heldOut],] %>% as.matrix %>% as.integer,
     obs = train$targets[-train$heldOut,] %>% select(-sig_id) %>% as.matrix %>% as.integer)
score(prediction = predictions[train$features$sig_id[-train$heldOut],],
      reference = train$targets[-train$heldOut,] %>% select(-sig_id))


# evaluate performance on held out samples
confusionMatrix(data = predictions[train$features$sig_id[train$heldOut],] %>% round(0) %>% as.factor,
                reference = train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix %>% as.factor,
                positive = "1")
RMSE(pred = predictions[train$features$sig_id[train$heldOut],] %>% as.matrix %>% as.integer,
     obs = train$targets[train$heldOut,] %>% select(-sig_id) %>% as.matrix %>% as.integer)
score(prediction = predictions[train$features$sig_id[train$heldOut],],
      reference = train$targets[train$heldOut,] %>% select(-sig_id))

# # Investigate results --------------------------------------------------------------------------------------------------
# 
# df.prediction <-
#   predictions[train$features$sig_id[train$heldOut],] %>% 
#   as.data.frame() %>% 
#   rownames_to_column("sig_id") %>% 
#   pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction")
# 
# df.reference <-
#   train$targets[train$heldOut,] %>% 
#   pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target")
#   
# df <-
#   df.reference %>% 
#   left_join(df.prediction, by = c("sig_id", "mechanism")) %>% 
#   mutate(logLoss = - target * log(limitRange(prediction, delta, 1-delta)) - 
#            (1-target) * log(1-limitRange(prediction,delta,1-delta)))
# 
# df %>%
#   arrange(desc(logLoss))
# 
# df %>% filter(prediction == 0 & target == 1)
# 
# df %>%
#   filter(mechanism == "nfkb_inhibitor") %>% 
#   summarise(mean(target))
# 
# df %>% 
#   mutate(mean = mean(df$target)) %>% 
#   mutate(logLossMean = - target * log(mean) - (1-target)*(log(1-mean))) %>% 
#   summarise(score.prediction = mean(logLoss),
#             score.mean = mean(logLossMean))
# 
# # where is the log loss?
# df %>% 
#   mutate(stratum = cut(prediction, include.lowest = TRUE, breaks = seq(0,1,.1))) %>% 
#   group_by(stratum) %>% 
#   summarise(count = n(),
#             activation = mean(target),
#             prediction = mean(prediction),
#             meanLogLoss = mean(logLoss),
#             totalLogLoss = mean(df$logLoss) * sum(logLoss) / sum(df$logLoss))
# 
# # what mechanisms contribute the most to the log loss?
# df.by.mechanism <-
#   df %>% 
#   group_by(mechanism) %>% 
#   summarise(activations = sum(target),
#             prediction = sum(prediction),
#             logLoss = mean(logLoss), .groups = 'drop') %>% 
#   arrange(desc(logLoss))
#   
# df.by.mechanism %>%   
#   ggplot(aes(x = prediction/activations, y = logLoss)) +
#   geom_point()
# 
# df.by.mechanism %>% filter(mechanism == "proteasome_inhibitor")
# 
# df.by.mechanism %>% data.table
# 
# # mechanisms by frequency
# train$targets %>% 
#   pivot_longer(cols = - "sig_id",
#                names_to = "mechanism",
#                values_to = "value") %>% 
#   filter(value == 1) %>% 
#   select(-value) %>% 
#   group_by(mechanism) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   arrange(desc(count))
# 
# # confusion matrix for one particular mechanism
# confusionMatrix(data = df %>% filter(mechanism == "proteasome_inhibitor") %>% pull(prediction) %>% round(0) %>% as.factor,
#                 reference = df %>% filter(mechanism == "proteasome_inhibitor") %>% pull(target) %>% as.factor,
#                 positive = "1")
# 
# Calculate confidence interval of score on held out samples -----------------------------------------------------------

temp <-
  train$targets %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target")

df <- 
  predictions[train$heldOut,] %>% 
  as.data.frame() %>% 
  rownames_to_column("sig_id") %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction") %>% 
  mutate(limitedPrediction = limitRange(prediction, delta, 1-delta)) %>% 
  left_join(temp, by = c("sig_id","mechanism"))

bootstrap.scores <- replicate(10, {
  df %>% 
    sample_frac(size = .25, replace = TRUE) %>% 
    summarise(logLoss = - target * log(limitedPrediction) - (1-target) * (log(1-limitedPrediction))) %>% 
    pull(logLoss) %>% 
    mean()
})




# Save results on the test set -----------------------------------------------------------------------------------------

test$predictions <- 
  predictions[test$features$sig_id %>% as.character,] %>% 
  format(scientific = FALSE) %>% 
  as.data.frame %>% 
  rownames_to_column("sig_id")

write_csv(x = test$predictions,
          path = file.path("files","submission.csv"),
          col_names = TRUE,
          append = FALSE)


# # Plot regions for one mechanism ---------------------------------------------------------------------------------------
# 
# # parameters
# mechanism <- "nfkb_inhibitor"
# mechanism <- df.by.mechanism %>% pull(mechanism) %>% magrittr::extract(2)
# method <- "glm"
# 
# # create model for the chosen mechanism
# fit <- train(method = method,
#              x = train$pca.features$x[,1:2],
#              y = train$targets[-train$heldOut,] %>% select(matches(mechanism)) %>% as.matrix %>% as.factor,
#              trControl = trainControl(method = 'none'))
# 
# # create grid
# res <- 5
# grid <-
#   expand.grid(PC1 = seq(-res + train$pca.features$x[,1] %>% min %>% floor,
#                          res + train$pca.features$x[,1] %>% max %>% ceiling,
#                          res),
#               PC2 = seq(-res + train$pca.features$x[,2] %>% min %>% floor,
#                          res + train$pca.features$x[,2] %>% max %>% ceiling,
#                          res)) %>%
#   as.matrix()
# 
# # make predictions on the grid
# probs <- 
#   predict(fit, grid, type = "prob") %>% 
#   pull('1')
# 
# # data frame with targets for held out samples
# df <- data.frame(pc.x = train$features[train$heldOut,] %>% select(-sig_id,-cp_type,-cp_time,-cp_dose) %>% 
#                    as.matrix() %>%  magrittr::multiply_by_matrix(train$pca.features$rotation[,1]),
#                  pc.y = train$features[train$heldOut,] %>% select(-sig_id,-cp_type,-cp_time,-cp_dose) %>% 
#                    as.matrix() %>% magrittr::multiply_by_matrix(train$pca.features$rotation[,2]),
#                  target = train$targets[train$heldOut,mechanism] %>% as.matrix %>% as.factor())
# 
# # plot
# data.frame(pc.x = grid[,1],
#            pc.y = grid[,2],
#            prediction = probs) %>% 
# ggplot(aes(x = pc.x, y = pc.y,
#            fill = prediction)) +
#   geom_raster(interpolate = T, alpha = .5) +
#   geom_point(df, mapping = aes(x = pc.x, y = pc.y, color = target), inherit.aes = FALSE) +
#   scale_color_manual(values = c('blue', 'red')) +
#   scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = .5) +
#   theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
#   labs(title = paste("Predictions with method =", method, "for mechanism =", mechanism),
#        x = "PC1", y = "PC2")
# 
# 
