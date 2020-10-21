# Model007.R

# using KNN to explore dataset

# Housekeeping ---------------------------------------------------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse") else library("tidyverse")
if (!require("caret")) install.packages("caret") else library("caret")
if (!require("rpart")) install.packages("rpart") else library("rpart")
# if (!require("MASS")) install.packages("MASS") else library("MASS")
if (!require("data.table")) install.packages("data.table") else library("data.table")
if (!require("kknn")) install.packages("kknn") else library("kknn")


options("datatable.print.topn" = 20)

# Declaration of some useful functions ---------------------------------------------------------------------------------

limitRange <- function(x){
  delta <- 10^(-15)
  pmin(pmax(x,delta),1-delta)}

score <- 
  function(reference,prediction){
    N <- nrow(reference)
    M <- ncol(reference)
    - sum(reference*log(limitRange(prediction)) + (1-reference)*log(1-limitRange(prediction))) / (N*M)
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
  filter(cp_time == "48", cp_type == "trt_cp", cp_dose == "D1") %>%
  select(- c(cp_time, cp_type, cp_dose))

train$targets <-
  train$targets %>% 
  filter(sig_id %in% train$features$sig_id)

# Run Principal Component Analysis on simplified set -------------------------------------------------------------------

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

# # Correlation between predictors and outcome ---------------------------------------------------------------------------
# 
# numDimensions <- 2
# 
# x <- train$features[-train$heldOut,] %>% select(-sig_id) %>% as.matrix()
# pc <- train$pca.features$x[-train$heldOut,] %>% as.matrix()
# y <- train$targets$nfkb_inhibitor[-train$heldOut] %>% as.matrix()
# 
# features_correlations <-
#   data.frame(feature = colnames(x),
#              correlation = cor(x,y)) %>% 
#   arrange(desc(correlation))
# 
# pc.correlations <-
#   data.frame(PC = colnames(pc),
#              correlation = cor(pc, y)) %>% 
#   arrange(desc(abs(correlation)))
# 
# # show most extreme correlations between PCs and outcomes
# pc.correlations %>% 
#   head(20)
# 
# highest.correlations <-
#   pc.correlations %>% 
#   head(numDimensions) %>% 
#   .$PC
# 
# # Explore nearest neighbors --------------------------------------------------------------------------------------------
# 
# # configure targets and predictors - configured for regression
# predictors <- train$pca.features$x[-train$heldOut,highest.correlations,drop=FALSE] # PCs with highest correlation with outcome
# targets <- train$targets$nfkb_inhibitor[-train$heldOut] # outcomes without held out samples
# 
# # configs for the model training
# # grid <- data.frame(k = c(1,3,5,10,20,50,100,200,400,800,1600,2005)) # very narrow to very wide neighborhood
# grid <- data.frame(k = c(1,3,5,15,25,75,175,255,425,925))
# # grid <- data.frame(k = c(1,3,5,10,20,50,100,200,400,2000))
# # grid <- data.frame(k = c(1,3,5,10,20,50,100,200,400)) # fit crashed at 800....
# 
# ctrl <- trainControl(method = 'cv', number = 5, verboseIter = TRUE, allowParallel = FALSE)
# 
# # train model
# fit_KNN <-
#   train(method = "knn",
#         x = predictors, 
#         y = targets, 
#         use.all = FALSE,
#         tuneGrid = grid,
#         trControl = ctrl)
# 
# 
# # calculate neighborhoods
# 
# N <- nrow(predictors)
# 
# distances <- sqrt(
#   (predictors[,1,drop=F] %>% matrix(nrow = N, ncol = N, byrow = F) - predictors[,1,drop=F] %>% matrix(nrow = N, ncol = N, byrow = T)) ^ 2 +
#     (predictors[,2,drop=F] %>% matrix(nrow = N, ncol = N, byrow = F) - predictors[,2,drop=F] %>% matrix(nrow = N, ncol = N, byrow = T)) ^ 2)
# 
# nearest.neighbors <-
#   distances %>% 
#   as.data.frame() %>% 
#   mutate(observation1 = row_number()) %>% 
#   pivot_longer(cols = -"observation1", names_to = "observation2", values_to = "distance") %>% 
#   mutate(observation2 = str_match(observation2,"\\d+") %>% as.integer) %>% 
#   filter(observation1 != observation2) %>% 
#   group_by(observation1) %>% 
#   # top_n(wt = -distance,n = 100) %>% 
#   mutate(target1 = targets[observation1],
#          target2 = targets[observation2]) %>% 
#   select(observation1,target1,observation2,target2,distance)
#   
#   
# nearest.neighbors %>% 
#   arrange(distance)
# 
# nearest.neighbors %>% 
#   ggplot(aes(x = distance, color = as.factor(target1-target2))) +
#   geom_density()
# 
# nearest.neighbors %>% 
#   mutate(distance = cut(distance,include.lowest = T,ordered_result = T,breaks = 0:170)) %>% 
#   group_by(distance) %>% 
#   summarise(rmse = RMSE(target1,target2)) %>% 
#   ggplot(aes(x = distance, y = rmse)) +
#   geom_point()
#   
# 
# 
# nearest.neighbors %>% 
#   # mutate(distance = cut(distance,include.lowest = T,ordered_result = T,breaks = 0:170)) %>% 
#   arrange(desc(distance)) %>% 
#   ggplot(aes(x = distance)) +
#   geom_density()
# 
# # number of points within a distance up to 5 (low RMSE)
# nearest.neighbors %>% 
#   filter(distance < 5) %>% 
#   group_by(observation1) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   arrange(count) %>% 
#   data.table()
# 
# # error in KNN as a function of number of points in the neighborhood
# nearest.neighbors %>% 
#   filter(distance < 5) %>% 
#   group_by(observation1) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   right_join(nearest.neighbors %>% filter(distance < 5), by = "observation1") %>% 
#   group_by(observation1) %>% 
#   summarise(target = first(target1),
#             neighborhood = mean(target2),
#             distance = mean(distance),
#             count = first((count))) %>% 
#   ggplot(aes(x = count, y = target - neighborhood)) +
#   geom_point() +
#   scale_x_continuous(trans = 'log10')
# 
# # prediction as average of points inside radius of 5
# nearest.neighbors %>% 
#   filter(distance < 5) %>% 
#   group_by(observation1) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   right_join(nearest.neighbors %>% filter(distance < 5), by = "observation1") %>% 
#   group_by(observation1) %>% 
#   summarise(target = first(target1),
#             prediction = mean(target2),
#             .groups = 'drop') %>% 
#   summarise(rmse.prediction = RMSE(target,prediction),
#             score.prediction = score(target %>% matrix(ncol=1),
#                                      prediction %>% matrix(ncol=1)),
#             rmse.mean = RMSE(target,mean(target)),
#             score.mean = score(target %>% matrix(ncol=1),
#                                mean(target)))
#   
# # investigate outliers from previous graph
# outliers <- 
#   nearest.neighbors %>% 
#   filter(distance < 5) %>% 
#   group_by(observation1) %>% 
#   summarise(count = n(), .groups = 'drop') %>% 
#   right_join(nearest.neighbors %>% filter(distance < 5), by = "observation1") %>% 
#   group_by(observation1) %>% 
#   summarise(target = first(target1),
#             neighborhood = mean(target2),
#             distance = mean(distance),
#             count = first((count))) %>% 
#   filter(target-neighborhood > .75) %>% 
#   arrange(distance) %>% 
#   pull(observation1)
#   
# nearest.neighbors %>% 
#   filter(observation1 %in% outliers) %>% 
#   group_by(observation1) %>% 
#   summarise(target = first(target1),
#             n_0 = sum(target2 == 0),
#             n_1 = sum(target2 == 1))
# 
# df <- data.frame(observation = 1:nrow(distances),
#                  target = targets %>% as.factor(),
#                  pc1 = predictors[,1],
#                  pc2 = predictors[,2])
# 
# df %>% 
# ggplot(aes(x = pc1, y = pc2, color = target)) +
#   geom_point(alpha = .5) +
#   geom_point(data = df %>% filter(observation %in% outliers),
#              aes(x = pc1, y = pc2),
#              alpha = 1, color = "black")
# 
# 
# 
# # LM -------------------------------------------------------------------------------------------------------------------
# 
# # train model
# fit_LM <-
#   train(method = "lm",
#         x = predictors, 
#         y = targets, 
#         trControl = ctrl)
# 
# # QDA ------------------------------------------------------------------------------------------------------------------
# 
# # train model
# fit_QDA <-
#   train(method = "qda",
#         x = predictors, 
#         y = targets %>% as.factor, 
#         trControl = ctrl)
# 
# # proof of concept - plot regions - most common mechanism --------------------------------------------------------------
# 
# # which PCs to plot
# pc.x = "PC1"
# pc.y = "PC2"
# 
# # supporting df
# df <- data.frame(observation = 1:nrow(train$features[-train$heldOut,]),
#                  sig_id = train$features$sig_id[-train$heldOut],
#                  target = train$targets$nfkb_inhibitor[-train$heldOut],
#                  train$pca.features$x[-train$heldOut,])
# 
# # create grid
# res <- 5
# grid <- 
#   expand.grid(pc.x = seq(df[,pc.x] %>% min %>% floor,
#                          df[,pc.x] %>% max %>% ceiling,
#                          res),
#               pc.y = seq(df[,pc.y] %>% min %>% floor,
#                          df[,pc.y] %>% max %>% ceiling,
#                          res)) %>% 
#   as.matrix()
# 
# # predictors
# predictors <- train$pca.features$x[-train$heldOut,c(pc.x,pc.y)]
# 
# # calculate distances
# N.grid <- nrow(grid)
# N.predictors <- nrow(predictors)
# 
# distances <-
#   sqrt((matrix(grid[,1], nrow = N.grid, ncol = N.predictors, byrow = F) -
#           matrix(predictors[,1], nrow = N.grid, ncol = N.predictors, byrow = T))^2 +
#          (matrix(grid[,2], nrow = N.grid, ncol = N.predictors, byrow = F) -
#             matrix(predictors[,2], nrow = N.grid, ncol = N.predictors, byrow = T))^2)
# 
# colnames(distances) <- train$targets$sig_id[-train$heldOut]
# 
# # # average of points in vicinity
# # vicinity <- 5
# # cbind(grid,
# #       distances %>% as.data.frame) %>% 
# #   pivot_longer(cols = -c(pc.x,pc.y), names_to = "sig_id", values_to = "distance") %>% 
# #   filter(distance < vicinity) %>% 
# #   left_join(select(train$targets,sig_id,nfkb_inhibitor), by = "sig_id") %>% 
# #   group_by(pc.x,pc.y) %>% 
# #   summarise(count = n(),
# #             mean = mean(nfkb_inhibitor),
# #             .groups = "drop") %>% 
# #   ggplot(aes(x = pc.x, y = pc.y, 
# #              color = round(mean,0) %>% as.factor)) +
# #   geom_point()
# 
# # average of k nearest neighbors
# k <- 10
# cbind(grid,
#       distances %>% as.data.frame) %>% 
#   pivot_longer(cols = -c(pc.x,pc.y), names_to = "sig_id", values_to = "distance") %>% 
#   left_join(select(train$targets,sig_id,nfkb_inhibitor), by = "sig_id") %>% 
#   group_by(pc.x,pc.y) %>% 
#   slice_min(n = k, order_by = distance) %>% 
#   summarise(prediction = mean(nfkb_inhibitor),
#             mean.distance = mean(distance), .groups = 'drop') %>% 
#   ggplot(aes(x = pc.x, y = pc.y,
#              fill = prediction)) +
#   geom_raster(interpolate = T, alpha = .5) +
#   geom_point(data = df %>% select(pc.x = all_of(pc.x),
#                                   pc.y = all_of(pc.y),target),
#              mapping = aes(x = pc.x, y = pc.y, color = target %>% as.factor)) +
#   scale_color_manual(values = c('blue', 'red')) +
#   scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = .5) +
#   theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# 
# # plot regions for any mechanism ---------------------------------------------------------------------------------------
# 
# # mechanisms by frequency
# mechanisms <- 
#   train$targets %>% 
#   select(-sig_id) %>% 
#   summarise_all(sum) %>% 
#   pivot_longer(cols = everything(), names_to = "mechanism", values_to = "count") %>% 
#   arrange(desc(count))
# 
# # mechanism under analysis
# mechanism <- mechanisms$mechanism[1]
# mechanism
# 
# # PCs with stronger correlation outside held out samples
# correlations <-
#   cor(train$targets[-train$heldOut,] %>% select(all_of(mechanism)) %>% as.matrix(),
#       train$pca.features$x[-train$heldOut,] %>% as.matrix()) %>% 
#   as.data.frame() %>%
#   pivot_longer(cols = everything(), names_to = "PC", values_to = "correlation") %>% 
#   arrange(desc(abs(correlation)))
# 
# # chose top 2 PCs with stronger correlation
# top.correlations <- correlations %>% head(2) %>% pull(PC)
# pc.x <- top.correlations[1]
# pc.y <- top.correlations[2]
# 
# # supporting df with point to plot - not held out
# df <- data.frame(observation = 1:nrow(train$features[-train$heldOut,]),
#                  sig_id = train$features$sig_id[-train$heldOut],
#                  target = train$targets$nfkb_inhibitor[-train$heldOut],
#                  train$pca.features$x[-train$heldOut,])
# 
# # supporting df with point to plot - held out
# df <- data.frame(observation = 1:nrow(train$features[train$heldOut,]),
#                  sig_id = train$features$sig_id[train$heldOut],
#                  target = train$targets$nfkb_inhibitor[train$heldOut],
#                  train$pca.features$x[train$heldOut,])
# 
# # create grid
# res <- 5
# grid <- 
#   expand.grid(pc.x = seq(-res + df[,pc.x] %>% min %>% floor,
#                          res + df[,pc.x] %>% max %>% ceiling,
#                          res),
#               pc.y = seq(-res + df[,pc.y] %>% min %>% floor,
#                          res + df[,pc.y] %>% max %>% ceiling,
#                          res)) %>% 
#   as.matrix()
# 
# # predictors
# predictors <- train$pca.features$x[-train$heldOut,c(pc.x,pc.y)]
# 
# # calculate distances
# N.grid <- nrow(grid)
# N.predictors <- nrow(predictors)
# 
# distances <-
#   sqrt((matrix(grid[,1], nrow = N.grid, ncol = N.predictors, byrow = F) -
#           matrix(predictors[,1], nrow = N.grid, ncol = N.predictors, byrow = T))^2 +
#          (matrix(grid[,2], nrow = N.grid, ncol = N.predictors, byrow = F) -
#             matrix(predictors[,2], nrow = N.grid, ncol = N.predictors, byrow = T))^2)
# 
# colnames(distances) <- train$targets$sig_id[-train$heldOut]
# 
# # average of k nearest neighbors
# k <- 10
# cbind(grid,
#       distances %>% as.data.frame) %>% 
#   pivot_longer(cols = -c(pc.x,pc.y), names_to = "sig_id", values_to = "distance") %>% 
#   left_join(select(train$targets,sig_id,nfkb_inhibitor), by = "sig_id") %>% 
#   group_by(pc.x,pc.y) %>% 
#   slice_min(n = k, order_by = distance) %>% 
#   summarise(prediction = mean(nfkb_inhibitor),
#             mean.distance = mean(distance), .groups = 'drop') %>% 
#   ggplot(aes(x = pc.x, y = pc.y,
#              fill = prediction)) +
#   geom_raster(interpolate = T, alpha = .5) +
#   geom_point(data = df %>% select(pc.x = all_of(pc.x),
#                                   pc.y = all_of(pc.y),target),
#              mapping = aes(x = pc.x, y = pc.y, color = target %>% as.factor)) +
#   scale_color_manual(values = c('blue', 'red')) +
#   scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = .5) +
#   theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# 
# # create predictions for each mechanism --------------------------------------------------------------------------------
# 
# numMechanisms <- (ncol(train$targets)-1)
# # numMechanisms <- 1
# 
# numNeighbors <- 10
# numNeighbors <- 300
# 
# predictions <-
#   sapply(1:numMechanisms, simplify = TRUE,
#          function(mec){
#            y <- train$targets %>% 
#              select(-sig_id) %>% 
#              magrittr::extract(-train$heldOut,mec) %>% 
#              matrix(ncol = 1)
#            colnames(y) <- "target"
# 
#            top.correlations <-
#              cor(y,
#                  train$pca.features$x[-train$heldOut,] %>% as.matrix()) %>% 
#              as.data.frame() %>%
#              pivot_longer(cols = everything(), names_to = "PC", values_to = "correlation") %>% 
#              arrange(desc(abs(correlation))) %>% 
#              head(2) %>% 
#              pull(PC)
#            
#            pc.x <- top.correlations[1]
#            pc.y <- top.correlations[2]
#            
#            x <- train$pca.features$x[-train$heldOut,c(pc.x,pc.y)]
#            
#            # fit <- knnreg(x = x, y = y, k = numNeighbors)
#            
#            fit <- train(method = "kknn",
#                         x = x,
#                         y = y %>% as.numeric(),
#                         metric = "RMSE",
#                         tuneGrid = data.frame(kmax = numNeighbors,
#                                               distance = 100,
#                                               kernel = 'gaussian'),
#                         trControl = trainControl(method = 'none'))
#            
#            prediction <- predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y)])
#          })
# 
# targets <-
#   train$targets[train$heldOut,] %>% 
#   select(-sig_id) %>% 
#   as.matrix
# 
# # overall results
# RMSE(predictions, targets)
# score(targets,predictions)
# 
# confusionMatrix(predictions %>% round(0) %>% as.factor(),
#                 targets %>% as.factor(),
#                 positive = "1")
# 
# # results for nfkb_inhibitor
# nfkb <- which(colnames(targets) == "nfkb_inhibitor")
# RMSE(targets[,nfkb], predictions[,nfkb])
# score(targets[,nfkb, drop = F], predictions[,nfkb, drop = F])
# 
# targets[,nfkb] %>% mean
# predictions[,nfkb] %>% mean
# 
# nfkb <-
#   data.frame(target = targets[,nfkb],
#              prediction = predictions[,nfkb]) %>%
#   mutate(score = (- target*log(limitRange(prediction)) - (1-target)*log(1-limitRange(prediction)))/(n()*206),
#          prediction.round = round(prediction,0),
#          score.round = (- target*log(limitRange(prediction.round)) - (1-target)*log(1-limitRange(prediction.round)))/(n()*206))
# 
# nfkb %>%
#   arrange(desc(abs(score))) %>%
#   head(50)
# 
# nfkb$score %>% sum
# nfkb$score.round %>% sum()
# nfkb %>% str

# # same as above, nfkb only ---------------------------------------------------------------------------------------------
# 
# # numNeighbors <- 10
# # numNeighbors <- 2000
# 
# y <- train$targets$nfkb_inhibitor[-train$heldOut] %>%
#   matrix(ncol = 1)
# colnames(y) <- "target"
# 
# top.correlations <-
#   cor(y,
#       train$pca.features$x[-train$heldOut,] %>% as.matrix()) %>%
#   as.data.frame() %>%
#   pivot_longer(cols = everything(), names_to = "PC", values_to = "correlation") %>%
#   arrange(desc(abs(correlation))) %>%
#   head(2) %>%
#   pull(PC)
# 
# pc.x <- top.correlations[1]
# pc.y <- top.correlations[2]
# 
# x <- train$pca.features$x[-train$heldOut,c(pc.x,pc.y)]
# 
# fit <- train(method = "kknn",
#              x = x,
#              y = y %>% as.numeric(),
#              # metric = "RMSE",
#              # metric = "logLoss",
#              # tuneGrid = data.frame(kmax = numNeighbors,
#              #                       distance = 100,
#              #                       kernel = 'gaussian'),
#              # tuneGrid = expand.grid(kmax = c(10,50,100,300,600),
#              #                        distance = c(1,10,100),
#              tuneGrid = expand.grid(kmax = 300,
#                                     distance = 100,
#                                     kernel = 'gaussian'),
#              trControl = trainControl(method = 'none',verboseIter = TRUE, allowParallel = FALSE))
# 
# # fit <- train(method = 'qda',
# #              x = x,
# #              y = y,
# #              trControl = trainControl(method = 'none', verboseIter = T, allowParallel = F))
# 
# 
# predictions <- predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y)])
# 
# targets <-
#   train$targets$nfkb_inhibitor[train$heldOut] %>%
#   as.matrix
# 
# # overall results
# RMSE(targets,predictions)
# score(targets,predictions)
# 
# confusionMatrix(predictions %>% round(0) %>% as.factor(),
#                 targets %>% as.factor(),
#                 positive = "1")
# 
# # comparison: results as the average
# avg <- train$targets$nfkb_inhibitor[-train$heldOut] %>% mean
# RMSE(targets, avg)
# score(targets,avg)
# 
# confusionMatrix(avg %>% matrix(ncol = 1, nrow = nrow(targets)) %>% round(0) %>% as.factor(),
#                 targets %>% as.factor(),
#                 positive = "1")
# 
# # investigate results
# 
# results <-
#   data.frame(target = train$targets$nfkb_inhibitor[train$heldOut],
#              prediction = predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y)])) %>%
#   mutate(prediction.score = - (target*log(limitRange(prediction)) + (1-target)*log(1-limitRange(prediction))),
#          average = mean(train$targets$nfkb_inhibitor[-train$heldOut]),
#          average.score = - (target*log(limitRange(average)) + (1-target)*log(1-limitRange(average))))
# 
# mean(results$prediction.score)
# mean(results$average.score)
# 
# # results %>% head(30)
# 
# # results %>%
# #   filter(prediction.score > average.score)
# 
# results %>%
#   mutate(prediction.bin = cut(prediction, include.lowest = T, breaks = seq(0,1,.1))) %>%
#   ggplot(aes(x = prediction.bin, y = target)) +
#   geom_boxplot()
# 
# # where are the most costly mistakes?
# results.bins <-
#   results %>%
#   mutate(prediction.bin = cut(prediction, include.lowest = T, breaks = seq(0,1,.1))) %>%
#   group_by(prediction.bin) %>%
#   summarise(count = n(),
#             mean.prediction = mean(prediction),
#             mean.score = mean(prediction.score),
#             total.score = sum(prediction.score)/nrow(results),
#             mean.target = mean(target),
#             .groups = 'drop')
# 
# results.bins %>%
#   pivot_longer(cols = -c('prediction.bin','mean.prediction'), values_to = 'value', names_to = 'stat') %>%
#   ggplot(aes(x = prediction.bin, y = value)) +
#   facet_wrap(. ~ stat, scales = 'free') +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 45))
# 
# results.bins
# 
# 
# # tweeking the predictions
# results %>%
#   # mutate(prediction = if_else(prediction < .1, mean(train$targets$nfkb_inhibitor[-train$heldOut]), prediction)) %>%
#   mutate(prediction = if_else(prediction < .1, .00525, prediction)) %>%
#   mutate(prediction.score = - (target*log(limitRange(prediction)) + (1-target)*log(1-limitRange(prediction))),
#          average = mean(train$targets$nfkb_inhibitor[-train$heldOut]),
#          average.score = - (target*log(limitRange(average)) + (1-target)*log(1-limitRange(average)))) %>%
#   summarise(mean(prediction.score))
# 
# mean(results$prediction.score)

# same as above, 3 dimensions, other mechanisms ------------------------------------------------------------------------

mechanism <- "serotonin_receptor_antagonist"

rm(fit)

y <- 
  train$targets[-train$heldOut,] %>% 
  select(matches(mechanism)) %>% 
  as.matrix
colnames(y) <- "target"

top.correlations <-
  cor(y,
      train$pca.features$x[-train$heldOut,] %>% as.matrix()) %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "PC", values_to = "correlation") %>% 
  arrange(desc(abs(correlation))) %>%
  head(3) %>% 
  pull(PC)

pc.x <- top.correlations[1]
pc.y <- top.correlations[2]
pc.z <- top.correlations[3]
# pc.w <- top.correlations[4]

x <- train$pca.features$x[-train$heldOut,c(pc.x,pc.y)]
# x <- train$pca.features$x[-train$heldOut,c(pc.x,pc.y,pc.z)]


# KKNN - more flexible
# fit <- train(method = "kknn",
#              x = x,
#              y = y %>% as.numeric(),
#              tuneGrid = expand.grid(kmax = 300, distance = 100, kernel = 'gaussian'),
#              trControl = trainControl(method = 'none',verboseIter = TRUE, allowParallel = FALSE))
# predictions <- predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y)], type = 'prob')
# predictions <- predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y,pc.z)], type = 'prob')

# QDA - a lot faster
fit <- train(method = 'qda',
             x = x,
             y = y %>% as.factor(),
             trControl = trainControl(method = 'none', verboseIter = T, allowParallel = F))

predictions <- predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y)], type = 'prob')
# predictions <- predict(fit, train$pca.features$x[train$heldOut,c(pc.x,pc.y,pc.z)], type = 'prob')
predictions <- predictions[,2] # for classification model (qda)

targets <- 
  train$targets[train$heldOut,] %>% 
  select(matches(mechanism)) %>% 
  as.matrix

# overall results
confusionMatrix(predictions %>% round(0) %>% as.factor(),
                targets %>% as.factor(),
                positive = "1")

RMSE(targets,predictions)
score(targets,predictions)

# comparison: results as the average
avg <- train$targets[-train$heldOut,] %>% select(matches(mechanism)) %>% as.matrix() %>% mean
RMSE(targets, avg)
score(targets,avg)

# tweaking the predictions
results %>%
  mutate(prediction = pmin(pmax(prediction, avg), .98)) %>% 
  mutate(prediction.score = - (target*log(limitRange(prediction)) + (1-target)*log(1-limitRange(prediction))),
         average = mean(train$targets$nfkb_inhibitor[-train$heldOut]),
         average.score = - (target*log(limitRange(average)) + (1-target)*log(1-limitRange(average)))) %>%
  summarise(mean(prediction.score))


# # unrelated - correlations between targets and predictors (instead of PCs) ---------------------------------------------
# 
# 
# correlations <-
#   cor(train$targets[-train$heldOut,] %>% select(-sig_id) %>% as.matrix(),
#       train$features[-train$heldOut,] %>% select(-sig_id) %>% as.matrix()) %>% 
#   as.data.frame() %>%
#   rownames_to_column("mechanism") %>% 
#   as.tibble %>% 
#   pivot_longer(cols = -"mechanism", names_to = "feature", values_to = "correlation") %>% 
#   arrange(desc(abs(correlation)))
# 
# highest.correlations <- 
#   correlations %>% 
#   group_by(mechanism) %>% 
#   summarise(correlation = max(abs(correlation)) * first(sign(correlation)))
#   ◘
# highest.correlations %>% 
#   ggplot(aes(x = correlation)) +
#   geom_density()
# 
# highest.correlations %>% 
#   arrange(desc(abs(correlation)))
# how many mechanisms per sample?



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