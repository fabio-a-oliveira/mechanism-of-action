# Instructions ---------------------------------------------------------------------------------------------------------

# If you are reading this code, here are some tips to make your life slightly easier:
# - This code is divided in sections roughly corresponding to the sections in the report.
# - CTRL + SHIFT + O shows/hides the section outline of the script, making it easier to navigate.
# - You can show/hide individual sections by clicking in the down arrow on the left.
# - CTRL + ALT + T runs an entire section.
# - CTRL + ENTER runs a line of code or the current selection and moves the cursor to the next line.
# - Some parts of this code can take more than a few minutes to run; there will be comments to alert 
#   you when appropriate.
# - The entire code takes 3-5 hours to run on a typical personal laptop.
# - Some parts of this code many require giving R access to additional memory and use of the hard drive to run; 
#   there will be comments to alert you when appropriate.

# HOUSEKEEPING ---------------------------------------------------------------------------------------------------------
# Load libraries, define some global options for the script

# logicals determining how to load and save data

runningOnKaggle <- FALSE
runningOnPC <- !runningOnKaggle

# install and load these packages

if (!require("tidyverse")) {install.packages("tidyverse"); library("tidyverse")}
if (!require("caret")) {install.packages("caret"); library("caret")}
if (!require("data.table")) {install.packages("data.table"); library("data.table")}
if (!require("matrixStats")) {install.packages("matrixStats"); library("matrixStats")}
if (!require("lubridate")) {install.packages("lubridate"); library("lubridate")}
if (!require("mda")) {install.packages("mda"); library("mda")}

# make sure these packages are installed but do not load

if (!("fs" %in% installed.packages())) install.packages("fs")
if (!("plyr" %in% installed.packages())) install.packages("plyr")
if (!("magrittr" %in% installed.packages())) install.packages("magrittr")
if (!("e1071" %in% installed.packages())) install.packages("e1071")
if (!("mda" %in% installed.packages())) install.packages("mda")
if (!("penalized" %in% installed.packages())) install.packages("penalized")

# READ FILES -----------------------------------------------------------------------------------------------------------
# Read the files required for the modeling/analysis; file names depend on whether script in running locally or
# hosted on kaggle

# read files - RUNNING ON KAGGLE NOTEBOOK

if (runningOnKaggle) {
  files <-U
    list(train_features = read_csv("../input/lish-moa/train_features.csv", col_types = cols()),
         train_targets_scored = read_csv("../input/lish-moa/train_targets_scored.csv", col_types = cols()),
         train_targets_nonscored = read_csv("../input/lish-moa/train_targets_nonscored.csv", col_types = cols()),
         test_features = read_csv("../input/lish-moa/test_features.csv", col_types = cols()))
}

# read files - RUNNING ON PC

if (runningOnPC) {
  
  files <-
    list(train_features = 
           fs::dir_ls("files") %>%
           str_subset("train_features_\\d{2}") %>%
           map_dfr(read_csv, col_types = cols(), ),
         train_targets_scored = 
           fs::dir_ls("files") %>%
           str_subset("train_targets_scored_\\d{2}") %>%
           map_dfr(read_csv, col_types = cols()),
         train_targets_nonscored = 
           fs::dir_ls("files") %>%
           str_subset("train_targets_nonscored_\\d{2}") %>%
           map_dfr(read_csv, col_types = cols()),
         test_features = 
           fs::dir_ls("files") %>%
           str_subset("test_features_\\d{2}") %>%
           map_dfr(read_csv, col_types = cols()))
}

# PREPARE TRAIN, TEST AND VALIDATE SETS --------------------------------------------------------------------------------
# Create the data partitions

# prepare train, test and predictions objects

train <-
  list(features = files$train_features,
       outcomes = files$train_targets_scored,
       outcomes_nonscored = files$train_targets_nonscored)

test <-
  list(features = files$test_features)

# create data frame to hold the predictions from each of the models

predictions <-
  expand.grid(sig_id = c(train$features$sig_id, test$features$sig_id),
              mechanism = files$train_targets_scored %>% select(-sig_id) %>% names) %>% 
  left_join(train$outcomes %>% 
              pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target"), 
            by = c("sig_id","mechanism"))

# parameters for partitioning

proportionTrainSingle <- .80
proportionTrainEnsemble <- .10
proportionValidate <- 1 - proportionTrainSingle - proportionTrainEnsemble

# create partitions

numSamplesTrainSingle <- (proportionTrainSingle * nrow(train$features)) %>% round(0)
numSamplesTrainEnsemble <- (proportionTrainEnsemble * nrow(train$features)) %>% round(0)
numSamplesValidate <- nrow(train$features) - numSamplesTrainEnsemble - numSamplesTrainSingle

# create partitions

partitions <- 
  list(controlSamples = 
         rbind(train$features,test$features) %>% 
         filter(cp_type == "ctl_vehicle") %>% 
         pull(sig_id))
partitions$trainSingle <-
  train$features %>%
  sample_n(numSamplesTrainSingle) %>% 
  pull(sig_id)
partitions$trainEnsemble <-
  train$features %>%
  filter(!(sig_id %in% partitions$trainSingle)) %>% 
  sample_n(numSamplesTrainEnsemble) %>% 
  pull(sig_id)
partitions$validate <-
  train$features %>% 
  filter(!(sig_id %in% partitions$trainSingle) & !(sig_id %in% partitions$trainEnsemble)) %>% 
  sample_n(numSamplesValidate) %>% 
  pull(sig_id)
partitions$test <-
  test$features %>% 
  pull(sig_id)

# remove temporary objects

rm(files,
   numSamplesTrainEnsemble,
   numSamplesTrainSingle,
   numSamplesValidate,
   proportionTrainEnsemble,
   proportionTrainSingle,
   proportionValidate)

# USEFUL STUFF ---------------------------------------------------------------------------------------------------------
# Define useful functions and objects

# epsilon - tolerance required for calculation of scores

eps <- 1e-15

# vector function that limits the range of an input

limitRange <- function(x, xmin, xmax){x %>%  pmin(xmax) %>% pmax(xmin)}

# function used to calculate the performance metric

logLoss <- function(obs,pred){
  pred <- limitRange(pred,eps,1-eps)
  mean(-obs*log(pred) - (1-obs)*log(1-pred))
}

# prevalence of each mechanism on the trainSingle set (not considering control samples)

prevalence <-
  train$outcomes %>% 
  left_join(train$outcomes_nonscored, by = "sig_id") %>% 
  filter(sig_id %in% partitions$trainSingle & !(sig_id %in% partitions$controlSamples)) %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  group_by(mechanism) %>% 
  summarise(prevalence = mean(target), .groups = 'drop')

# full set of mechanisms

mechanisms <- 
  data.frame(mechanism = train$outcomes %>% select(-sig_id) %>% names(),
             source = "scored") %>% 
  rbind(data.frame(mechanism = train$outcomes_nonscored %>% select(-sig_id) %>% names(),
                   source = "nonscored")) %>% 
  left_join(prevalence, by = "mechanism")

# create vector with sig_id of all samples without any mechanism

partitions$noMechanism <-
  train$outcomes %>% 
  left_join(train$outcomes_nonscored, by = "sig_id") %>% 
  pivot_longer(cols = c(-"sig_id"),
               names_to = "mechanism", values_to = "target") %>% 
  group_by(sig_id) %>% 
  summarise(totalMechanisms = sum(target), .groups = "drop") %>% 
  filter(totalMechanisms == 0) %>% 
  pull(sig_id)

# remove temporary objects

rm(prevalence)


# ENSURE TRAINING SETS CONTAIN POSITIVE ENTRIES FOR ALL MECHANISMS -----------------------------------------------------
# trainSingle needs to hold at least one sample with each mechanism present in the other sets
# trainEnsemble needs to hold at least one sample with each mechanism present in the validate set

# if needed, bring samples to trainSingle from trainEnsemble or validate

missing.mechanisms <-
  train$outcomes %>% 
  filter(sig_id %in% partitions$trainSingle) %>% 
  pivot_longer(cols = -"sig_id",
               names_to = "mechanism",
               values_to = "target") %>% 
  group_by(mechanism) %>% 
  summarise(zeros = sum(target == 0),
            ones = sum(target == 1), .groups = "drop") %>% 
  filter(ones == 0) %>% 
  pull(mechanism)

if (length(missing.mechanisms) >= 1) {
  
  samples.to.exchange <- 
    train$outcomes %>% 
    filter(sig_id %in% partitions$trainEnsemble | sig_id %in% partitions$validate) %>%
    pivot_longer(cols = -"sig_id",
                 names_to = "mechanism",
                 values_to = "target") %>% 
    filter(mechanism %in% missing.mechanisms & target == 1) %>% 
    group_by(mechanism) %>% 
    summarise(id = first(sig_id), .groups = 'drop') %>% 
    pull(id)
  
  partitions$trainSingle <- c(partitions$trainSingle, samples.to.exchange)
  partitions$trainEnsemble <- partitions$trainEnsemble[!(partitions$trainEnsemble %in% samples.to.exchange)]
  partitions$validate <- partitions$validate[!(partitions$validate %in% samples.to.exchange)]
  
}

# if needed, bring samples to trainEnsemble from validate

missing.mechanisms <-
  train$outcomes %>% 
  filter(sig_id %in% partitions$trainEnsemble) %>% 
  pivot_longer(cols = -"sig_id",
               names_to = "mechanism",
               values_to = "target") %>% 
  group_by(mechanism) %>% 
  summarise(zeros = sum(target == 0),
            ones = sum(target == 1), .groups = "drop") %>% 
  filter(ones == 0) %>% 
  pull(mechanism)

if (length(missing.mechanisms) >= 1) {
  
  samples.to.exchange <- 
    train$outcomes %>% 
    filter(sig_id %in% partitions$validate) %>%
    pivot_longer(cols = -"sig_id",
                 names_to = "mechanism",
                 values_to = "target") %>% 
    filter(mechanism %in% missing.mechanisms & target == 1) %>% 
    group_by(mechanism) %>% 
    summarise(id = first(sig_id), .groups = 'drop') %>% 
    pull(id)
  
  partitions$trainEnsemble <- c(partitions$trainEnsemble, samples.to.exchange)
  partitions$validate <- partitions$validate[!(partitions$validate %in% samples.to.exchange)]
  
}

# remove temporary variables

rm(missing.mechanisms,
   samples.to.exchange)


# FIND CENTROIDS FOR EACH MECHANISM ON trainSingle SET, RUN PCA AND GET ROTATION MATRIX --------------------------------
# Find a rotation to a coordinate system that concentrates the separation between classes in the first components

# centroids for each scored mechanism

train$centroids <- 
  train$outcomes %>% 
  filter(sig_id %in% partitions$trainSingle & !(sig_id %in% partitions$controlSamples)) %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  filter(target == 1) %>%
  select(-target) %>% 
  left_join(train$features, by = "sig_id") %>% 
  select(-cp_type,-cp_time,-cp_dose) %>%
  group_by(mechanism) %>% 
  summarise(across(-sig_id, mean), .groups = 'drop')

# centroids for samples with no scored mechanism

noMechanisms <-
  train$outcomes %>% 
  filter(sig_id %in% partitions$trainSingle & !(sig_id %in% partitions$controlSamples)) %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  group_by(sig_id) %>% 
  summarise(numMechanisms = sum(target), .groups = 'drop') %>% 
  filter(numMechanisms == 0) %>% 
  pull(sig_id)

centroidNoMechanisms <-
  train$features %>% 
  filter(sig_id %in% noMechanisms) %>% 
  select(-sig_id, -cp_type , -cp_dose, -cp_time) %>% 
  summarise(across(everything(),mean)) %>% 
  mutate(mechanism = 'none', .before = everything())

# join centroids for each mechanism and centroid for entries with no mechanism

train$centroids <- 
  train$centroids %>% 
  rbind(centroidNoMechanisms)

# perform PCA on centroids

train$pca.centroids <-
  train$centroids %>% 
  as.data.frame() %>% 
  column_to_rownames("mechanism") %>% 
  as.matrix %>% 
  prcomp()

# calculate principal components for all samples and include in features objects

train$features <- 
  train$features %>% 
  select(-sig_id, -cp_type, -cp_time, -cp_dose) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(train$pca.centroids$rotation) %>% 
  cbind(train$features) %>% 
  select(sig_id, cp_type, cp_time, cp_dose, starts_with("g-"), starts_with("c-"), starts_with("PC"))

test$features <- 
  test$features %>% 
  select(-sig_id, -cp_type, -cp_time, -cp_dose) %>% 
  as.matrix() %>% 
  magrittr::multiply_by_matrix(train$pca.centroids$rotation) %>% 
  cbind(test$features) %>% 
  select(sig_id, cp_type, cp_time, cp_dose, starts_with("g-"), starts_with("c-"), starts_with("PC"))

# remove temporary objects

rm(centroidNoMechanisms,
   noMechanisms)

# MODEL 1 - BENCHMARK --------------------------------------------------------------------------------------------------

# - Predict 0 for control samples and the mechanism prevalence (on trainSingle set) for all others

predictions$benchmark <- 
  predictions %>% 
  left_join(mechanisms, 
            by = "mechanism") %>% 
  mutate(benchmark = if_else(sig_id %in% partitions$controlSamples, 0 , prevalence)) %>% 
  pull(benchmark)

# MODEL 2 - LOGISTIC REGRESSION ----------------------------------------------------------------------------------------
# ~ 8min to run this section

# - Logistic regression applied cp_time, cp_dose and principal components
# - Prediction on control samples explicitly set to 0
# - Pre-processing: g- and c- predictors rotated to centroids hyperplane

# parameters

method = "logit"

ctrl = trainControl(method = 'none',
                    verboseIter = FALSE,
                    classProbs = TRUE,
                    summaryFunction = mnLogLoss)

numPrincipalComponents = 20

listPrincipalComponents = paste("PC", 1:numPrincipalComponents, sep = "")

# train one-versus-all models for each mechanism and calculate predictions

results <- 
  sapply(mechanisms$mechanism[mechanisms$source == "scored"],
         simplify = FALSE,
         function(selectedMechanism){
           
           # print status
           
           paste("Time: ", now(), 
                 " | Method: ",method,
                 " | Mechanism: ",selectedMechanism, 
                 " (",which(mechanisms == selectedMechanism),"/",
                 sum(mechanisms$source == "scored"),")", 
                 sep = "") %>% 
             print()
           
           # pre-process features and outcomes on trainSingle set
           
           df <- 
             rbind(train$features, test$features) %>% 
             left_join(train$outcomes, by = "sig_id") %>% 
             rename(outcome = all_of(selectedMechanism)) %>% 
             mutate(outcome = factor(outcome, levels = c("0","1"), labels = c("negative","positive"))) %>% 
             select(sig_id, cp_time, cp_dose, any_of(listPrincipalComponents), outcome)
           
           set <- (df$sig_id %in% partitions$trainSingle) & !(df$sig_id %in% partitions$controlSamples)
           
           # fit model
           
           fit <- train(outcome ~ .,
                        data = select(df, -sig_id),
                        subset = set,
                        method = "glm",
                        family = "binomial",
                        trControl = ctrl)

           # make predictions
           
           df %>% 
             mutate(mechanism = selectedMechanism,
                    prediction = predict(fit, df, type = 'prob') %>% pull(positive)) %>% 
             mutate(prediction = if_else(sig_id %in% partitions$controlSamples, 0, prediction)) %>%
             select(sig_id,mechanism,prediction)
           
         }) %>% 
  plyr::rbind.fill()

# introduce cap on minimum and maximum predictions and include in predictions object

predictions$logit <-
  predictions %>% 
  left_join(results, by = c("sig_id","mechanism")) %>%
  mutate(prediction = limitRange(prediction, benchmark/3, .99)) %>% 
  pull(prediction)

# remove temporary variables

rm(method, ctrl, numPrincipalComponents, listPrincipalComponents, results)

# MODEL 3 - K NEAREST NEIGHBORS ----------------------------------------------------------------------------------------
# ~ 1h00min minutes to run this section

# - K nearest neighbors applied to all predictors
# - Prediction on control samples explicitly set to 0
# - Pre-processing: g- and c- predictors rotated to centroids hyperplane, predictors centered and scaled

# parameters

method = "KNN"

numNeighbors = 200

ctrl = trainControl(method = 'none',
                    verboseIter = FALSE,
                    classProbs = TRUE,
                    summaryFunction = mnLogLoss)

numPrincipalComponents = 10

listPrincipalComponents = paste("PC", 1:numPrincipalComponents, sep = "")

# train one-versus-all models for each mechanism and calculate predictions

results <- 
  sapply(mechanisms$mechanism[mechanisms$source == "scored"],
         simplify = FALSE,
         function(selectedMechanism){
           
           # print status
           
           paste("Time: ", now(), 
                 " | Method: ",method,
                 " | Mechanism: ",selectedMechanism, 
                 " (",which(mechanisms == selectedMechanism),"/",
                 sum(mechanisms$source == "scored"),")", 
                 sep = "") %>% 
             print()
           
           # pre-process features and outcomes on trainSingle set
           
           df <- 
             rbind(train$features, test$features) %>% 
             select(sig_id, cp_time, cp_dose, any_of(listPrincipalComponents)) %>% 
             mutate(cp_dose = if_else(cp_dose == "D1", 1, 2)) %>% 
             left_join(train$outcomes, by = "sig_id") %>% 
             rename(outcome = all_of(selectedMechanism)) %>% 
             mutate(outcome = factor(outcome, levels = c("0","1"), labels = c("negative","positive"))) %>% 
             select(sig_id, any_of(listPrincipalComponents), outcome)
           
           set <- (df$sig_id %in% partitions$trainSingle) & !(df$sig_id %in% partitions$controlSamples)
           
           # fit model
           
           fit <- train(outcome ~ .,
                        data = select(df, -sig_id),
                        subset = set,
                        method = "knn",
                        tuneGrid = data.frame(k = numNeighbors),
                        trControl = ctrl)
           
           # make predictions
           
           df %>% 
             mutate(mechanism = selectedMechanism,
                    prediction = predict(fit, df, type = 'prob') %>% pull(positive)) %>% 
             mutate(prediction = if_else(sig_id %in% partitions$controlSamples, 0, prediction)) %>%
             select(sig_id,mechanism,prediction)
           
         }) %>% 
  plyr::rbind.fill()

# introduce cap on minimum and maximum predictions and include in predictions object

predictions$knn <-
  predictions %>% 
  left_join(results, by = c("sig_id","mechanism")) %>%
  mutate(prediction = limitRange(prediction, benchmark/3, .99)) %>% 
  pull(prediction)

# remove temporary variables

rm(method, ctrl, numPrincipalComponents, listPrincipalComponents, results)

# MODEL 4 - NAIVE BAYES WITH LOESS SMOOTHING ---------------------------------------------------------------------------
# ~ 6min minutes to run this section

# - use selection of principal components as predictors (to avoid correlation between predictors)
# - for each selected predictor, calculate smooth conditional probability of positive outcome
# - combine conditional probability obtained separately for each predictor into final prediction
# - prediction on control samples explicitly set to 0

# parameters

method = "Naive Bayes with loess smoothing"
numBins = 15000
loessSpan = .3
numPrincipalComponents = 3
listPrincipalComponents = paste("PC", 1:numPrincipalComponents, sep = "")

# train one-versus-all models for each mechanism and calculate predictions

results <- 
  sapply(mechanisms$mechanism[mechanisms$source == "scored"],
         simplify = FALSE,
         function(selectedMechanism){
           
           # print status
           
           paste("Time: ", now(), 
                 " | Method: ",method,
                 " | Mechanism: ",selectedMechanism, 
                 " (",which(mechanisms == selectedMechanism),"/",
                 sum(mechanisms$source == "scored"),")", 
                 sep = "") %>% 
             print()
           
           # # pre-process features and outcomes on trainSingle set
          
           df <-
             rbind(train$features, test$features) %>% 
             select(sig_id, any_of(listPrincipalComponents)) %>% 
             left_join(train$outcomes %>% select(sig_id, all_of(selectedMechanism)), 
                       by = "sig_id") %>% 
             rename("target" = all_of(selectedMechanism))
           
           prevalence <- mechanisms$prevalence[mechanisms$mechanism == selectedMechanism]
           
           # make prediction based on conditional probability for each PC independently
           
           pc.results <- 
             sapply(listPrincipalComponents,
                    simplify = TRUE,
                    function(PC){
                      
                      df.pc <-
                        df %>% 
                        filter(sig_id %in% partitions$trainSingle &
                                 !(sig_id %in% partitions$controlSamples)) %>% 
                        rename(PC = all_of(PC)) %>% 
                        mutate(bin = cut(PC, numBins, include.lowest = TRUE)) %>% 
                        group_by(bin) %>% 
                        summarise(PC = mean(PC), prevalence_bin = mean(target), .groups = "drop")
                      
                      fit <-
                        loess(prevalence_bin ~ PC,
                              data = df.pc,
                              family = "gaussian",
                              span = loessSpan)
                      
                      df %>% 
                        rename(PC = all_of(PC)) %>% 
                        mutate(prediction = predict(fit, PC)) %>% 
                        replace_na(list(prediction = 0)) %>% 
                        mutate(prediction = limitRange(prediction, 0, 1)) %>% 
                        mutate(prediction = 
                                 if_else(sig_id %in% partitions$controlSamples, 0, prediction)) %>% 
                        pull(prediction)
                      
                    })
           
           # assemble results for each sig_id for the relevant mechanism
           
           data.frame(sig_id = df$sig_id,
                      mechanism = selectedMechanism,
                      pc.results) %>% 
             pivot_longer(cols = any_of(listPrincipalComponents), names_to = "PC", values_to = "prediction") %>% 
             group_by(sig_id, mechanism) %>% 
             summarise(prediction = mean(prediction), .groups = "drop")
             # summarise(prediction = prod(prediction)/prevalence^(numPrincipalComponents-1), .groups = "drop")
           
         }) %>% 
  plyr::rbind.fill()

# introduce cap on minimum and maximum predictions and include in predictions object

predictions$naiveBayes <-
  predictions %>% 
  left_join(results, by = c("sig_id","mechanism")) %>%
  mutate(prediction = limitRange(prediction, benchmark/10, .99)) %>% 
  pull(prediction)

# remove temporary variables

rm(method, numBins, loessSpan, numPrincipalComponents, listPrincipalComponents, results)



# MODEL 5 - SUPPORT VECTOR CLASSIFIER ----------------------------------------------------------------------------------
# ~ 3h30min minutes to run this section

# - use cp_time, cp_dose, and selection of principal components as predictors
# - for each mechanism, tune a linear SVM algorithm to find optimal separating hyperplane and calculate class probs
# - prediction on control samples explicitly set to 0

# parameters

method <- "svmLinear2"
numPrincipalComponents <- 50
listPrincipalComponents = paste("PC", 1:numPrincipalComponents, sep = "")
grid <- data.frame(cost = c(.01, .02, .05, .10, .20, .50))
ctrl <- trainControl(method = 'LGOCV', p = .9, number = 1, 
                     verboseIter = FALSE,
                     classProbs = TRUE, summaryFunction = mnLogLoss)

# train models for each mechanism

results <- sapply(mechanisms$mechanism[mechanisms$source == "scored"],
                  simplify = FALSE,
                  function(selectedMechanism){
                    
                    # print status
                    
                    paste("Time: ", now(), 
                          " | Method: ",method,
                          " | Mechanism: ",selectedMechanism, 
                          " (",which(mechanisms == selectedMechanism),"/",
                          sum(mechanisms$source == "scored"),")", 
                          sep = "") %>% 
                      print()
                    
                    # define training set
                    
                    df <-
                      rbind(train$features, test$features) %>% 
                      left_join(train$outcomes, by = "sig_id") %>% 
                      select(sig_id, cp_type, cp_time, cp_dose,
                             any_of(listPrincipalComponents),
                             all_of(selectedMechanism)) %>% 
                      rename(outcome = all_of(selectedMechanism)) %>% 
                      mutate(outcome = factor(outcome, levels = c("0","1"), labels = c("negative","positive")),
                             cp_dose = if_else(cp_dose == "D1", 1, 2))
                    
                    # fit model
                    
                    fit <- train(outcome ~ .,
                                 df %>% 
                                   filter(sig_id %in% partitions$trainSingle &
                                            !(sig_id %in% partitions$controlSamples)) %>% 
                                   select(-sig_id, -cp_type),
                                 method = method,
                                 trControl = ctrl,
                                 tuneGrid = grid)
                    
                    # calculate predictions on all samples
                    
                    df %>% 
                      mutate(mechanism = selectedMechanism,
                             prediction = predict(fit, df, type = "prob")$positive,
                             cost = fit$bestTune$cost) %>% 
                      mutate(prediction = if_else(sig_id %in% partitions$controlSamples, 0 , prediction)) %>% 
                      select(sig_id, mechanism, prediction, cost)
                    
                  }) %>% 
  plyr::rbind.fill()

# introduce cap on minimum and maximum predictions and include in predictions object

predictions$svmLinear <-
  predictions %>% 
  left_join(results, by = c("sig_id","mechanism")) %>%
  mutate(prediction = limitRange(prediction, benchmark/10, .99)) %>% 
  pull(prediction)

# remove temporary variables

rm(method, grid, ctrl, numPrincipalComponents, listPrincipalComponents, results)

# MODEL 6 - PENALIZED MIXTURE DISCRIMINANT ANALYSIS --------------------------------------------------------------------
# ~ 30min minutes to run this section

# - use selection of principal components as predictors
# - convert multi-label problem to multi-class by combining all outcomes observed on trainSingle partition
# - prediction on control samples explicitly set to 0

# parameters

method <- "PMDA"
numPrincipalComponents <- 10
subClasses <- 2
listPrincipalComponents = paste("PC", 1:numPrincipalComponents, sep = "")

# set individual id for each combination of mechanisms found in trainSingle set

outcomeIds <-
  # train$outcomes %>% # if using only scored mechanisms
  left_join(train$outcomes, train$outcomes_nonscored, by = "sig_id") %>% # if using scored and nonscored mechanisms
  mutate(none = if_else(sig_id %in% partitions$noMechanism, 1, 0)) %>% 
  filter(sig_id %in% partitions$trainSingle &
           !(sig_id %in% partitions$controlSamples)) %>%
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  filter(target == 1) %>% 
  group_by(sig_id) %>% 
  summarise(outcome = str_c(mechanism, collapse = "|"),
            numMechanisms = sum(target),
            .groups = 'drop') %>% 
  group_by(outcome) %>% 
  summarise(count = n(),
            numMechanisms = first(numMechanisms),
            .groups = 'drop') %>% 
  arrange(desc(count)) %>% 
  mutate(outcomeId = paste("id_", str_pad(row_number(),3,side = "left",pad = "0"), sep = ""))

# match outcome ids to each sample in trainSingle set

y <- 
  # train$outcomes %>%  # if using only scored mechanisms
  left_join(train$outcomes, train$outcomes_nonscored, by = "sig_id") %>%  # if using scored and nonscored mechanisms
  mutate(none = if_else(sig_id %in% partitions$noMechanism, 1, 0)) %>% 
  filter(sig_id %in% partitions$trainSingle &
           !(sig_id %in% partitions$controlSamples)) %>%
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  filter(target == 1) %>% 
  group_by(sig_id) %>% 
  summarise(outcome = str_c(mechanism, collapse = "|"),
            .groups = 'drop') %>% 
  left_join(outcomeIds, by = "outcome") %>% 
  select(sig_id,outcomeId)

# conversion matrix between mechanisms and outcomes

conversionMatrix <-
  y %>% 
  group_by(outcomeId) %>% 
  sample_n(1) %>% 
  # left_join(train$outcomes, by = "sig_id") %>% # if using only scored mechanisms
  left_join(left_join(train$outcomes, train$outcomes_nonscored, by = "sig_id"), by = "sig_id") %>% # if using scored and nonscored mechanisms
  select(-sig_id) %>% 
  column_to_rownames("outcomeId") %>% 
  as.matrix()

# training set

x <-
  y %>% 
  left_join(train$features, by = "sig_id") %>% 
  select(any_of(listPrincipalComponents))

df <- cbind(y,x) %>% select(-sig_id) 

# PMDA model

fit <- mda(outcomeId ~ .,
           df,
           subClasses = subClasses,
           method = gen.ridge)

# make predictions

x.full <-
  rbind(train$features, test$features) %>% 
  select(any_of(listPrincipalComponents))

results <- 
  predict(fit, x.full, type = "posterior")  %>% 
  magrittr::multiply_by_matrix(conversionMatrix) %>% 
  as.data.frame() %>% 
  mutate(sig_id = rbind(train$features, test$features)$sig_id) %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "prediction") %>% 
  right_join(predictions, by = c("mechanism","sig_id")) %>% 
  mutate(prediction = if_else(sig_id %in% partitions$controlSamples, 0, prediction)) %>% 
  select(sig_id,mechanism,prediction)

# introduce cap on minimum and maximum predictions and include in predictions object

predictions$pmda <-
  predictions %>% 
  left_join(results, by = c("sig_id","mechanism")) %>%
  mutate(prediction = limitRange(prediction, benchmark/10, .99)) %>% 
  pull(prediction)

# remove temporary variables

rm(method, numPrincipalComponents, listPrincipalComponents, results, x, y, outcomeIds, x.full, conversionMatrix)

# MODEL 7 - ENSEMBLE (MEAN) --------------------------------------------------------------------------------------------

# make predictions as the simple average of the predictions from the previous methods

predictions$ensembleMean <-
  predictions %>% 
  mutate(ensembleMean = (benchmark + knn + logit + naiveBayes + pmda + svmLinear) / 6) %>% 
  pull(ensembleMean)

# MODEL 8 - ENSEMBLE (LINEAR REGRESSION) -------------------------------------------------------------------------------

# - calculate residuals from prediction made by simple ensemble
# - fit linear regression model using predictions from other models as predictors
# - apply L1 and L2 regularization to avoid large coefficients due to highly correlated predictors

# parameters

method <- "penalized"
grid <- expand.grid(lambda1 = .1, lambda2 = c(.1, .2, .5, 1, 2, 5, 10))
ctrl <- trainControl(method = 'cv', number = 10, verboseIter = FALSE)

# train models for each mechanism

results <- 
  sapply(mechanisms$mechanism[mechanisms$source == "scored"],
         simplify = FALSE,
         function(selectedMechanism){
           
           # print status
           
           paste("Time: ", now(), 
                 " | Method: ",method,
                 " | Mechanism: ",selectedMechanism, 
                 " (",which(mechanisms == selectedMechanism),"/",
                 sum(mechanisms$source == "scored"),")", 
                 sep = "") %>% 
             print()
           
           # define training set
           
           df <-
             predictions %>% 
             filter(mechanism == selectedMechanism) %>% 
             mutate(residual = target - ensembleMean) %>% 
             select(-target)
           
           # fit model
           
           fit <- train(residual ~ .,
                        data = 
                          df %>%
                          filter(sig_id %in% partitions$trainEnsemble &
                                   !(sig_id %in% partitions$controlSamples)) %>% 
                          select(-sig_id, -mechanism, -ensembleMean),
                        method = method,
                        tuneGrid = grid,
                        trControl = ctrl,
                        trace = FALSE)
           
           # calculate predictions on all samples
           
           df %>% 
             mutate(prediction = predict(fit, df) + df$ensembleMean) %>% 
             mutate(prediction = limitRange(prediction, 0, 1)) %>% 
             mutate(prediction = if_else(sig_id %in% partitions$controlSamples, 0, prediction),
                    lambda1 = fit$bestTune$lambda1,
                    lambda2 = fit$bestTune$lambda2) %>% 
             select(sig_id, mechanism, prediction, lambda1, lambda2)
           
         }) %>% 
  plyr::rbind.fill()

# introduce cap on minimum and maximum predictions and include in predictions object

predictions$ensembleWeighted <-
  predictions %>% 
  left_join(results, by = c("sig_id","mechanism")) %>%
  mutate(prediction = limitRange(prediction, benchmark/10, .99)) %>% 
  pull(prediction)

# remove temporary variables

rm(method, grid, ctrl, results)

# MEASURE PERFORMANCE --------------------------------------------------------------------------------------------------

predictions %>% 
  filter(sig_id %in% partitions$validate) %>% 
  pivot_longer(cols = -c("sig_id","mechanism","target"),
               values_to = "prediction", names_to = "method") %>% 
  group_by(method) %>% 
  summarise(logLoss = logLoss(target,prediction),
            RMSE = RMSE(target,prediction),
            accuracy = mean(round(prediction,0) == target),
            balancedAccuracy = .5*(1-mean(round(prediction[target == 0],0)) + 
                                     .5*mean(round(prediction[target == 1],0))),
            .groups = "drop")

# SAVE RESULTS ---------------------------------------------------------------------------------------------------------
# if running on a kaggle notebook, save predictions on test set on format required for submission
# if running on PC, save all files required for the RMarkdown report

# save results - RUNNING ON KAGGLE NOTEBOOK

if (runningOnKaggle) {
  predictions %>% 
    mutate(prediction = ensembleWeighted) %>% 
    select(sig_id, mechanism, prediction) %>% 
    pivot_wider(names_from = "mechanism", values_from = "prediction") %>% 
    write_csv(path = "./submission.csv",
              col_names = TRUE,
              append = FALSE)
}

# save results - RUNNING ON PC

if (runningOnPC) {
  
  # save mechanisms and partitions objects
  
  write_csv(x = mechanisms,
            file = file.path("files","mechanisms.csv"),
            append = FALSE)
  
  save(list = c("partitions"),
       file = file.path("files", "partitions.RData"))
  
  # partition the predictions object and save each chunk in a separate file
  
  numPartitions = 200
  
  for (i in 1:numPartitions) {
    
    fileName = paste("predictions_",
                     str_pad(i, side = "left", width = 3, pad = "0"),
                     ".csv", sep = "")
    
    predictions %>% 
      mutate(fileNum = sort(rep(x = 1:numPartitions, 
                                length.out = nrow(predictions)))) %>% 
      filter(fileNum == i) %>% 
      select(-fileNum) %>% 
      write_csv(file = file.path("files",fileName),
                append = FALSE)
  }
  
  rm(numPartitions)
}


