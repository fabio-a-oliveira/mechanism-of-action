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
         mutate(sig_id = as.factor(sig_id)))

train$outcomes <- 
  files$train_targets_scored %>% 
  filter(sig_id %in% train$features$sig_id) %>% 
  mutate(sig_id = as.factor(sig_id))
train$outcomes_nonscored <-
  files$train_targets_nonscored %>% 
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
validate$outcomes_nonscored <- 
  files$train_targets_nonscored %>% 
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


# Explore nonscored targets --------------------------------------------------------------------------

# scored and nonscored combined
targets <-
  left_join(train$outcomes, train$outcomes_nonscored, by = "sig_id")
  
# mechanisms for the combined table
targets %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  group_by(sig_id) %>% 
  summarise(mechanisms = sum(target),
            .groups = 'drop') %>% 
  pull(mechanisms) %>% 
  table()

# mechanisms for the scored targets only
train$outcomes %>% 
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  group_by(sig_id) %>% 
  summarise(mechanisms = sum(target),
            .groups = 'drop') %>% 
  pull(mechanisms) %>% 
  table()

# frequency of scored and nonscored targets
frequency <-
  targets %>%
  pivot_longer(cols = -"sig_id", names_to = "mechanism", values_to = "target") %>% 
  group_by(mechanism) %>% 
  summarise(count = sum(target), .groups = 'drop') %>% 
  mutate(scored = if_else(mechanism %in% names(train$outcomes),
                          "scored", "nonscored")) %>% 
  arrange(desc(count))

frequency %>% 
  ggplot(aes(x = count, color = scored)) +
  geom_density()

frequency %>% 
  group_by(scored) %>% 
  summarise(min = min(count),
            mean = mean(count),
            max = max(count),
            .groups = 'drop')

frequency %>% filter(count == 0) %>% as.data.frame()
