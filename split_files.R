# HEADER ---------------------------------------------------------------------------------------------------------------

# This script loads the files provided for the competition and splits their contents into smaller files which are 
# small enough to be uploaded to GitHub (< 10MB). It then loads the contents from the smaller files and compares
# them to the originals, to make sure no sample has been lost.

# LOAD LIBRARIES AND READ ORIGINAL FILES -------------------------------------------------------------------------------

# load required libraries

library("tidyverse")
library("readr")
library("fs")
library("purrr")

# read original files

files <-
  list(train_features = read_csv(file.path("files","train_features.csv")),
       train_targets_scored = read_csv(file.path("files","train_targets_scored.csv")),
       train_targets_nonscored = read_csv(file.path("files","train_targets_nonscored.csv")),
       test_features = read_csv(file.path("files","test_features.csv"))) 

# SPLIT CONTENT OF ORIGINAL FILES INTO MULTIPLE SMALLER FILES ----------------------------------------------------------

# split contents of train_features.csv file into multiple smaller files

numPartitions = 16

for (i in 1:numPartitions) {
  
  fileName = paste("train_features_",
                   str_pad(i, side = "left", width = 2, pad = "0"),
                   ".csv", sep = "")
  
  files$train_features %>% 
    mutate(fileNum = sort(rep(x = 1:numPartitions, 
                              length.out = nrow(files$train_features)))) %>% 
    filter(fileNum == i) %>% 
    select(-fileNum) %>% 
    write_csv(file = file.path("files",fileName),
              append = FALSE)
}

# split contents of test_features.csv file into multiple smaller files

numPartitions = 2

for (i in 1:numPartitions) {
  
  fileName = paste("test_features_",
                   str_pad(i, side = "left", width = 2, pad = "0"),
                   ".csv", sep = "")
  
  files$test_features %>% 
    mutate(fileNum = sort(rep(x = 1:numPartitions, 
                              length.out = nrow(files$test_features)))) %>% 
    filter(fileNum == i) %>% 
    select(-fileNum) %>% 
    write_csv(file = file.path("files",fileName),
              append = FALSE)
}

# split contents of train_targets_scored.csv file into multiple smaller files

numPartitions = 2

for (i in 1:numPartitions) {
  
  fileName = paste("train_targets_scored_",
                   str_pad(i, side = "left", width = 2, pad = "0"),
                   ".csv", sep = "")
  
  files$train_targets_scored %>% 
    mutate(fileNum = sort(rep(x = 1:numPartitions, 
                              length.out = nrow(files$train_targets_scored)))) %>% 
    filter(fileNum == i) %>% 
    select(-fileNum) %>% 
    write_csv(file = file.path("files",fileName),
              append = FALSE)
}

# split contents of train_targets_nonscored.csv file into multiple smaller files

numPartitions = 2

for (i in 1:numPartitions) {
  
  fileName = paste("train_targets_nonscored_",
                   str_pad(i, side = "left", width = 2, pad = "0"),
                   ".csv", sep = "")
  
  files$train_targets_nonscored %>% 
    mutate(fileNum = sort(rep(x = 1:numPartitions, 
                              length.out = nrow(files$train_targets_nonscored)))) %>% 
    filter(fileNum == i) %>% 
    select(-fileNum) %>% 
    write_csv(file = file.path("files",fileName),
              append = FALSE)
}

# LOAD REDUCED FILES AND COMPARE CONTENTS TO THE ORIGINAL ONES ---------------------------------------------------------

# is the combination of reduced files identical to the original?

identical(files$train_features,
          dir_ls("files") %>% 
            str_subset("train_features_\\d{2}") %>% 
            map_dfr(read_csv, col_types = cols()))

identical(files$test_features,
          dir_ls("files") %>% 
            str_subset("test_features_\\d{2}") %>% 
            map_dfr(read_csv, col_types = cols()))

identical(files$train_targets_scored,
          dir_ls("files") %>% 
            str_subset("train_targets_scored_\\d{2}") %>% 
            map_dfr(read_csv, col_types = cols()))

identical(files$train_targets_nonscored,
          dir_ls("files") %>% 
            str_subset("train_targets_nonscored_\\d{2}") %>% 
            map_dfr(read_csv, col_types = cols()))
