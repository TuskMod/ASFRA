#load libraries for targets script
library(targets)
library(tarchetypes)
# library(crew)
# library(crew.cluster)
library(clustermq)
library(Rcpp)
library(data.table)
library(geosphere)
library(stringr)
library(tidyverse)
library(dplyr)
library(sf)
library(terra)

## for each variable type output the following
# Detection pronanility
# Time to detection
# Number of 

# What do we want to see?
# I think comparing the effectiveness of different surveillance strateiges
# ....I think we need a comparison between the homog vs heterog landscape one day...
# For now, let us just replicate parts of Madison's work
# I think the stat modeling is the next part...

detection_folder <- "Output/detections/"

file_paths <- list.files(path = detection_folder,
                         pattern= "\\.csv",
                         full.names = TRUE)
is_detected <- c()
time_detected <- c()
total_sampled <- c()
all_sampled <- c()
infec_sampled <- c()

for (f in 1:length(file_paths)){
  split_file <- strsplit(file_paths[f],"_")
  if(split_file[[1]][[4]] == "v6.csv"){
    file_data <- read.csv(file_paths[f])
    total_sampled <-c(total_sampled,sum(file_data$detections))
    detect_df <- file_data %>% filter(sample_type == 1)
    detect_sum <- sum(detect_df$detections)
    infec_sampled <- c(infec_sampled,detect_sum)
    if(detect_sum > 0){
      is_detected <- c(is_detected,1)
    }
    else{
      is_detected <-c(is_detected,0) 
    }
    just_detect_instances <- detect_df %>% filter(detections > 0)
    if (length(just_detect_instances) >  0){
      time_val <- head(just_detect_instances,n=1)$timestep
    }
    else{
      time_val <- 0
    }
    time_detected <- c(time_detected,time_val)
  }
}

print("Prob Detected: ")
print(sum(is_detected)/length(is_detected))
print("Average Time to detection: ")
print(mean(time_detected))
print("Num Sampled :")
print(sum(total_sampled)/length(total_sampled))
print("Infected Sampled")
print(sum(infec_sampled)/length(infec_sampled))
title_string <- paste0("Average Week Detected")
hist(time_detected,xlab="Week Detected")
#main = (paste("Average Week Detected: ",round(mean(time_detected),2))))
hist(total_sampled)
hist(infec_sampled)
