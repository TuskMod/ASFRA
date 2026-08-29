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

detection_folder <- "Output/Hancock-MS-detections/detections/"

file_paths <- list.files(path = detection_folder,
                         pattern= "\\.csv",
                         full.names = TRUE)
is_detected <- c()
time_detected <- c()
total_sampled <- c()
all_sampled <- c()
infec_sampled <- c()

hancock_ms_variables <- read.csv("hancock_ms_vars.csv")

plotIncidences <- function(file_paths,variables){
  
  hi_vars <- variables %>% filter(contact == "Hi")
  hi_vars <- hi_vars$var
  lo_vars <- variables %>% filter(contact == "Lo")
  lo_vars <- lo_vars$var
  print(hi_vars)
  print(lo_vars)
  lo_vars <- map(lo_vars,as.numeric())
  hi_vars <- map(hi_vars,as.numeric())
  infec_vals <- c()
  class_types <- c()
  weeks <- seq(1,78,1)
  all_times <- c()
  detect_weeks <- c()
  for (f in 1:length(file_paths)){
    file_data <- read.csv(file_paths[f])
    split_file <- strsplit(file_paths[f],"_")
    var_str <- split_file[[1]][[4]]
    var_list <- strsplit(var_str,"[.]")
    var_val <- as.numeric(substr(var_list[[1]][[1]],2,2))
    all_I <- rep(0,78)
    file_data <- file_data %>% filter(sample_type == 1)
    first_detect <- which(file_data$detections > 0)
    if(length(first_detect) > 0){
      detect_weeks <- c(detect_weeks,first_detect[[1]])
    }
    
    I_vals <- (file_data$I)
    if(length(I_vals) != 78){
     for (k in 1:length(I_vals)){
       all_I[k] <- I_vals[k]
     }
    }
    else{
      all_I <- I_vals
    }
    infec_vals <- c(infec_vals,all_I)
    if(var_val %in% lo_vars){
      class_types <- c(class_types,rep("Lo",78))
    }
    else{
      class_types <- c(class_types,rep("Hi",78))
    }
    all_times <- c(all_times,weeks)
  }
 
  average_week <- mean(detect_weeks)
 print(average_week)
 all_data <- data.frame(Incidence= infec_vals,
                        Time= all_times,
                        ContactRate = class_types
                        )
 ggplot(all_data,aes(x = Time,y=Incidence,color=ContactRate)) +
   stat_summary(geom = "point",fun="mean") +
   stat_summary(fun.data = "mean_se", geom="errorbar") +
   geom_vline(xintercept=round(average_week),linetype="dashed",color="red") + 
   ggtitle("Hi Contact Scenarios only ones that cause Outbreaks, average week of detection is 33.4 ")
 #return()
}

## plot effect of variances now just on Hi Contact

plotVariantIncidences <- function(file_paths,variables){
  
  hi_vars <- variables %>% filter(contact == "Hi")
  hi_vars <- hi_vars$var
  lo_vars <- variables %>% filter(contact == "Lo")
  lo_vars <- lo_vars$var
  print(hi_vars)
  print(lo_vars)
  lo_vars <- map(lo_vars,as.numeric())
  hi_vars <- map(hi_vars,as.numeric())
  infec_vals <- c()
  variant_type <- c()
  weeks <- seq(1,78,1)
  all_times <- c()
  detect_weeks_pol <- c()
  detect_weeks_dr <- c()
  all_vars <- c()
  for (f in 1:length(file_paths)){
    file_data <- read.csv(file_paths[f])
    split_file <- strsplit(file_paths[f],"_")
    var_str <- split_file[[1]][[4]]
    var_list <- strsplit(var_str,"[.]")
    var_val <- as.numeric(substr(var_list[[1]][[1]],2,2))
    if(var_val %in% lo_vars){
      next
    }
    all_I <- rep(0,78)
    file_data <- file_data %>% filter(sample_type == 1)
    first_detect <- which(file_data$detections > 0)
   
    
    if(var_val %in% c(1,2,5,6)){
      variant_type <- "Pol"
      if(length(first_detect) > 0){
        detect_weeks_pol <- c(detect_weeks_pol,first_detect[[1]])
      }
    }
    else{
      variant_type <- "DR"
      if(length(first_detect) > 0){
        detect_weeks_dr <- c(detect_weeks_dr,first_detect[[1]])
      }
      
    }
    variants <- rep(variant_type,78)
    
    I_vals <- (file_data$I)
    if(length(I_vals) != 78){
      for (k in 1:length(I_vals)){
        all_I[k] <- I_vals[k]
      }
    }
    else{
      all_I <- I_vals
    }
    infec_vals <- c(infec_vals,all_I)
  
    all_times <- c(all_times,weeks)
    all_vars <- c(all_vars,variants)
  }
  
  #average_week <- mean(detect_weeks)
  print(infec_vals)
  print(all_times)
  print(all_vars)
  
  all_data <- data.frame(Incidence= infec_vals,
                         Time= all_times,
                         Variant = all_vars
  )
  print(mean(detect_weeks_pol))
  print(mean(detect_weeks_dr))
  ggplot(all_data,aes(x = Time,y=Incidence,color=(Variant))) +
    stat_summary(geom = "point",fun="mean") +
    scale_color_manual(values=c("#B10DC9","#85144b"))+
    stat_summary(fun.data = "mean_se", geom="errorbar") +
    geom_vline(xintercept=round(mean(detect_weeks_dr)),linetype="dashed",color="red") + 
    geom_vline(xintercept=round(mean(detect_weeks_pol)),linetype="dashed",color="blue") +
    ggtitle("Hi contact only: DR Variant detected more quickly than Poland Variant ")
  #return()
}

plotDetectionProbs <- function (file_paths,variables){
  
  hi_vars <- variables %>% filter(contact == "Hi")
  hi_vars <- hi_vars$var
  lo_vars <- variables %>% filter(contact == "Lo")
  lo_vars <- lo_vars$var
  print(hi_vars)
  print(lo_vars)
  lo_vars <- map(lo_vars,as.numeric())
  hi_vars <- map(hi_vars,as.numeric())
  infec_vals <- c()
  variant_type <- c()
  weeks <- seq(1,78,1)
  all_weeks <- c()
  detect_weeks_pol <- c()
  detect_weeks_dr <- c()
  detect_probs <- c()
  all_vars <- c()
  all_types <- c()
  all_densities <- c()
  for (f in 1:length(file_paths)){
    file_data <- read.csv(file_paths[f])
    split_file <- strsplit(file_paths[f],"_")
    var_str <- split_file[[1]][[4]]
    var_list <- strsplit(var_str,"[.]")
    var_val <- as.numeric(substr(var_list[[1]][[1]],2,2))
    if(var_val %in% lo_vars){
      next
    }
    file_data <- file_data %>% filter(sample_type == 1)
    first_detect <- which(file_data$detections > 0)
    
    N <- file_data$S + file_data$E + file_data$I + file_data$R + file_data$C +file_data$Z
    
    detection_probability <- 1- ((file_data$S + file_data$R) / N)
    all_detect <- rep(0,78)
    all_N <- rep(0,78)
    all_I <- rep(0,78)
    all_Z <- rep(0,78)
    all_C <- rep(0,78)
    all_R <- rep(0,78)
    if(N[1] > 20000){
      density <- rep(3,78)
      print(detection_probability)
    }
    else{
      density <- rep(1.5,78)
    }
    if(length(detection_probability) != 78){
      for (k in 1:length(detection_probability)){
        all_detect[k] <- detection_probability[k]
        all_N[k] <- N[k]/N[1]
        all_I[k] <- (file_data$I[k] + file_data$E[k])/N[1]
        all_C[k] <- (file_data$C[k])/N[1]
        all_Z[k] <- (file_data$Z[k])/N[1]
        all_R[k] <- (file_data$R[k])/N[1]
      }
    }
    else{
      all_detect <- detection_probability
      all_N <- N/N[1]
      all_I <- file_data$I/N[1]
      all_C <- file_data$C/N[1]
      all_Z <- file_data$Z/N[1]
      all_R <- file_data$R/N[1]
    }
    if(var_val %in% c(1,2,5,6)){
      variant_type <- "Pol"
      if(length(first_detect) > 0){
        detect_weeks_pol <- c(detect_weeks_pol,first_detect[[1]])
      }
    }
    else{
      variant_type <- "DR"
      if(length(first_detect) > 0){
        detect_weeks_dr <- c(detect_weeks_dr,first_detect[[1]])
      }
      
    }
    variants <- rep(variant_type,78)
    all_vars <- c(all_vars,variants,variants,variants,variants,variants,variants)
    all_weeks <- c(all_weeks,weeks,weeks,weeks,weeks,weeks,weeks)
    all_types <- c(all_types,rep("D",78),rep("I+E",78),rep("R",78),rep("C",78),rep("Z",78),rep("N",78))
    all_densities <- c(all_densities,density,density,density,density,density,density)
    detect_probs <- c(detect_probs,all_detect,all_I,all_R,all_C,all_Z,all_N)
  }
  
  all_data <- data.frame(Probability= detect_probs,
                         Time= all_weeks,
                         Class= all_types,
                         Strain = all_vars,
                         Density = all_densities)
  ggplot(all_data,aes(x = Time,y=Probability,color=Class)) +
    stat_summary(geom = "point",fun="mean") +
    stat_summary(fun.data = "mean_se", geom="errorbar") +
    facet_grid(Density ~ Strain, axes="all", axis.labels = "all_x") +
    geom_vline(xintercept=round(mean(detect_weeks_dr)),linetype="dashed",color="red") + 
   geom_vline(xintercept=round(mean(detect_weeks_pol)),linetype="dashed",color="blue") 
  # ggtitle("Highest probability of detection occurs after Incidence Peak ")
  
}




plotIncidences(file_paths,hancock_ms_variables)
plotVariantIncidences(file_paths,hancock_ms_variables)
plotDetectionProbs(file_paths,hancock_ms_variables)
#plotInfectionCurvesOverTime(file_paths,hancock_ms_variables)

for (f in 1:length(file_paths)){
  #split_file <- strsplit(file_paths[f],"_")
  #if(split_file[[1]][[4]] == "v8.csv"){
 
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
    time_detected <- c(time_detected,time_val)
  }
#}



print("Prob Detected: ")
print(sum(is_detected)/length(is_detected))
print("Average Time to detection: ")
print(mean(time_detected))
print("Num Sampled :")
print(sum(total_sampled)/length(total_sampled))
print("Infected Sampled")
print(sum(infec_sampled)/length(infec_sampled))
title_string <- paste0("Average Week Detected")
#hist(time_detected,xlab="Week Detected")
#main = (paste("Average Week Detected: ",round(mean(time_detected),2))))
#hist(total_sampled)
#hist(infec_sampled)
