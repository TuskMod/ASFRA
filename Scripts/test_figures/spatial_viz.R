library(tidyverse)
library(ggplot2)
plotSpatialSpread <- function (file_paths,variables){
  
  hi_vars <- variables %>% filter(contact == "Hi")
  hi_vars <- hi_vars$var
  lo_vars <- variables %>% filter(contact == "Lo")
  lo_vars <- lo_vars$var

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
  cell_x <- rep(0,40000)
  cell_y <- rep(seq(1,200,1),200)
  counter <- 1
  for (x in 1:200){
    for(y in 1:200){
      cell_x[counter] <- x
      counter <- counter + 1
    }
  }
  counts <- 0
  # let's only analyze files which have full time lengths...
  for (f in 1:length(file_paths)){
    file_data <- read.csv(file_paths[f])
    #file_data <- file_data %>% filter(sample_type == 1)
    cell_vals <- rep(0,200*200)
    
    split_file <- strsplit(file_paths[f],"_")
    var_str <- split_file[[1]][[4]]
    var_list <- strsplit(var_str,"[.]")
    var_val <- as.numeric(substr(var_list[[1]][[1]],2,2))
    time_amount <- unique(file_data$time)
    if(var_val %in% lo_vars | length(time_amount) != 78){
      next
    }
    counts = counts + 1
    for(k in 1:length(time_amount)){
      relv_data <- file_data %>% filter(time == k) %>% filter(I > 0)
      infec_cells <- relv_data$cell
      print(infec_cells)
      for(i in 1:length(infec_cells)){
        print(infec_cells[i])
        cell_vals[infec_cells[i]] <- cell_vals[infec_cells[i]] + 1
      }
      new_data_frame <- data.frame(X = cell_x,
                                   Y = cell_y,
                                   Infec = cell_vals)
      ggplot(new_data_frame, aes(x = X, y= Y,fill = Infec)) +
        geom_tile()
      save_name <- paste0("Figures/spat_maps/",k,"_",counts,".png")
      ggsave(save_name)
    }
    
  
    }
}


detection_folder <- "Output/solocs.all/"

file_paths <- list.files(path = detection_folder,
                         pattern= "\\.gz",
                         full.names = TRUE)
is_detected <- c()
time_detected <- c()
total_sampled <- c()
all_sampled <- c()
infec_sampled <- c()

hancock_ms_variables <- read.csv("MS-hancock-vars.csv")

plotSpatialSpread(file_paths, hancock_ms_variables)