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

county_df = read.csv("county_data.csv")

county_shapefile <- county_df$county_shp[[7]]
sample.design <- read.csv(county_df$sample_design[[7]])


# by about the 10th or 11th week, should have lots of ASF spreading...

solocs_dat <- read.csv("Output/solocs.all/solocs.all_r19_l4_v4.gz")
solocs_df <- as.data.frame(solocs_dat)
i = 11

bad_rows <- solocs_df %>%
            filter(pref == 0) %>%
            filter(S > 0)
 
fiscal_week_current <- ((i - 1) %% 52) + 1
# so here, rows_to_sample comes from real data
rows_to_sample <- sample.design[sample.design$fiscal_week == fiscal_week_current, ]