library(spdep)
library(dplyr)
library(sf)
library(ggplot2)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

n_country_tmp <- n_country 
n_country_tmp[which(n_country_tmp=="Cote d'Ivoire")] <- "Côte D'Ivoire"

# Read in data
my_sf <- read_sf("data/shape/world_geolev1_2024/world_geolev1_2024.shp") %>%
  filter(CNTRY_NAME %in% n_country_tmp)

my_sf <- st_make_valid(my_sf)

ggplot(my_sf) +
  geom_sf(aes(fill=CNTRY_NAME), show.legend = FALSE) +
  theme_void()

# Removing Afghanistan, Niger and Madagascar as they do not have data available
countries <- unique(my_sf$CNTRY_NAME)[-c(1, 12, 17)]
queen_list <- list()
for(i in 1:length(countries)) {
  tmp_sf <- my_sf %>% filter(CNTRY_NAME==countries[i])
  tmp_sf <- st_make_valid(tmp_sf) %>% filter(!st_is_empty(.))
  queen_hoods <- st_geometry(tmp_sf) %>% 
    poly2nb(queen = TRUE) 
  queen_list[[countries[i]]] <- queen_hoods
}

mm_list <- list()
for(i in 1:length(countries)) {
  mm = nb2mat(queen_list[[i]],zero.policy=TRUE, style="B")
  mm_list[[countries[i]]] <- mm
}



