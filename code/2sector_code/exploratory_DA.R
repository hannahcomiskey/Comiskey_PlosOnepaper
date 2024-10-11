library(mclust)
library(factoextra)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

# Source code --------------------------------------
source('code/load_functions.R')
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# Set up rates of change -------------
lagged_data <- FP_source_data_wide %>%
  left_join(area_classification) %>%
  group_by(Country, Region, Method) %>%
  arrange(average_year) %>%
  mutate(logit.public = log(Public/(1-Public)),
         lag.logit.public = lag(logit.public),
         rate_of_change.public = logit.public - lag.logit.public,
         logit.private = log(Private/(1-Private)),
         lag.logit.private = lag(logit.private),
         rate_of_change.private = logit.private - lag.logit.private) %>%
  select(!c(Public_n, Private_n, index_subnat, index_country, index_method, index_year, check_total))

# Create dataframe for clustering
cluster_df <- lagged_data %>% 
  select(!c(Public, Private, logit.public, logit.private, lag.logit.public, lag.logit.private)) %>% 
  group_by(Country, Region, Method) %>%
  filter(average_year==max(average_year)) %>% 
  filter(is.na(rate_of_change.public)==FALSE | is.na(rate_of_change.private)==FALSE)

glimpse(cluster_df)

plotting_data <- cluster_df %>%
  pivot_longer(cols = c(rate_of_change.public, rate_of_change.private), names_to = 'Rate_type', values_to = 'value')


ggplot(plotting_data) +
  geom_density(aes(x=value, fill=Rate_type), alpha=0.5) +
  facet_wrap(~Method, scales='free')

glimpse(plotting_data)
