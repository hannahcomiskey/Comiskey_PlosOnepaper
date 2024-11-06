country_subnat_tbl <- FP_source_data_wide %>% 
  group_by(Country, Region) %>% 
  dplyr::select(Country, Region) %>% 
  distinct() # table of country and regions (repeats in region names)
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix
n_country <- unique(country_subnat_tbl$Country)
n_subnat <- country_subnat_tbl$Region

FP_source_data_wide <- subnat_index_fun(FP_source_data_wide, n_subnat, country_subnat_tbl$Country)
FP_source_data_wide <- country_index_fun(FP_source_data_wide, n_country)
FP_source_data_wide <- method_index_fun(FP_source_data_wide, n_method)

# Time indexing - important for splines -----------------------------------
all_years <- seq(from = 1990, to = 2023.5, by=0.5) # shorter due to memory issues
n_all_years <- length(all_years)

FP_source_data_wide <- FP_source_data_wide %>%
  mutate(index_year = match(average_year,all_years)) 

#################################################
# setup for JAGS data ---------------------------
#################################################

t_seq_2 <- floor(FP_source_data_wide$index_year) # Time sequence for countries
country_seq <- FP_source_data_wide$Country
n_sector <- c("Public", "Commercial_medical", "Other") # Names of categories
n_obs <- nrow(FP_source_data_wide) # Total number of observations
year_seq <- seq(min(t_seq_2),max(t_seq_2), by=1)
n_years <- length(year_seq)

# Find the observation year indexes in the prediction years
country_index_tbl <- FP_source_data_wide %>% 
  group_by(Country, index_country) %>% 
  dplyr::select(Country, index_country) %>%
  distinct() # table of country and regions (repeats in region names)

index_country_subnat_tbl <- FP_source_data_wide %>% 
  group_by(index_country, index_subnat) %>% 
  dplyr::select(Region, index_subnat, Country, index_country) %>% #Super_region, index_superregion) %>%
  distinct() # table of country and regions (repeats in region names)

count_provinces <- index_country_subnat_tbl %>% # count number of provinces in each country
  group_by(index_country) %>% 
  count(index_country) %>%
  dplyr::select(index_country, n) %>%
  distinct()

match_country <- index_country_subnat_tbl$index_country

match_years <- FP_source_data_wide$index_year

match_method <- FP_source_data_wide$index_method

match_subnat <-  FP_source_data_wide$index_subnat

# Get T_star and match_Tstar -----------------------------
T_star <- FP_source_data_wide %>%
  group_by(Country, Region) %>%
  dplyr::filter(index_year==max(index_year)) %>%
  dplyr::select(Country, Region, index_country, index_subnat, average_year, index_year) %>%
  arrange(index_subnat) %>%
  ungroup() %>%
  dplyr::select(index_country, index_subnat, average_year, index_year) %>%
  distinct()
