####################################
# Using only FPET countries for now 
####################################
area_classification <- mcmsupply::Country_and_area_classification %>% 
  dplyr::select(`Country or area`, Region) %>% 
  dplyr::rename(Country = `Country or area`) %>%
  dplyr::rename(Super_region = Region) %>%
  dplyr::mutate(Country = case_when(Country== "Bolivia (Plurinational State of)" ~ "Bolivia",
                                    Country== "Republic of Moldova" ~ "Moldova",
                                    Country== "Viet Nam" ~ "Vietnam",
                                    Country== "Kyrgyzstan" ~ "Kyrgyz Republic",
                                    .default = as.character(Country)))

FP_source_data_wide <- FP_source_data_wide %>%
  left_join(area_classification) 

# Time indexing - important for splines -----------------------------------
all_years <- seq(from = 1990, to = 2030.5, by=0.5) # shorter due to memory issues
n_all_years <- length(all_years)

FP_source_data_wide <- FP_source_data_wide %>%
  mutate(index_year = match(average_year,all_years))


#################################################
# Get test and training data -------------------
#################################################

training_data <- FP_source_data_wide %>% 
  filter(average_year <= 2015)

country_subnat_tbl <- training_data %>% 
  group_by(Country, Region) %>% 
  dplyr::select(Country, Region) %>% 
  distinct() # table of country and regions (repeats in region names)
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix
n_country <- unique(country_subnat_tbl$Country)
n_subnat <- country_subnat_tbl$Region

test_data <- FP_source_data_wide %>% 
  filter(average_year > 2015 & Country %in% n_country & Region %in% n_subnat)


#################################################
# Adding indexes to  data ----------------------
#################################################

training_data <- subnat_index_fun(training_data, n_subnat, country_subnat_tbl$Country)
training_data <- country_index_fun(training_data, n_country)
training_data <- method_index_fun(training_data, n_method)

test_data <- subnat_index_fun(test_data, n_subnat, country_subnat_tbl$Country)
test_data <- country_index_fun(test_data, n_country)
test_data <- method_index_fun(test_data, n_method)


#################################################
# setup for JAGS data ---------------------------
#################################################

t_seq_2 <- floor(training_data$index_year) # Time sequence for countries
n_sector <- c("Public", "Private") # Names of categories
n_obs <- nrow(training_data) # Total number of observations
year_seq <- seq(min(t_seq_2),max(t_seq_2), by=1)
n_years <- length(year_seq)

# Find the observation year indexes in the prediction years
# country_index_tbl <- FP_source_data_wide %>% 
#   group_by(Country, index_country) %>% 
#   dplyr::select(Country, index_country) %>%
#   distinct() # table of country and regions (repeats in region names)

index_country_subnat_tbl <- training_data %>% 
  group_by(index_country, index_subnat) %>% 
  dplyr::select(Region, index_subnat, Country, index_country) %>%
  distinct() # table of country and regions (repeats in region names)

count_provinces <- index_country_subnat_tbl %>% # count number of provinces in each country
  group_by(index_country) %>% 
  count(index_country) %>%
  dplyr::select(index_country, n) %>%
  distinct()

match_country <- index_country_subnat_tbl$index_country

match_years <- training_data$index_year

match_method <- training_data$index_method

match_subnat <-  training_data$index_subnat

# Get T_star and match_Tstar -----------------------------
T_star <- training_data %>%
  group_by(Country, Region) %>%
  dplyr::filter(index_year==max(index_year)) %>%
  dplyr::select(Country, Region, index_country, index_subnat, average_year, index_year) %>%
  arrange(index_subnat) %>%
  ungroup() %>%
  dplyr::select(index_country, index_subnat, average_year, index_year) %>%
  distinct()

#################################################
# Splines --------------------------------------
#################################################
nseg=10
Kstar <- vector()
B.ik <- array(dim = c(length(n_subnat), length(all_years),nseg+3))
knots.all <- matrix(nrow=length(n_subnat), ncol=nseg+3)

for(i in 1:nrow(T_star)) {
  index_p <- T_star$index_subnat[i]
  index_msn <- T_star$average_year[i]
  res <- bs_bbase_precise(all_years, lastobs=index_msn, nseg = nseg) 
  B.ik[index_p,,] <- res$B.ik
  Kstar[index_p] <- res$Kstar
  knots.all[index_p,] <- res$knots.k
}

K <- dim(res$B.ik)[2]
H <- K-1

## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,res$B.ik[,1], type= "n", xaxt="n",
     xlab = "Year",
     ylim = c(0,1), ylab ="Basis Function",
     xlim = range(all_years))
axis(1, at = min(all_years):max(all_years))
abline(v=res$knots.k, col = seq(1, K), lwd = 1)
for (k in 1:K){
  lines(all_years,res$B.ik[,k], type= "l", col = k, lwd = 1)
}


