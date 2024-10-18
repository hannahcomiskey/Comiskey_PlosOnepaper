##############################################
# Read in SE data
##############################################
subnat_FPsource_data <- readRDS("subnat_bivar_SE_source_data_20.RDS") # All countries

subnat_FPsource_data <- subnat_FPsource_data %>% 
  ungroup() %>%
  dplyr::mutate(Region = stringr::str_to_title(Region))

###############################################
# Updated cleaning approach -------------------
###############################################
subnat_FPsource_data <- subnat_FPsource_data %>%
  dplyr::filter(Region!="NA")

# Filter where the sample size is smaller than 2 people across all three sectors ---------

FP_source_data_wide <- subnat_FPsource_data %>% # Proportion data
  dplyr::ungroup() %>%
  dplyr::rename(Public.SE = se.Public, Private.SE = se.Private) %>%
  dplyr::select(Country, Region, Method,  average_year, Public, Private, Public.SE, Private.SE, Public_n, Private_n) %>%
  dplyr::arrange(Country) %>%
  dplyr::filter(Public_n>=10 | Private_n>=10) 

# Make sure proportions add to 1 ----------
FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::rowwise() %>%
  dplyr::mutate(check_total = sum(Public, Private, na.rm = TRUE))

col_index <- which(colnames(FP_source_data_wide)=="Public")-1 # column index before CM column, as CM=1
# When check_Total=1, replace missing values with 0 ----------
for (i in 1:nrow(FP_source_data_wide)) {
  if(FP_source_data_wide$check_total[i]>0.99) {
    na_cols <- which(is.na(FP_source_data_wide[i, c("Public", "Private")])==TRUE) # NA values
    FP_source_data_wide[i, na_cols+col_index] <- as.list(rep(0, length(na_cols)))
  }
}

# Transform exactly 1 and 0 values away from boundary using lemon-squeezer approach ---------
FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::mutate(Public = (Public*(nrow(FP_source_data_wide)-1)+0.5)/nrow(FP_source_data_wide)) %>%
  dplyr::mutate(Private = (Private*(nrow(FP_source_data_wide)-1)+0.5)/nrow(FP_source_data_wide)) %>%
  dplyr::select(Country, Region, Method, average_year, Public, Private, Public.SE, Private.SE, Public_n, Private_n, check_total) #, count_NA, remainder)

# Clean SE values --------------
FP_source_data_wide$count_SE.NA <- rowSums(is.na(FP_source_data_wide %>% dplyr::select(Public.SE, Private.SE))) # count NAs

SE_source_data_wide_norm <- FP_source_data_wide %>% dplyr::filter(count_SE.NA==0 & Private.SE>=0.005 & Public.SE>=0.005) # Normal obs. No action needed.
SE_source_data_wide_X <- FP_source_data_wide %>% dplyr::filter(Private.SE<0.005 | Public.SE<0.005) # Get obs with two missing sectors

# Replacing SE=0 and no NAs: Leontine's suggestion 
SE_source_data_wide_X <- SE_source_data_wide_X %>%
  dplyr::filter(Public_n>=20 | Private_n>=20) # Remove small sample sizes (DHS has 10 units sampled per cluster as min., 20 as average)
col_index <- which(colnames(SE_source_data_wide_X)=="Public.SE")-1 # column index before CM column, as CM=1
DEFT_data <- readxl::read_xlsx("DEFT_DHS_database.xlsx") %>%
  rename(average_year = Year)
SE_source_data_wide_X <- SE_source_data_wide_X %>% left_join(DEFT_data)

# https://onlinestatbook.com/2/sampling_distributions/samp_dist_p.html

for(i in 1:nrow(SE_source_data_wide_X)) {
  num.SE0 <- which(SE_source_data_wide_X[i,c("Public.SE", "Private.SE")]<0.005)
  num.SEna <- which(is.na(SE_source_data_wide_X[i,c("Public.SE", "Private.SE")])==TRUE)
  DEFT <- ifelse(is.na(SE_source_data_wide_X$DEFT[i])==TRUE, 1.5, SE_source_data_wide_X$DEFT[i])
  N1 <- sum(SE_source_data_wide_X[i, c('Public_n', 'Private_n')], na.rm=TRUE) # Number of women surveyed
  phat <- 0.5/(N1+1) # Posterior mean of p under Jefferys prior for true prevalence of 0s.
  SE.hat <- sqrt((phat*(1-phat))/N1)
  SE_source_data_wide_X[i,c(col_index+num.SEna,col_index+num.SE0)] <- SE.hat*DEFT
}

FP_source_data_wide <- bind_rows(SE_source_data_wide_norm, SE_source_data_wide_X) # Put data back together again

# Remove proportions with two sectors still missing 
FP_source_data_wide <- FP_source_data_wide %>% 
  filter(is.na(Public)==FALSE)

FP_source_data_wide <- FP_source_data_wide %>% dplyr::arrange(Country, Region, Method, average_year)

#################################################
# Issues with subnational names -----------------
#################################################
# Replace issues with Burkina Faso names
FP_source_tmp <- FP_source_data_wide %>%
  dplyr::filter(Country=="Burkina Faso")  %>%
  dplyr::mutate(Region = case_when(Region == "North" ~ "Nord",
                                   Region == "East" ~ "Est",
                                   Region == "West" ~ "Centre-Ouest",
                                   Region ==  "Central/South" ~ "Centre-Sud",
                                   Region == "Ouagadougou" ~ "Centre Including Ouagadougou",
                                   TRUE ~ as.character(Region)))

FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::filter(Country!="Burkina Faso")

FP_source_data_wide <- FP_source_data_wide %>%
  merge(FP_source_tmp, all = TRUE)
