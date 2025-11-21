##############################################
# Read in SE data
##############################################
load("data/subnat_bivar_data.rda") # All countries

FP_2030_countries <- c(
  # "Afghanistan", 
  "Benin", "Burkina Faso", "Cameroon",
  "Cote d'Ivoire", "Ethiopia", "Ghana",
  "Guinea", "Kenya", "Liberia", "Madagascar",
  "Malawi", "Mali", "Mozambique", "Myanmar",
  "Nepal", "Niger", "Nigeria", "Pakistan",
  "Rwanda", "Senegal", "Tanzania", "Uganda",
  "Zimbabwe"
)

subnat_bivar_data <- subnat_bivar_data %>% 
  filter(Country %in% FP_2030_countries) %>%
  ungroup() %>%
  dplyr::mutate(Region = stringr::str_to_title(Region))

###############################################
# Updated cleaning approach -------------------
###############################################
subnat_bivar_data <- subnat_bivar_data %>%
  dplyr::filter(Region!="NA")

# Filter where the sample size is smaller than 10 people across all three sectors ---------

FP_source_data_wide <- subnat_bivar_data %>% # Proportion data
  dplyr::ungroup() %>%
  dplyr::rename(Public.SE = se.Public, Private.SE = se.Private) %>%
  dplyr::select(Country, Region, Method,  average_year, Public, Private, Public.SE, Private.SE, Public_n, Private_n) %>%
  dplyr::arrange(Country) %>%
  dplyr::filter(Public_n>=20 | Private_n>=20) # Remove small sample sizes (DHS has 10 units sampled per cluster as min., 20 as average)

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
FP_source_data_wide$count_SE.NA <- rowSums(is.na(FP_source_data_wide %>% dplyr::select(Public.SE))) # count NAs

SE_source_data_wide_norm <- FP_source_data_wide %>% dplyr::filter(count_SE.NA==0 & Public.SE>0.005) # Normal obs. No action needed.
SE_source_data_wide_X <- FP_source_data_wide %>% dplyr::filter(count_SE.NA>0 | Public.SE<=0.005) # Get obs with two missing sectors

# Replacing SE=0 and no NAs: Leontine's suggestion 
SE_source_data_wide_X <- SE_source_data_wide_X %>%
  dplyr::filter(Public_n>=20 | Private_n>=20) # Remove small sample sizes (DHS has 10 units sampled per cluster as min., 20 as average)
col_index <- which(colnames(SE_source_data_wide_X)=="Public.SE")-1 # column index before CM column, as CM=1
load("data/DEFT_DHS_database.rda") 
DEFT_data <- DEFT_DHS_database %>%
  dplyr::rename(average_year = Year)
SE_source_data_wide_X <- SE_source_data_wide_X %>%  dplyr::left_join(DEFT_data)

# # https://onlinestatbook.com/2/sampling_distributions/samp_dist_p.html

if(nrow(SE_source_data_wide_X)>0) {
  for(i in 1:nrow(SE_source_data_wide_X)) {
    num.SE0 <- which(SE_source_data_wide_X[i,c("Public.SE")]==0)
    num.SEnon0 <- which(SE_source_data_wide_X[i,c("Public.SE")]!=0)
    num.SEverytiny <- which(SE_source_data_wide_X[i,c("Public.SE")]>0 & SE_source_data_wide_X[i,c("Public.SE")]<0.001)
    num.SEna <- which(is.na(SE_source_data_wide_X[i,c("Public.SE")])==TRUE)
    DEFT <- ifelse(is.na(SE_source_data_wide_X$DEFT[i])==TRUE, 1.5, SE_source_data_wide_X$DEFT[i])
    N1 <- sum(SE_source_data_wide_X[i, c('Public_n', 'Private_n')], na.rm=TRUE) # Number of women surveyed
    phat <- 0.5/(N1+1) # Posterior mean of p under Jefferys prior for true prevalence of 0s.
    SE.hat <- sqrt((phat*(1-phat))/N1)
    SE_source_data_wide_X[i,c(col_index+num.SEna, col_index+num.SE0, col_index+num.SEverytiny)] <- SE.hat*DEFT
    
  }
  FP_source_data_wide <- bind_rows(SE_source_data_wide_norm, SE_source_data_wide_X) # Put data back together again
}


# Remove proportions with two sectors still missing 
FP_source_data_wide <- FP_source_data_wide %>% 
  filter(is.na(Public)==FALSE)

FP_source_data_wide <- FP_source_data_wide %>% dplyr::arrange(Country, Region, Method, average_year)


# #################################################
# # Issues with subnational names -----------------
# #################################################
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

# # Replace issues with Rwanda names
FP_source_tmp <- FP_source_data_wide %>%
  dplyr::filter(Country=="Rwanda") %>%
  dplyr::mutate(Region = case_when(Region == "Ouest" ~ "West",
                                   Region == "Nord" ~ "North",
                                   Region == "Sud" ~ "South",
                                   Region == "Est" ~ "East",
                                   Region == "Ville de Kigali" ~ "Kigali",
                                   Region == "Ville De Kigali" ~ "Kigali",
                                   Region == "Kigali City" ~ "Kigali",
                                   Region == "Butare, Gitarama (Central, South)" ~ "South",
                                   Region == "Cyangugu, Gikongoro (Southwest)" ~ "South",
                                   Region == "Byumba, Kibungo, Umutara (Northeast)" ~ "North",
                                   Region == "Gisenyi, Kibuye, Ruhengeri (Northwest)" ~ "West",
                                   TRUE ~ as.character(Region)))

FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::filter(Country!="Rwanda")

FP_source_data_wide <- FP_source_data_wide %>%
  merge(FP_source_tmp, all = TRUE)

# Replace issues with Nigeria names
FP_source_tmp <- FP_source_data_wide %>%
  dplyr::filter(Country=="Nigeria") %>%
  dplyr::mutate(Region = case_when(Region == "Northeast" ~ "North East",
                                   Region == "Northwest" ~ "North West",
                                   Region == "Southeast" ~ "South East",
                                   Region == "Southwest" ~ "South West",
                                   Region == "Central" ~ "North Central",
                                   TRUE ~ as.character(Region)))

FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::filter(Country!="Nigeria")

FP_source_data_wide <- FP_source_data_wide %>%
  merge(FP_source_tmp, all = TRUE)

# Replace issues with Cote d'Ivoire names
FP_source_tmp <- FP_source_data_wide %>%
  dplyr::filter(Country=="Cote d'Ivoire") %>%
  dplyr::mutate(Region = case_when(Region == "Center-East" ~ "Center East",
                                   Region == "Center-North" ~ "Center North",
                                   Region == "Center-West" ~ "Center West",
                                   Region == "Center-South" ~ "Center South",
                                   TRUE ~ as.character(Region)))

FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::filter(Country!="Cote d'Ivoire")

FP_source_data_wide <- FP_source_data_wide %>%
  merge(FP_source_tmp, all = TRUE)


#################################################
# Plot raw data -------------------------------
#################################################

# tmp1 <- FP_source_data_wide %>% dplyr::filter(Country=="Kenya") %>%
#   select(Country, Method, Region, average_year, Public, Private) %>%
#   pivot_longer(cols =c(Public, Private), names_to = "sector_category", values_to = "proportion")
# 
# tmp2 <- FP_source_data_wide %>% dplyr::filter(Country=="Kenya") %>%
#   select(Country, Method, Region, average_year, Public.SE, Private.SE) %>%
#   pivot_longer(cols =c(Public.SE, Private.SE), names_to = "sector_category", values_to = "SE_proportion") %>%
#   mutate(sector_category = substr(sector_category, 1, nchar(sector_category)-3))
# 
# plot_data <- merge(tmp1, tmp2) %>%
#   mutate(prop_min = ifelse(proportion - 2*SE_proportion<0, 0, proportion - 2*SE_proportion)) %>%
#   mutate(prop_max = ifelse(proportion + 2*SE_proportion>1, 1, proportion + 2*SE_proportion))
# 
# unique(plot_data$Region)
# 
# ggplot(data = plot_data %>% filter(Region=="Nairobi")) +
#   labs(y="Proportion of contraceptives supplied", x = "Year", title = paste0(unique(tmp2$Region),",",unique(tmp2$Country))) +
#   geom_point(aes(x=average_year, y=proportion, colour=sector_category)) +
#   geom_errorbar(aes(ymin = prop_min, ymax = prop_max, x=average_year, colour=sector_category), width = 1.5) +
#   geom_line(aes(x=average_year, y=proportion, colour=sector_category)) +
#   scale_y_continuous(limits=c(0,1))+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 9)) +
#   theme(legend.position = "bottom")+
#   labs(fill = "Sector") +
#   facet_wrap(~Method)


