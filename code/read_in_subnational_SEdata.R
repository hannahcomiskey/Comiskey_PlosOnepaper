##############################################
# Read in SE data
##############################################
subnat_FPsource_data <- readRDS("data/subnat_SE_source_data_20.RDS") # All countries

FP_2030_countries <- c("Kenya","Cameroon", 'Benin', 'India')

subnat_FPsource_data <- subnat_FPsource_data %>% 
  #filter(Country %in% FP_2030_countries) %>%
  ungroup() %>%
  dplyr::select(Country, Region, Method, sector_category, average_year, proportion, SE.proportion, n) %>%
  dplyr::mutate(Region = stringr::str_to_title(Region))

###############################################
# Updated cleaning approach -------------------
###############################################
subnat_FPsource_data <- subnat_FPsource_data %>%
  dplyr::filter(Region!="NA")

# Filter where the sample size is smaller than 2 people across all three sectors ---------

FP_source_data_wide <- subnat_FPsource_data %>% # Proportion data
  dplyr::ungroup() %>%
  dplyr::select(Country, Region, Method,  average_year, sector_category, proportion, SE.proportion, n) %>%
  tidyr::pivot_wider(names_from = sector_category, values_from = c(proportion,SE.proportion,n)) %>% # separate data into columns for each sector
  dplyr::rename(Commercial_medical = proportion_Commercial_medical,
                Public = proportion_Public,
                Other = proportion_Other,
                Public.SE = SE.proportion_Public,
                Commercial_medical.SE = SE.proportion_Commercial_medical,
                Other.SE = SE.proportion_Other) %>%
  dplyr::arrange(Country) %>%
  dplyr::filter(n_Commercial_medical>=10 | n_Public>=10 | n_Other>=10) 

# Make sure proportions add to 1 ----------
FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::rowwise() %>%
  dplyr::mutate(check_total = sum(Commercial_medical, Other, Public, na.rm = TRUE)) %>% 
  dplyr::select(Country, Region, Method, average_year, Commercial_medical, Other, Public, Commercial_medical.SE, Other.SE, Public.SE, n_Commercial_medical, n_Other, n_Public, check_total)

col_index <- which(colnames(FP_source_data_wide)=="Commercial_medical")-1 # column index before CM column, as CM=1
# When check_Total=1, replace missing values with 0 ----------
for (i in 1:nrow(FP_source_data_wide)) {
  if(FP_source_data_wide$check_total[i]>0.99) {
    na_cols <- which(is.na(FP_source_data_wide[i, c("Commercial_medical", "Other", "Public")])==TRUE) # NA values
    FP_source_data_wide[i, na_cols+col_index] <- as.list(rep(0, length(na_cols)))
  }
}

# Transform exactly 1 and 0 values away from boundary using lemon-squeezer approach ---------
FP_source_data_wide <- FP_source_data_wide %>%
  dplyr::mutate(Commercial_medical = (Commercial_medical*(nrow(FP_source_data_wide)-1)+0.33)/nrow(FP_source_data_wide)) %>%   # Y and SE transformation to account for (0,1) limits (total in sector)
  dplyr::mutate(Other = (Other*(nrow(FP_source_data_wide)-1)+0.33)/nrow(FP_source_data_wide)) %>%
  dplyr::mutate(Public = (Public*(nrow(FP_source_data_wide)-1)+0.33)/nrow(FP_source_data_wide)) %>%
  dplyr::select(Country, Region, Method, average_year, Commercial_medical, Other, Public, Commercial_medical.SE, Other.SE, Public.SE, n_Other, n_Public, n_Commercial_medical, check_total) #, count_NA, remainder)

# Clean SE values --------------
FP_source_data_wide$count_SE.NA <- rowSums(is.na(FP_source_data_wide %>% dplyr::select(Public.SE, Commercial_medical.SE, Other.SE))) # count NAs

SE_source_data_wide_norm <- FP_source_data_wide %>% dplyr::filter(count_SE.NA==0 & Commercial_medical.SE>0 & Public.SE>0) # Normal obs. No action needed.
SE_source_data_wide_X <- FP_source_data_wide %>% dplyr::filter(Commercial_medical.SE<0.001 | Public.SE<0.001) # Get obs with two missing sectors

# Replacing SE=0 and no NAs: Leontine's suggestion 
SE_source_data_wide_X <- SE_source_data_wide_X %>%
  dplyr::filter(n_Public>=20 | n_Commercial_medical >=20 | n_Other >=20) # Remove small sample sizes (DHS has 10 units sampled per cluster as min., 20 as average)
col_index <- which(colnames(SE_source_data_wide_X)=="Commercial_medical.SE")-1 # column index before CM column, as CM=1
DEFT_data <- readxl::read_xlsx("data/DEFT_DHS_database.xlsx") %>%
  rename(average_year = Year)
SE_source_data_wide_X <- SE_source_data_wide_X %>% left_join(DEFT_data)

if(nrow(SE_source_data_wide_X)>0) {
  for(i in 1:nrow(SE_source_data_wide_X)) {
    num.SE0 <- which(SE_source_data_wide_X[i,c("Commercial_medical.SE","Other.SE","Public.SE")]==0)
    num.SEnon0 <- which(SE_source_data_wide_X[i,c("Commercial_medical.SE","Other.SE","Public.SE")]!=0)
    num.SEverytiny <- which(SE_source_data_wide_X[i,c("Commercial_medical.SE","Other.SE","Public.SE")]>0 & SE_source_data_wide_X[i,c("Commercial_medical.SE","Other.SE","Public.SE")]<0.001) 
    num.SEna <- which(is.na(SE_source_data_wide_X[i,c("Commercial_medical.SE","Other.SE","Public.SE")])==TRUE)
    if(length(num.SE0)==1 & length(num.SEnon0)>0) {
      mean.SE <- mean(as.vector(unlist(SE_source_data_wide_X[i,col_index+num.SEnon0]))) # Noticed that other two columns have identical SE. Assign SE to third column.
      SE_source_data_wide_X[i,col_index+num.SE0] <- mean.SE
    } else {
      DEFT <- ifelse(is.na(SE_source_data_wide_X$DEFT[i])==TRUE, 1.5, SE_source_data_wide_X$DEFT[i])
      N1 <- length(which(SE_source_data_wide_X$Public>0.99)) + length(which(SE_source_data_wide_X$Commercial_medical>0.99)) # Number of obs=1
      phat <- (N1 + 1/2)/(nrow(FP_source_data_wide)+1)
      SE.hat <- sqrt((phat*(1-phat))/(N1+1))
      SE_source_data_wide_X[i,c(col_index+num.SEna, col_index+num.SE0, col_index+num.SEverytiny)] <- SE.hat*DEFT
    }
  }
}

FP_source_data_wide <- bind_rows(SE_source_data_wide_norm, SE_source_data_wide_X) # Put data back together again

# Remove proportions with two sectors still missing 
FP_source_data_wide <- FP_source_data_wide %>% 
  filter(is.na(Public)==FALSE & is.na(Other)==FALSE | is.na(Public)==FALSE & is.na(Commercial_medical)==FALSE | is.na(Commercial_medical)==FALSE & is.na(Other)==FALSE)

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
# FP_source_tmp <- FP_source_data_wide %>%
#   dplyr::filter(Country=="Rwanda") %>%
#   dplyr::mutate(Region = case_when(Region == "Ouest" ~ "West",
#                                    Region == "Nord" ~ "North",
#                                    Region == "Sud" ~ "South",
#                                    Region == "Est" ~ "East",
#                                    Region == "Ville de Kigali" ~ "Kigali",
#                                    Region == "Ville De Kigali" ~ "Kigali",
#                                    Region == "Kigali City" ~ "Kigali",
#                                    Region == "Butare, Gitarama (Central, South)" ~ "South",
#                                    Region == "Cyangugu, Gikongoro (Southwest)" ~ "South",
#                                    Region == "Byumba, Kibungo, Umutara (Northeast)" ~ "North",
#                                    Region == "Gisenyi, Kibuye, Ruhengeri (Northwest)" ~ "West",
#                                    TRUE ~ as.character(Region))) 
# 
# FP_source_data_wide <- FP_source_data_wide %>%
#   dplyr::filter(Country!="Rwanda")
# 
# FP_source_data_wide <- FP_source_data_wide %>%
#   merge(FP_source_tmp, all = TRUE)
# 
# # Replace issues with Nigeria names
# FP_source_tmp <- FP_source_data_wide %>%
#   dplyr::filter(Country=="Nigeria") %>%
#   dplyr::mutate(Region = case_when(Region == "Northeast" ~ "North East",
#                                    Region == "Northwest" ~ "North West",
#                                    Region == "Southeast" ~ "South East",
#                                    Region == "Southwest" ~ "South West",
#                                    TRUE ~ as.character(Region)))
# 
# FP_source_data_wide <- FP_source_data_wide %>%
#   dplyr::filter(Country!="Nigeria")
# 
# # Replace issues with Cote d'Ivoire names
# FP_source_tmp <- FP_source_data_wide %>%
#   dplyr::filter(Country=="Cote d'Ivoire") %>%
#   dplyr::mutate(Region = case_when(Region == "Center-East" ~ "Center East",
#                                    Region == "Center-North" ~ "Center North",
#                                    Region == "Center-West" ~ "Center West",
#                                    Region == "Center-South" ~ "Center South",
#                                    TRUE ~ as.character(Region)))
# 
# FP_source_data_wide <- FP_source_data_wide %>%
#   dplyr::filter(Country!="Cote d'Ivoire")
# 
# FP_source_data_wide <- FP_source_data_wide %>%
#   merge(FP_source_tmp, all = TRUE)

#################################################
# Plot raw data -------------------------------
#################################################

# tmp1 <- FP_source_data_wide %>% dplyr::filter(Country=="Congo Democratic Republic" & Region =="Katanga") %>%
#   select(Country, Method, Region, average_year, Other, Public, Commercial_medical) %>%
#   pivot_longer(cols =c(Other, Public, Commercial_medical), names_to = "sector_category", values_to = "proportion")
# 
# tmp2 <- FP_source_data_wide %>% dplyr::filter(Country=="Congo Democratic Republic" & Region =="Katanga") %>%
#   select(Country, Method, Region, average_year, Other.SE, Public.SE, Commercial_medical.SE) %>%
#   pivot_longer(cols =c(Other.SE, Public.SE, Commercial_medical.SE), names_to = "sector_category", values_to = "SEproportion") %>%
#   mutate(sector_category = substr(sector_category, 1, nchar(sector_category)-3))
# 
# plot_data <- merge(tmp1, tmp2) %>%
#   mutate(prop_min = ifelse(proportion - 2*SEproportion<0, 0, proportion - 2*SEproportion)) %>%
#   mutate(prop_max = ifelse(proportion + 2*SEproportion>1, 1, proportion + 2*SEproportion))
#   
# ggplot(data = plot_data) +
#   labs(y="Proportion of contraceptives supplied", x = "Year", title = paste0(unique(tmp2$Region),",",unique(tmp2$Country))) +
#   geom_point(aes(x=average_year, y=proportion, colour=sector_category)) +
#   geom_errorbar(data=plot_data, aes(ymin = prop_min, ymax = prop_max, x=average_year, colour=sector_category), width = 1.5) +
#   geom_line(aes(x=average_year, y=proportion, colour=sector_category)) +
#   scale_y_continuous(limits=c(0,1))+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 9)) +
#   theme(legend.position = "bottom")+
#   labs(fill = "Sector") +
#   facet_wrap(~Method)