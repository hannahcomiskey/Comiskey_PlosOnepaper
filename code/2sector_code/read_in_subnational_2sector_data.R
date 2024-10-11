##############################################
# Read in SE data
##############################################
subnat_FPsource_data <- readRDS("data/subnat_bivar_SE_source_data_20.RDS") # All countries

FP_2030_countries <- c("Kenya","Cameroon", 'Benin', 'India')

subnat_FPsource_data <- subnat_FPsource_data %>% 
  #filter(Country %in% FP_2030_countries) %>%
  ungroup() %>%
  dplyr::select(Country, Region, Method, average_year, Public, Private, Public_n, Private_n) %>%
  dplyr::mutate(Region = stringr::str_to_title(Region))

###############################################
# Updated cleaning approach -------------------
###############################################
subnat_FPsource_data <- subnat_FPsource_data %>%
  dplyr::filter(Region!="NA")

# Filter where the sample size is smaller than 2 people across all three sectors ---------

FP_source_data_wide <- subnat_FPsource_data %>% # Proportion data
  dplyr::ungroup() %>%
  dplyr::select(Country, Region, Method,  average_year, Public, Private, Public_n, Private_n) %>%
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
  dplyr::select(Country, Region, Method, average_year, Public, Private, Public_n, Private_n, check_total) #, count_NA, remainder)

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