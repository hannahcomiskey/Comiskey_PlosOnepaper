library(survey)
library(dplyr)
library(tidyr)
library(haven)
library(labelled)
library(tibble)
library(ipumsr)

########################################################
# Analysis of IPUMS data 
########################################################
# Read data in
ddi <- read_ipums_ddi("data/idhs_00013.xml")
ipums_data <- read_ipums_micro(ddi)

# Convert the labels to factors (and drop the unused levels)
ipums_data <- ipums_data %>%
  mutate(Country = as_factor(lbl_clean(COUNTRY)))

table(ipums_data$Country, useNA = "always") %>% nrow()

# Count surveys for each country - influences choice of region column
n_surveys <- ipums_data %>%
  group_by(Country) %>%
  select(Country, YEAR) %>%
  distinct() %>%
  count()

country_codes <- readxl::read_excel("data/country_codes_DHS.xlsx") %>%
  select(Code, `Country Name`) %>%
  rename(Country = `Country Name`)

surveys_country_codes <- merge(n_surveys, country_codes) 

##############################################
# Replace numeric codes with district names
##############################################

refined_data <- ipums_data %>%
  mutate(Afghanistan_subnat = as_factor(lbl_clean(GEO_AF2015))) %>%
  mutate(Benin_subnat = as_factor(lbl_clean(GEO_BJ1996_2017))) %>%
  mutate(BFaso_subnat_93 = as_factor(lbl_clean(GEO_BF1993))) %>%
  mutate(BFaso_subnat_98 = as_factor(lbl_clean(GEO_BF1998))) %>%
  mutate(BFaso_subnat = as_factor(lbl_clean(GEO_BF2003_2018))) %>%
  mutate(Cameroon_subnat = as_factor(lbl_clean(GEO_CM1991_2022))) %>%
  mutate(DRCongo_subnat = as_factor(lbl_clean(GEO_CD2007_2013))) %>%
  mutate(CdIvoire_subnat_98 = as_factor(lbl_clean(GEO_CI1998))) %>%
  mutate(CdIvoire_subnat_11 = as_factor(lbl_clean(GEO_CI2011))) %>%
  mutate(CdIvoire_subnat_94 = as_factor(lbl_clean(GEO_CI1994))) %>%
  mutate(Ethiopia_subnat = as_factor(lbl_clean(GEO_ET2000_2019))) %>%
  mutate(Ghana_subnat = as_factor(lbl_clean(GEO_GH1988_2019))) %>%
  mutate(Guinea_subnat = as_factor(lbl_clean(GEO_GN1999_2021))) %>%
  mutate(India_subnat = as_factor(lbl_clean(GEO_IA1992_2019))) %>%
  mutate(Kenya_subnat = as_factor(lbl_clean(GEO_KE1989_2014))) %>%
  mutate(Liberia_subnat = as_factor(lbl_clean(GEO_LR2007_2019))) %>%
  mutate(Madagascar_subnat = as_factor(lbl_clean(GEO_MG1992_2021))) %>%
  mutate(Malawi_subnat = as_factor(lbl_clean(GEO_MW1992_2017))) %>%
  mutate(Mali_subnat = as_factor(lbl_clean(GEO_ML1987_2021))) %>%
  mutate(Mozambique_subnat = as_factor(lbl_clean(GEO_MZ1997_2018))) %>%
  mutate(Myanmar_subnat = as_factor(lbl_clean(GEO_MM2015))) %>%
  mutate(Nepal_subnat = as_factor(lbl_clean(GEO_NP1996_2016))) %>%
  mutate(Niger_subnat = as_factor(lbl_clean(GEO_NE1992_2021))) %>%
  mutate(Nigeria_subnat_99 = as_factor(lbl_clean(GEO_NG1999))) %>%
  mutate(Nigeria_subnat_90 = as_factor(lbl_clean(GEO_NG1990))) %>%
  mutate(Nigeria_subnat = as_factor(lbl_clean(GEO_NG2003_2021))) %>%
  mutate(Pakistan_subnat = as_factor(lbl_clean(GEO_PK1991_2017))) %>%
  mutate(Rwanda_subnat = as_factor(lbl_clean(GEO_RW1992_2005))) %>%
  mutate(Rwanda_subnat_0519 = as_factor(lbl_clean(GEO_RW2005_2019))) %>%
  mutate(Senegal_subnat = as_factor(lbl_clean(GEO_SN1986_2021))) %>%
  mutate(Tanzania_subnat = as_factor(lbl_clean(GEO_TZ1991_2017))) %>%
  mutate(Uganda_subnat = as_factor(lbl_clean(GEO_UG1995_2018))) %>%
  mutate(Zimbabwe_subnat = as_factor(lbl_clean(GEO_ZW1994_2015))) 

##############################################
# Replace source columns with labels
##############################################

refined_data <- refined_data %>%
  mutate(Current_source_standard = as_factor(lbl_clean(FPLASTSRCS))) %>%
  mutate(Current_source_detailed = as_factor(lbl_clean(FPLASTSRCD))) %>%
  mutate(Religion = as_factor(lbl_clean(RELIGION))) %>%
  mutate(Wealth_quintiles = as_factor(lbl_clean(WEALTHQ))) %>%
  mutate(Method = as_factor(lbl_clean(FPMETHNOW))) %>%
  select(PSU, DOMAIN, SAMPLE, Country, YEAR, PERWEIGHT, Current_source_standard, Current_source_detailed, Method, Afghanistan_subnat:Zimbabwe_subnat)

##############################################
# Reduce to 5 methods
##############################################

subnat_source_data <- refined_data %>% 
  filter(Method %in% c("Pill", "Norplant/Implants", "Female Sterilization", "Injections", "IUD"))

##############################################
# Clean up 
##############################################
# Mapping detailed columns to Current_source labels
Current_source_map <- subnat_source_data %>%
  ungroup() %>%
  select(Current_source_standard, Current_source_detailed) %>%
  distinct()

subnat_source_data <- subnat_source_data %>% 
  mutate(Current_source_standard = as.character(Current_source_standard)) %>% 
  mutate(Current_source_standard = 
           ifelse(is.na(Current_source_standard),
                  purrr::map_chr(.x = .$Current_source_detailed, ~ as.character(Current_source_map$Current_source_standard[match(.x, Current_source_map$Current_source_detailed)])), 
                  Current_source_standard))

subnat_source_data <- subnat_source_data %>% mutate(sector_standard = case_when(Current_source_standard=="Govt Clinic/Pharm" | Current_source_standard=="Govt Home/Comm delivery" ~ "Public",
                                                                                Current_source_standard=="NGO" ~ "Private",
                                                                                Current_source_standard=="Private Clin/Deliv" ~ "Private",
                                                                                Current_source_standard=="Private Pharmacy" ~ "Private",
                                                                                Current_source_standard=="Church, Shop, friends, books" ~ "Private",
                                                                                Current_source_standard=="Other" ~ "Private",
                                                                                Current_source_standard=="Don't know" | Current_source_standard=="Missing" | Current_source_standard=="NIU (not in universe)"  ~  NA_character_))


#saveRDS(subnat_source_data, file = "data/IPUMS/draft_subnat_source_data_1_ipums.RDS")

# Address issue with Rwanda subnational multi-names 
subnat_source_data <- subnat_source_data %>%
  mutate(Rwanda_subnat = case_when(Rwanda_subnat== "Butare, Gitarama (Central, South)" ~ "South",
                                   Rwanda_subnat== "Cyangugu, Gikongoro (Southwest)" ~ "South",
                                   Rwanda_subnat== "Byumba, Kibungo, Umutara (Northeast)" ~ "North",
                                   Rwanda_subnat== "Gisenyi, Kibuye, Ruhengeri (Northwest)" ~ "West",
                                   TRUE ~ as.character(Rwanda_subnat))) 

# Collapse regions into one column 
Region <- vector()
colindex1 <- which(colnames(subnat_source_data)=="Afghanistan_subnat")
colindex2 <- which(colnames(subnat_source_data)=="Zimbabwe_subnat")
for(i in 1:nrow(subnat_source_data)) {
  reg_index <- which(!is.na(subnat_source_data[i,colindex1:colindex2]))+(colindex1-1)
  if(length(reg_index)==0) {
    Region[i] <- NA
  } else {
    Region[i] <- as.vector(unlist(subnat_source_data[i,reg_index]))
  }
}

subnat_source_data_2 <- subnat_source_data
subnat_source_data_2$Region <- as.factor(Region)
subnat_source_data_2$sector <- as.factor(subnat_source_data_2$sector_standard)

saveRDS(subnat_source_data_2, file="data/IPUMS/draft_subnat_source_data_2_ipums.RDS")

###########################
# Select required columns
###########################
subnat_source_data_3 <- subnat_source_data_2 %>%
  select(!Afghanistan_subnat:Zimbabwe_subnat)

# Append columns together for pivots
subnat_source_data_3$region_country_year <- as.factor(paste(subnat_source_data_3$Region, "_", subnat_source_data_3$Country, "_", subnat_source_data_3$YEAR, sep=""))
subnat_source_data_3$method_region_country_year <- as.factor(paste(subnat_source_data_3$Method, "_", subnat_source_data_3$Region, "_", subnat_source_data_3$Country, "_", subnat_source_data_3$YEAR, sep=""))
subnat_source_data_3$method_region_country_sector <- as.factor(paste(subnat_source_data_3$Method, "_", subnat_source_data_3$Region,"_", subnat_source_data_3$Country,"_", subnat_source_data_3$sector,sep=""))
subnat_source_data_3$num <- 1

##############################################
# Remove any stratum that only have one PSU
##############################################
one_psu_strata <- subnat_source_data_3 %>% count(DOMAIN) %>% filter(n == 1)
subnat_source_data_3 <- subnat_source_data_3 %>%
  filter(!(DOMAIN %in% one_psu_strata$DOMAIN))

saveRDS(subnat_source_data_3, file="data/IPUMS/draft_subnat_source_data_3_ipums.RDS")

#############################################################
# Set up survey design - filter out NAs from design columns
#############################################################
n_methods <- c("Pill", "Norplant/Implants", "Female Sterilization", "Injections", "IUD")

subnat_source_data_3 <- subnat_source_data_3 %>% 
  filter(is.na(PSU)==FALSE & is.na(method_region_country_year)==FALSE & is.na(sector)==FALSE & is.na(Method)==FALSE & is.na(DOMAIN)==FALSE & is.na(PERWEIGHT)==FALSE & is.na(YEAR)==FALSE)

IPUMSdesign<-svydesign(id= subnat_source_data_3$PSU, strata=subnat_source_data_3$DOMAIN, weights=subnat_source_data_3$PERWEIGHT, data=subnat_source_data_3, nest=TRUE)
options(survey.lonely.psu="adjust")

##################################
## Calculate SE  ---------------
##################################

folder <- "data/IPUMS/SEdf"

method.sector <- as.data.frame(prop.table(svytable(~method_region_country_year +sector, IPUMSdesign), 1))
method.sector <- method.sector %>% separate(method_region_country_year, c("Method","Region", "Country", "Year"), "_") %>%
  rename(Subnat_Freq=Freq)

# Subset design for sectors
country_year_comb <- subnat_source_data_3 %>%
  group_by(Country) %>%
  select(Country, YEAR) %>%
  distinct()

rwanda_combs <- country_year_comb %>% filter(Country=="Rwanda")

for(i in 1:nrow(country_year_comb)) {
  country_code <- country_year_comb$Country[i]
  year <- country_year_comb$YEAR[i]
  tmp <- subnat_source_data_3 %>% filter(Country==country_code & YEAR==year)
  tmp$method_region <- as.factor(paste(tmp$Method, "_", tmp$Region, sep=""))
  tmpdesign<-svydesign(id= tmp$PSU, strata=tmp$DOMAIN, weights=tmp$PERWEIGHT, data=tmp, nest=TRUE)
  
  counts <- tmp %>% 
    group_by(sector) %>% 
    count(method_region) %>% 
    separate(method_region, c("Method","Region"), "_") %>%
    rename(sector_categories = sector)
  
  my_SEdf_tmp1 <- as_tibble(svyby(~I(sector=="Commercial_medical"), ~I(method_region), design=tmpdesign, svyciprop)) %>%
    separate(`I(method_region)`, c("Method","Region"), "_") %>%
    mutate(sector_categories="Commercial_medical") %>%
    mutate(Country = country_code) %>%
    mutate(Year = year) %>%
    rename(Subnat_Freq=`I(sector == \"Commercial_medical\")`) %>%
    rename(Subnat_Freq_SE = `se.as.numeric(I(sector == \"Commercial_medical\"))`)
  
  my_SEdf_tmp2 <- as_tibble(svyby(~I(sector=="Public"), ~I(method_region), design=tmpdesign, svyciprop)) %>%
    separate(`I(method_region)`, c("Method","Region"), "_") %>%
    mutate(sector_categories="Public") %>%
    mutate(Country = country_code) %>%
    mutate(Year = year) %>%
    rename(Subnat_Freq=`I(sector == \"Public\")`) %>%
    rename(Subnat_Freq_SE = `se.as.numeric(I(sector == \"Public\"))`)
  
  my_SEdf_tmp3 <- as_tibble(svyby(~I(sector=="Other"), ~I(method_region), design=tmpdesign, svyciprop)) %>%
    separate(`I(method_region)`, c("Method","Region"), "_") %>%
    mutate(sector_categories="Other") %>%
    mutate(Country = country_code) %>%
    mutate(Year = year) %>%
    rename(Subnat_Freq=`I(sector == \"Other\")`) %>%
    rename(Subnat_Freq_SE = `se.as.numeric(I(sector == \"Other\"))`)
  
  my_SEdfall <- rbind(my_SEdf_tmp3, my_SEdf_tmp2)
  my_SEdfall <- rbind(my_SEdfall, my_SEdf_tmp1) 
  my_SEdfall <- left_join(my_SEdfall, counts) 
  
  
  saveRDS(as.data.frame(my_SEdfall), paste(folder, "/",country_code ,"_", year, "_SEdf_ipums.RDS" , sep=""))
}
