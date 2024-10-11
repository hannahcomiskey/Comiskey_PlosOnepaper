# Set working directory and load function file

library(survey)
library(dplyr)
library(tidyr)
library(haven)
library(xlsx)
library(labelled)
library(tibble)

calculate_subnatSE_data <- function(filepath, myresultsfolder, strata_null = FALSE) {
  # Set Results Location
  folder <- myresultsfolder
  
  # Analysis of most recent survey 
  country_codes_DHS <- readxl::read_excel("data/country_codes_DHS.xlsx") # country codes and names
  
  data <- read_dta(file=filepath)
  year <-  min(data$v007)
  survey <- "DHS"
  population <- "AW"
  country_code <- unique(data$v000)
  year
  val_labels(data$v024)

  data$region <- as_factor(data$v024, levels = "labels")
  
  # Percentage currently using a modern method
  data <- data %>% mutate(mcpr= case_when(v313==3 ~ 1, TRUE ~ 0))
  
  # Percentage married
  data <- data %>% mutate(married= case_when(v502==1 ~ 1, TRUE ~ 0))
  
  # Religion
  val_labels(data$v130)
  data <- data %>% mutate(religion = case_when(v130==1 ~ "traditional",
                                                   v130==2 ~ "spiritual",
                                                   v130==3 ~ "christian",
                                                   v130==6 ~ "other"))
  
  # Wealth quintiles
  val_labels(data$v190)
  data <- data %>% mutate(wealth_quintiles = case_when(v190==1 ~ "poorest",
                                                   v190==2 ~ "poorer",
                                                   v190==3 ~ "middle",
                                                   v190==4 ~ "richer",
                                                   v190==5 ~ "richest"))
  # Methods
  val_labels(data$v312)
  data <- data %>% mutate(modern_method= case_when(v312==6 ~ "Sterilization (female)",
                                                   v312==7 ~ "Sterilization (male)",
                                                   v312==2 ~ "IUD",
                                                   v312==11 ~ "Implant",
                                                   v312==3 ~ "Injectable",
                                                   v312==1 ~ "OC Pills",
                                                   v312==5 ~ "Condom (m)"  ,
                                                   v312==13 ~ "LAM" ,
                                                   v312==4 | v312==14 | v312==15 | v312==16 | v312==17 | v312==18 | v312==19 | v312==20 ~ "Other Modern Methods"  ,
                                                   v312==0 | v312==8 | v312==9 | v312==10 | v312==12 ~ "None" ))
  
  data <- data %>% mutate(modern_method_source= case_when(v312==6 ~ "Sterilization (F)",
                                                          v312==7 ~ "Sterilization (M)",
                                                          v312==2 ~ "IUD",
                                                          v312==11 ~ "Implant",
                                                          v312==3 ~ "Injectable",
                                                          v312==1 ~ "OC Pills",
                                                          v312==5 ~ "Condom (M)",
                                                          v312==14 ~ "Condom (F)",
                                                          v312==4 | v312==15 | v312== 17 | v312== 18 | v312== 19 | v312== 20 ~ "Other Modern Methods",
                                                          v312==16 ~ "Emergency contraception",
                                                          v312==0 | v312==8 | v312==9 | v312==10 | v312==12 ~ "None"))
  
  
  val_labels(data$v326) # detailed breakdown of sectors
  val_labels(data$v327) # overview breakdown of sectors
  
  data <- data %>% mutate(sector = case_when(v327==1 | v327==2 ~ "Public",
                                             v327==3 ~ "NGO",
                                             v327==4 ~ "Private Clinic",
                                             v327==5 ~ "Pharmacy",
                                             v327==6 ~ "Shop Church Friend",
                                             v327==7 ~ "Other",
                                             v327==8 ~  NA_character_))
  
  data <- data %>% mutate(sector_categories= case_when(v327==1 | v327==2 ~ "Public",
                                                       v327==3 ~ "Commercial_medical",
                                                       v327==4 ~ "Commercial_medical",
                                                       v327==5 ~ "Commercial_medical",
                                                       v327==6 ~ "Other",
                                                       v327==7 ~ "Other",
                                                       v327==8 ~  NA_character_))
  
  data$num <- 1
  #sel <- select(data, modern_method_source,  region,sector, method_region, method_region_sector )
  data$sampleweights <- data$v005/1000000
  
  n_methods <- c("Implant", "Injectable", "IUD", "Condom (M)", "Sterilization (F)", "OC Pills", "Emergency contraception", "Other Modern Methods")
  data <- data %>% filter(is.na(sector_categories)==FALSE & modern_method_source %in% n_methods)
  data$method_region <- paste(data$modern_method_source, "_", data$region,  sep="")
  data$method_region_sector <- paste(data$modern_method_source, "_", data$region,"_", data$sector,  sep="")
  data$num <- 1
  
  data %>% count(sector)
  data %>% count(sector_categories)
  data %>% count(modern_method_source)
  
  ########################################################
  # Counts of observations in each sector, method and region
  ########################################################
  
  counts <- data %>% group_by(sector_categories) %>% count(method_region) %>% separate(method_region, c("modern_method_source","region"), "_")
  
  
  ########################################################
  # Calculate Percent distribution of current users of modern methods, by most recent source of method
  ########################################################
  
  test <- data %>% select(v021, v023, sampleweights, modern_method_source, sector_categories)
  
  # Complex sample design parameters
  if(strata_null==TRUE) {
    DHSdesign <- svydesign(ids=~v021, strata=NULL, weights=~sampleweights, data=data, nest=TRUE)
    options(survey.lonely.psu="adjust")
  } else {
    DHSdesign <- svydesign(ids=~v021, strata=~v023, weights=~sampleweights, data=data, nest=TRUE)
    options(survey.lonely.psu="adjust")
  }
  
  ##############################################
  # Remove any stratum that only have one PSU
  ##############################################
  
  one_psu_strata <- data %>% count(v023) %>% filter(n == 1)
  data <- data %>%
    filter(!(v023 %in% one_psu_strata$v023))
  
  ##################################
  ## Calculate SE  ---------------
  ##################################
  method.sector <- as.data.frame(prop.table(svytable(~method_region +sector, DHSdesign), 1))
  method.sector <- method.sector %>% separate(method_region, c("modern_method_source","region"), "_") %>%
    rename(Subnat_Freq=Freq)
  
  my_SEdf_public <- as_tibble(svyby(~I(sector_categories=="Public"), ~I(method_region), design=DHSdesign, svyciprop)) %>% 
    separate(`I(method_region)`, c("modern_method_source","region"), "_") %>%
    rename(Subnat_Freq=`I(sector_categories == "Public")`) %>%
    rename(Subnat_Freq_SE = `se.as.numeric(I(sector_categories == "Public"))`) %>%
    mutate(sector_categories="Public")
  
  my_SEdf_comm_med <- as_tibble(svyby(~I(sector_categories=="Commercial_medical"), ~I(method_region), design=DHSdesign, svyciprop)) %>% 
    separate(`I(method_region)`, c("modern_method_source","region"), "_") %>%
    rename(Subnat_Freq=`I(sector_categories == "Commercial_medical")`) %>%
    rename(Subnat_Freq_SE = `se.as.numeric(I(sector_categories == "Commercial_medical"))`) %>%
    mutate(sector_categories="Commercial_medical")
  
  my_SEdf_other <- as_tibble(svyby(~I(sector_categories=="Other"), ~I(method_region), design=DHSdesign, svyciprop)) %>% 
    separate(`I(method_region)`, c("modern_method_source","region"), "_") %>%
    rename(Subnat_Freq=`I(sector_categories == "Other")`) %>%
    rename(Subnat_Freq_SE = `se.as.numeric(I(sector_categories == "Other"))`) %>%
    mutate(sector_categories="Other")
  
  my_SEdf <- rbind(my_SEdf_comm_med, my_SEdf_other)
  my_SEdf <- rbind(my_SEdf, my_SEdf_public) 
  my_SEdf$country_code <- country_code
  my_SEdf$year <- year
  my_SEdf <- merge(my_SEdf, counts)
  
  write.xlsx(as.data.frame(my_SEdf), paste(folder, "/", country_code,"_", year, "_SEdf.xlsx" , sep=""))
  
  return(head(my_SEdf))
}


# List files in directory

my_files <- list.files("data/dta_files/")

dud_list <- vector()
for(i in my_files) {
  print(i)
  test_data <- read_dta(file= paste0("data/dta_files/", i))
  if(is.null(val_labels(test_data$v024))==TRUE) {
    dud_list <- c(i, dud_list)
    next
  } else {
    tryCatch({calculate_subnational_data(i,   "data/results")}, error=function(e){
      dud_list <- c(i, dud_list)
      cat("ERROR :",conditionMessage(e), "\n")
    })
  }
  saveRDS(dud_list, "data/dud_list_jobrun.RDS")
}

