# # Calculation for sector data
P_median_global_2files <- function(subnat_index_table,
                                   method_index_table,
                                   sector_type,
                                   year_index_table,
                                   file_chain1,
                                   file_chain2) {
  
  mod_P1 <- readRDS(file_chain1)
  mod_P2 <- readRDS(file_chain2)
  
  P_dims <-  dim(mod_P1)
  
  time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
    dplyr::filter(average_year>floored_year) %>%
    dplyr::select(index_year) %>%
    unlist() %>%
    as.vector()
  
  averageyear_index_table <- year_index_table %>%
    dplyr::filter(average_year>floored_year) %>%
    select(average_year, index_year) %>%
    dplyr::mutate(index_year = 1:n())
  
  P_s_med <- array(dim=c(P_dims[3], length(time_index), P_dims[2])) # subnat, time, method
  
  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_s_med[j,k,r] <- stats::median(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]))
      }
    }
  }
  
  # P_dims <- dim(Psamps)
  #
  # #### P median
  # P_s_med <- array(dim=c(length(time_index),P_dims[2],P_dims[3])) # method, year, subnat
  #
  # # Create a table for storing individual true country public data
  # for(k in 1:length(time_index)) { # time loop
  #   for(s in 1:P_dims[3]) { # subnat
  #     for (m in 1:P_dims[2]) { # method
  #       P_s_med[k,m,s] <- stats::median(Psamps[,m,s,time_index[k]])
  #     }
  #   }
  # }
  
  P_s_med <- plyr::adply(P_s_med, c(1,2,3))
  colnames(P_s_med) <- c("index_subnat", "index_year", "index_method",  "median_p")
  
  P_s_med <- P_s_med %>%
    group_by(index_method, index_subnat, index_year) %>%
    # summarise(median_p = stats::median(estimate)) %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(averageyear_index_table) %>%
    dplyr::mutate(Sector = sector_type)
  
  P_Q <- array(dim=c(P_dims[3], length(time_index), P_dims[2], 4)) # subnat, time, method, stats::quantile(95, 80)
  
  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_Q[j,k,r, 1:2] <- as.vector(unlist(stats::quantile(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]), prob=c(0.025, 0.975))))
        P_Q[j,k,r, 3:4] <- as.vector(unlist(stats::quantile(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]), prob=c(0.1, 0.9))))
      }
    }
  }
  
  P_Q <- plyr::adply(P_Q, c(1,2,3))
  colnames(P_Q) <- c("index_subnat", "index_year", "index_method",  "lower_95", "uppper_95", "lower_80", "upper_80")
  
  P_s_med <- P_Q %>%
    group_by(index_method, index_subnat, index_year) %>%
    # summarise(median_p = stats::median(estimate)) %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    left_join(P_s_med)
  
  return(P_s_med)
}


calculate_bivariateSE_data <- function(filepath, myresultsfolder) {

  ##############################################
  # Set Results Location
  ##############################################
  
  folder <- myresultsfolder
  
  ########################################################
  # Analysis of most recent survey 
  ########################################################
  data <- read_dta(file=filepath)
  year <- min(data$v007)
  survey <- "DHS"
  population <- "AW"
  country_code <- unique(data$v000)
  
  data$region <- as_factor(data$v024, levels = "labels")
  
  data <- data %>% mutate(mcpr= case_when(v313==3 ~ 1, TRUE ~ 0))
  
  data <- data %>% mutate(married= case_when(v502==1 ~ 1, TRUE ~ 0))
  
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
  
  
  data <- data %>% mutate(sector= case_when(v327==1 | v327==2 ~ "Public",
                                            v327==3 ~ "NGO",
                                            v327==4 ~ "Private Clinic",
                                            v327==5 ~ "Pharmacy",
                                            v327==6 ~ "Shop Church Friend",
                                            v327==7 ~ "Other",
                                            v327==8 ~  NA_character_))
  
  
  # data <- data %>% mutate(sector_categories= case_when(v327==1 | v327==2 ~ "Public",
  #                                                      v327==3 ~ "Commercial_medical",
  #                                                      v327==4 ~ "Commercial_medical",
  #                                                      v327==5 ~ "Commercial_medical",
  #                                                      v327==6 ~ "Other",
  #                                                      v327==7 ~ "Other",
  #                                                      v327==8 ~  NA_character_))
  
  data <- data %>% mutate(sector_categories= case_when(v327==1 | v327==2 ~ "Public",
                                                       v327==3 ~ "Private",
                                                       v327==4 ~ "Private",
                                                       v327==5 ~ "Private",
                                                       v327==6 ~ "Private",
                                                       v327==7 ~ "Private",
                                                       v327==8 ~  NA_character_))
  
  data$num <- 1
  data$sampleweights <- data$v005/1000000
  
  data <- data %>% filter(is.na(sector_categories)==FALSE & modern_method_source %in% c("Implant", "Injectable", "OC Pills", "Sterilization (F)", "Sterilization (M)", "Condom (M)", "Condom (F)", "IUD","Other Modern Methods", "Emergency contraception"))
  
  counts <- data %>% group_by(sector_categories) %>% count(modern_method_source) %>% rename(method = modern_method_source)
  
  #################################
  # Calculate SE  ---------------
  #################################
  DHSdesign<-svydesign(id=data$v021, strata=data$v023, weights=data$sampleweights, data=data, nest=TRUE)
  options(survey.lonely.psu="adjust")
  
  summary(DHSdesign)
  
  # Get proportions and var-covar
  d.s <- update(DHSdesign, modern_method_source = factor(modern_method_source)) 
  d.s <- update(d.s, sector_categories = factor(sector_categories)) 
  
  # Pivot counts for incorporation into data
  counts_wide <- counts %>% tidyr::pivot_wider(names_from='sector_categories', values_from = 'n')
  colnames(counts_wide)[-1] <- paste0(colnames(counts_wide)[-1], '_n')
  
  # Calculate  
  prop_mat <- svyby(~I(sector_categories), ~I(modern_method_source), design=d.s, svymean, covmat=TRUE)
  vcov_matrix <- vcov(prop_mat)
  
  # Clean up proportion matrix
  colnames(prop_mat) <- gsub("I\\(sector_categories\\)", "", colnames(prop_mat))
  colnames(prop_mat) <- gsub("I\\(modern_method_source\\)", "method", colnames(prop_mat))
  prop_mat$country_code <- country_code
  prop_mat$year <- year
  prop_mat <- merge(prop_mat, counts_wide)
  xlsx::write.xlsx(prop_mat, paste(folder, "/proportions/prop_", country_code,"_", year, "_SEdf.xlsx" , sep=""))
  
  # # Clean up covariance matrix
  colnames(vcov_matrix) <- gsub("I\\(sector_categories\\)", "", colnames(vcov_matrix))
  rownames(vcov_matrix) <- gsub("I\\(sector_categories\\)", "", rownames(vcov_matrix))
  vcov_matrix <- as_tibble(vcov_matrix) %>% mutate(Method_sector = rownames(vcov_matrix))
  vcov_matrix$country_code <- country_code
  vcov_matrix$year <- year
  xlsx::write.xlsx(vcov_matrix, paste(folder, "/varcov/varcov_",  country_code,"_", year, "_SEdf.xlsx" , sep=""))
  
  return(head(prop_mat))
}

# Calculate subnational data from DHS microdata --------------------------------
calculate_subnational_data <- function(filepath, myresultsfolder) {
  # Subnational Inputs for SStoEMU
  # Kristin Bietsch, PhD
  # Avenir Health
  # 08/08/19
  
  ##############################################
  # Set Results Location
  folder <- myresultsfolder
  
  ########################################################
  # Analysis of most recent survey 
  ########################################################
  data <- read_dta(file=filepath)
  year <- min(data$v007)
  survey <- "DHS"
  population <- "AW"
  country_code <- unique(data$v000)
  
  val_labels(data$v024)
  
  data$region <- as_factor(data$v024, levels = "labels")
  
  data <- data %>% mutate(mcpr= case_when(v313==3 ~ 1, TRUE ~ 0))
  
  data <- data %>% mutate(married= case_when(v502==1 ~ 1, TRUE ~ 0))
  
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
  
  
  #val_labels(data$v327)
  data <- data %>% mutate(sector= case_when(v327==1 | v327==2 ~ "Public",
                                            v327==3 ~ "NGO",
                                            v327==4 ~ "Private Clinic",
                                            v327==5 ~ "Pharmacy",
                                            v327==6 ~ "Shop Church Friend",
                                            v327==7 ~ "Other",
                                            v327==8 ~  NA_character_))
  
  
  data$method_region <- paste(data$modern_method_source, "_", data$region,  sep="")
  data$method_region_sector <- paste(data$modern_method_source, "_", data$region,"_", data$sector,  sep="")
  data$num <- 1
  #sel <- select(data, modern_method_source,  region,sector, method_region, method_region_sector )
  
  ########################################################
  data$sampleweights <- data$v005/1000000
  
  design <- svydesign(ids=~v021, strata=~v025, weights=~sampleweights, data=data)
  options(survey.lonely.psu="adjust")
  
  married <- data %>% filter(married==1)
  design.married <- svydesign(ids=~v021, strata=~v025, weights=~sampleweights, data=married)
  options(survey.lonely.psu="adjust")
  
  modern <- data %>% filter(modern_method!="None")
  design.modern <- svydesign(ids=~v021, strata=~v025, weights=~sampleweights, data=modern)
  options(survey.lonely.psu="adjust")
  
  
  ########################################################
  regional.distribution <- as.data.frame(svymean(~region,
                                                 design=design,
                                                 weights=~sampleweights,
                                                 data=data, drop.empty.groups=FALSE))
  regional.distribution$region <- substring(row.names(regional.distribution),7)
  regional.distribution <- regional.distribution %>% select(region, mean)
   
  #############################################################################
  # method*region interaction
  # unweighted count
  # distribution of sector
  method.sector <-as.data.frame(prop.table(svytable(~method_region +sector, design), 1))
  method.sector <- method.sector %>% separate(method_region, c("modern_method_source","region"), "_") %>%
    rename(Subnat_Freq=Freq)
  
  method_region <- svyby(~num, by=~method_region,
                         design=design,
                         weights=~sampleweights,
                         unwtd.count,
                         data=data, drop.empty.groups=FALSE)
  method_region <- method_region %>% separate(method_region, c("modern_method_source","region"), "_") %>%
    rename(Subnat_count=counts) %>% select(-se)
  method.sector <- left_join(method.sector, method_region, by=c("modern_method_source", "region"))
  
  
  method.sector.nat <- as.data.frame(prop.table(svytable(~modern_method_source +sector, design), 1))
  method.sector.nat <- method.sector.nat %>% rename(Nat_Freq=Freq)
  
  method.sector <- left_join(method.sector, method.sector.nat, by=c("modern_method_source", "sector"))
  method.sector <- method.sector %>% mutate(Frequency= case_when(Subnat_count>=25 ~ Subnat_Freq,
                                                                 Subnat_count<25 ~ Nat_Freq ))
  
  method.sector$Source <- paste(year , survey, sep=" ")
  method.sector$Population <- population
  method.sector$Country.Method <- paste(method.sector$region , ":", method.sector$modern_method_source, sep="")
  
  method.sector <- method.sector %>% select(region, Source, Population, modern_method_source, Country.Method, sector, Frequency, Subnat_count)
  method.sector <- method.sector %>% spread(sector, Frequency) 
  
  method.sector <- method.sector %>% mutate(Notes= case_when(Subnat_count>=25 ~ "",
                                                             Subnat_count<25 ~ "Used National Average because samples too small @ regional level" ))
  
  # Add in NGO column if not already present
  cols <- c(NGO= NA_real_)
  method.sector <- add_column(method.sector, !!!cols[setdiff(names(cols), names(method.sector))])
  
  mysectorlist <- c("Public", "NGO", "Private Clinic", "Pharmacy"  ,  "Shop Church Friend",  "Other" )
  method.sector <- method.sector %>% select(region, Source, Population, modern_method_source, Country.Method, colnames(method.sector)[which(colnames(method.sector) %in% mysectorlist)], Notes, Subnat_count) %>% rename(N=Subnat_count)
  
  ###############################################################
  ###############################################################
  ###############################################################
  write.xlsx(as.data.frame(method.sector), paste(folder, "/", country_code,"_", year, "_SubnationalEMUPrep.xlsx" , sep=""), sheetName="Source Data", col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)
  write.xlsx(as.data.frame(regional.distribution), paste(folder, "/", country_code,"_", year, "_SubnationalEMUPrep.xlsx" , sep=""), sheetName="Regional Distribution", col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)

  return(head(method.sector))
}


# Calculate India subnational data from DHS microdata --------------------
calculate_subnational_india_data <- function(filepath) {
  # Subnational Inputs for SStoEMU
  # Kristin Bietsch, PhD
  # Avenir Health
  # 08/08/19
  
  ##############################################
  # Set Results Location
  folder <- "data/subnational/results/india/"
  
  ########################################################
  # Analysis of most recent survey 
  ########################################################
  data <- read_dta(file=filepath) 
  country_codes_DHS <- read_excel("data/subnational/country_codes_DHS.xlsx")
  
  year <- min(data$v007)
  survey <- "DHS"
  population <- "AW"
  country_code <- unique(data$v000)
  
  region_code <- substr(filepath, 30,31)
  region_name <- country_codes_DHS %>% filter(`India States`==region_code) %>% select(`State Name`)
  region_name <- as.vector(unlist(region_name))
  
  # use file name if region is unavailable
  data$region <- ifelse(is.null(val_labels(data$v024))==TRUE, region_name, as_factor(data$v024, levels = "labels"))
  
  data <- data %>% mutate(mcpr= case_when(v313==3 ~ 1, TRUE ~ 0))
  
  data <- data %>% mutate(married= case_when(v502==1 ~ 1, TRUE ~ 0))
  
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
  
  
  #val_labels(data$v327)
  data <- data %>% mutate(sector= case_when(v327==1 | v327==2 ~ "Public",
                                            v327==3 ~ "NGO",
                                            v327==4 ~ "Private Clinic",
                                            v327==5 ~ "Pharmacy",
                                            v327==6 ~ "Shop Church Friend",
                                            v327==7 ~ "Other",
                                            v327==8 ~  NA_character_))
  
  
  data$method_region <- paste(data$modern_method_source, "_", data$region,  sep="")
  data$method_region_sector <- paste(data$modern_method_source, "_", data$region,"_", data$sector,  sep="")
  data$num <- 1
  #sel <- select(data, modern_method_source,  region,sector, method_region, method_region_sector )
  
  ########################################################
  data$sampleweights <- data$v005/1000000
  
  design <- svydesign(ids=~v021, strata=~v025, weights=~sampleweights, data=data)
  options(survey.lonely.psu="adjust")
  
  married <- data %>% filter(married==1)
  design.married <- svydesign(ids=~v021, strata=~v025, weights=~sampleweights, data=married)
  options(survey.lonely.psu="adjust")
  
  modern <- data %>% filter(modern_method!="None")
  design.modern <- svydesign(ids=~v021, strata=~v025, weights=~sampleweights, data=modern)
  options(survey.lonely.psu="adjust")
  
  
  ########################################################
  #regional.distribution <- as.data.frame(svymean(~region,
  #design=design,
  #weights=~sampleweights,
  #data=data, drop.empty.groups=FALSE))
  #regional.distribution$region <- substring(row.names(regional.distribution),7)
  #regional.distribution <- regional.distribution %>% select(region, mean)
  
  ########################################################
  mcpr.married <-  svyby(~mcpr, by=~region,
                         design=design.married,
                         FUN=svymean,
                         na.rm=TRUE,
                         weights=~sampleweights,
                         data=married, drop.empty.groups=FALSE)
  mcpr.married <- mcpr.married %>% rename(mCPR_MW=mcpr, Country=region) %>% select(-se)
  
  mcpr.all <-  svyby(~mcpr, by=~region,
                     design=design,
                     FUN=svymean,
                     na.rm=TRUE,
                     weights=~sampleweights,
                     data=data, drop.empty.groups=FALSE)
  mcpr.all <- mcpr.all %>% rename(mCPR_AW=mcpr, Country=region) %>% select(-se)
  
  mcpr <- full_join(mcpr.married, mcpr.all, by="Country")
  mcpr$Year <- year
  mcpr$CountryYear <- paste(mcpr$Country, mcpr$Year , sep="")
  mcpr$Survey <- paste(year , survey, sep=" ")
  mcpr <- mcpr %>% select(Country, Year, CountryYear, Survey, mCPR_MW, mCPR_AW)
  
  #############################################################################
  # if want to check modern method coding
  #svymean(~modern_method, design=design,weights=~sampleweights,data=data, drop.empty.groups=FALSE)
  
  method.prev <-as.data.frame(prop.table(svytable(~modern_method+region, design = design.modern), margin = 2))
  method.prev <- method.prev %>% spread(modern_method, Freq) 
  method.prev$Source <- paste(year , survey, sep=" ")
  method.prev$Population <- population
  
  mymethodslist <- c("Sterilization (female)", "Sterilization (male)", "IUD", "Implant", "Injectable", "OC Pills" , "Condom (m)", "LAM", "Other Modern Methods" )
  
  method.prev <- method.prev %>% select(region, Source, Population, colnames(method.prev)[which(colnames(method.prev) %in% mymethodslist)])
  
  #############################################################################
  # method*region interaction
  # unweighted count
  # distribution of sector
  method.sector <-as.data.frame(prop.table(svytable(~method_region +sector, design), 1))
  method.sector <- method.sector %>% separate(method_region, c("modern_method_source","region"), "_") %>%
    rename(Subnat_Freq=Freq)
  
  method_region <- svyby(~num, by=~method_region,
                         design=design,
                         weights=~sampleweights,
                         unwtd.count,
                         data=data, drop.empty.groups=FALSE)
  method_region <- method_region %>% separate(method_region, c("modern_method_source","region"), "_") %>%
    rename(Subnat_count=counts) %>% select(-se)
  method.sector <- left_join(method.sector, method_region, by=c("modern_method_source", "region"))
  
  
  method.sector.nat <- as.data.frame(prop.table(svytable(~modern_method_source +sector, design), 1))
  method.sector.nat <- method.sector.nat %>% rename(Nat_Freq=Freq)
  
  method.sector <- left_join(method.sector, method.sector.nat, by=c("modern_method_source", "sector"))
  method.sector <- method.sector %>% mutate(Frequency= case_when(Subnat_count>=25 ~ Subnat_Freq,
                                                                 Subnat_count<25 ~ Nat_Freq ))
  
  method.sector$Source <- paste(year , survey, sep=" ")
  method.sector$Population <- population
  method.sector$Country.Method <- paste(method.sector$region , ":", method.sector$modern_method_source, sep="")
  
  method.sector <- method.sector %>% select(region, Source, Population, modern_method_source, Country.Method, sector, Frequency, Subnat_count)
  method.sector <- method.sector %>% spread(sector, Frequency) 
  
  method.sector <- method.sector %>% mutate(Notes= case_when(Subnat_count>=25 ~ "",
                                                             Subnat_count<25 ~ "Used National Average because samples too small @ regional level" ))
  
  # Add in NGO column if not already present
  cols <- c(NGO= NA_real_)
  method.sector <- add_column(method.sector, !!!cols[setdiff(names(cols), names(method.sector))])
  
  mysectorlist <- c("Public", "NGO", "Private Clinic", "Pharmacy"  ,  "Shop Church Friend",  "Other" )
  method.sector <- method.sector %>% select(region, Source, Population, modern_method_source, Country.Method, colnames(method.sector)[which(colnames(method.sector) %in% mysectorlist)], Notes, Subnat_count) %>% rename(N=Subnat_count)
  
  ###############################################################
  ###############################################################
  ###############################################################
  write.xlsx(as.data.frame(method.sector), paste(folder, "/", region_name,"_", year, "_SubnationalEMUPrep.xlsx" , sep=""), sheetName="Source Data", col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)
  #write.xlsx(as.data.frame(regional.distribution), paste(folder, "/", region_name,"_", year, "_SubnationalEMUPrep.xlsx" , sep=""), sheetName="Regional Distribution", col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
  #write.xlsx(as.data.frame(method.prev), paste(folder, "/", region_name,"_", year, "_SubnationalEMUPrep.xlsx" , sep=""), sheetName="Method Mix", col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
  
  return(head(method.sector))
}

# Collapse sectors into categories
collapse_sectors_ss2emu <- function(my_data) {
  my_data <- my_data %>% mutate(sector_category = ifelse(Source %in% c("Private.Clinic","Pharmacy", "NGO"),"commercial_medical",
                                                         ifelse(Source %in% c("Shop.Church.Friend", "Other"), "other","public")))
  return(my_data)
}

# Calculate category proportion
calculate_category_prop <- function(my_data, gps_cols = c("centroid_long","centroid_lat")) {
  
  tmp <- my_data %>%
    distinct(country_name, region, Year, Source, modern_method_source, proportion) 
  
 tmp <- collapse_sectors_ss2emu(tmp)
  
  tmp <- tmp %>%
    group_by(region, Year, modern_method_source, country_name, sector_category) %>%
    mutate(category_proportion = sum(proportion, na.rm=TRUE))
  
  updated_data <- merge(tmp, my_data)
  
  updated_data <- updated_data %>%
    dplyr::select(country_name, region, Year, sector_category,  modern_method_source, category_proportion, Population, Notes, N, gps_cols) %>%
    distinct()
  
  return(updated_data)
}

# Find centroids of GPS clusters

gps_to_cartesian <- function(my_data) {
  r = 6371 # earth radius
  
  updated_coords_data <- my_data %>%
    mutate(long_rad = long * pi/180) %>%
    mutate(lat_rad = lat * pi/180) %>%
    group_by(region, Country) %>%
    mutate(x_coord = r * cos(lat_rad) * cos(long_rad)) %>% # cartesian
    mutate(y_coord = r * cos(lat_rad) * sin(long_rad)) %>%
    mutate(z_coord = r * sin(lat_rad)) %>%
    mutate(x_av = mean(x_coord)) %>% # compute centoid
    mutate(y_av = mean(y_coord)) %>%
    mutate(z_av = mean(z_coord)) %>%
    mutate(centroid_lat = asin(z_av/r)* 180/pi) %>% # convert average x,y,z back to long/lat
    mutate(centroid_long = atan2(y_av, x_av)* 180/pi) %>% 
    distinct(region, Country, Year, centroid_long, centroid_lat)
  
  return(updated_coords_data)
}

# functions for creating basis functions-----------------------------------------------------
tpower <- function(x, t, p) {
  # Truncated p-th power function
  (x - t)^p * (x > t)
}

bbase <- function(x, xdat = NULL, xl = min(x), xr = max(x), deg = 3, dx = 0.1) {
  # Construct B-spline basis
  if (is.null(xdat)) 
    xdat = x
  
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  return(list(B.ik = B, knots.k = knots))
}

# Make data; function to create universal df of estimates -------------------
create_subnational_df <- function(myfile_names, mydf.list, mycountry_codes, mycountry_codes_DHS){
  updated.df <- tibble()
  for(i in 1:length(myfile_names)) {
    print(i)
    countrycode <- names(mydf.list[i])
    countrycode <- substr(countrycode, 1,2)
    mydata <- tibble(mydf.list[[i]])
    mydata <- mydata %>% mutate(country_code = mycountry_codes[i], 
                                year = as.numeric(year),
                                n = as.numeric(n))
    if(countrycode=="IA") {
      country_name <- "India" } else {
      if(countrycode=="PG") {
        country_name <- "Papua New Guinea" } else {
        country_name <- mycountry_codes_DHS %>% filter(Code==mycountry_codes[i]) %>% dplyr::select(`Country Name`) }
      }
    country_name <- as.vector(unlist(country_name))
    mydata <- mydata %>% mutate(Country = country_name)
    updated.df <- bind_rows(mydata, updated.df)
  }
  return(updated.df)
}

# Running jags model functions --------------------------------------------------------------
# Examine trends in lagged values
rates_of_change_lagged <- function(my_data) {
  new_data <- my_data %>%
    arrange(average_year) %>%
    group_by(Method,Country, region) %>%
    mutate(lag.logitcm_tp_ratio = lag(logitcm_tp_ratio, n = 1, default = NA)) %>% # lagged column for subtraction
    mutate(diff.logitcm_tp_ratio = logitcm_tp_ratio - lag.logitcm_tp_ratio) %>% # changes between years
    mutate(rate.logitcm_tp_ratio = diff.logitcm_tp_ratio/(average_year - lag(average_year))) # rates of change
  return(new_data)
}

# Examine trends in lagged values
prop_pub_time_point_diff <- function(my_data) {
  new_data <- my_data %>%
    group_by(Method, Country, region) %>%
    arrange(average_year) %>%
    mutate(lag.logitpublic = lag(logitpublic, n = 1, default = NA)) %>% # lagged column for subtraction
    mutate(diff.logitpublic = logitpublic - lag.logitpublic) %>% # changes between years
    mutate(rate.logitpublic = diff.logitpublic/(average_year - lag(average_year))) # rates of change
  return(new_data)
}

flat_cor_mat <- function(cor_r){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 3 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_r <- cor_r %>% distinct(cor, .keep_all = TRUE)
  cor_r$cor <- round(cor_r$cor,1)
  cor_r <- cor_r %>% filter(cor!=1)
  return(cor_r)
}

#### Collapsing Private sector breakdown into categories
collapse_sectors <- function(my_data) {
  my_data <- my_data %>% mutate(sector_category = ifelse(Sector %in% c("private hospital/clinic","pharmacy","private doctor"),"commerical_medical", 
                                                         ifelse(Sector %in% c("shop", "friend/relative", "non-medical sources", " non-medical sources"), "commercial_non-medical",
                                                                ifelse(Sector=="family planning clinic", "NGO",
                                                                       ifelse(Sector=="church", "NPO",
                                                                              ifelse(Sector %in% c("other private medical source", "other unspecified source"), "other","public"))))))
  my_data <- my_data %>%
    filter(is.na(Proportion)==FALSE) %>%
    group_by(Country, average_year, sector_category, Method) %>%
    mutate(total_sector_prop = sum(Proportion))
  
  return(my_data)
}

# Standard naming
standard_method_names <- function(my_data) {
  #levels(my_data$Method) <- c(levels(my_data$Method), "Condom", "Female Sterilization", "Male Sterilization", "Other Modern Methods", "OC Pills")
  my_data <- my_data %>%
    mutate(Method = replace(Method, Method == "Condom (m+f)", "Condom")) %>%
    mutate(Method = replace(Method, Method == "male condom", "Condom")) %>%
    mutate(Method = replace(Method, Method == "condom", "Condom")) %>%
    #mutate(Method = replace(Method, Method == "female condom", "Condom")) %>%
    mutate(Method = replace(Method, Method == "Sterilization (female)", "Female Sterilization")) %>%
    mutate(Method = replace(Method, Method == "female sterilization", "Female Sterilization")) %>%
    mutate(Method = replace(Method, Method == "Female sterilization", "Female Sterilization")) %>%
    mutate(Method = replace(Method, Method == "Sterilization (male)", "Male Sterilization")) %>%
    mutate(Method = replace(Method, Method == "Male sterilization", "Male Sterilization")) %>%
    mutate(Method = replace(Method, Method == "male sterilization", "Male Sterilization")) %>%
    mutate(Method = replace(Method, Method == "Pill", "OC Pills")) %>%
    mutate(Method = replace(Method, Method == "pill", "OC Pills")) %>%
    mutate(Method = replace(Method, Method == "injections", "Injectables")) %>%
    mutate(Method = replace(Method, Method == "implants", "Implants")) %>%
    mutate(Method = replace(Method, Method %in% c("Other", "Other Modern Methods", "LAM", "diaphragm", "female condom", "foam or jelly", "standard days method", "diaphragm, foam or jelly", "lactational amenorrhea", "emergency contraception", "other modern methods"), "Other Modern Methods")) #%>%
  #as.data.frame()
  
  return(my_data)
}

### Calculating super_region_index 
superregion_index_fun <- function(my_data, n_region) {
  my_data$index_superregion <- NA
  for (i in 1:length(n_region)) {
    for (j in 1:nrow(my_data)) {
      region_name <- n_region[i]
      if(my_data$Super_region[j]==region_name) {
        my_data$index_superregion[j] <- i
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}

### Calculating subnat_index 
subnat_index_fun <- function(my_data, my_subnat, my_country) {
  for (i in 1:length(my_subnat)) {
    for (j in 1:nrow(my_data)) {
      subnat_name <- my_subnat[i]
      country_name <- my_country[i]
      if(my_data$Country[j]==country_name & my_data$Region[j]==subnat_name) {
        my_data$index_subnat[j] <- i
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}

### Calculating country_index 
country_index_fun <- function(my_data, my_country) {
  for (i in 1:length(my_country)) {
    for (j in 1:nrow(my_data)) {
      country_name <- my_country[i]
      if(my_data$Country[j]==country_name) {
        my_data$index_country[j] <- i
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}

### Calculating continent_index 
continent_index_fun <- function(my_data, n_con) {
  my_data$index_continent <- NA
  for (i in 1:length(n_con)) {
    for (j in 1:nrow(my_data)) {
      con_name <- n_con[i]
      if(my_data$continent[j]==con_name) {
        my_data$index_continent[j] <- i
      } else {
        next
      }
    }
  }
  return(my_data)
}

### Calculating method_index 
method_index_fun <- function(my_data, my_methods) {
  my_data$index_method <- rep(NA, nrow(my_data))
  for (i in 1:length(my_methods)) {
    for (j in 1:nrow(my_data)) {
      method_name <- my_methods[i]
      if(my_data$Method[j]==method_name) {
        my_data$index_method[j] <- i
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}

### Calculating sector_index 
sector_index_fun <- function(my_data, my_sectors) {
  for (i in 1:length(my_sectors)) {
    for (j in 1:nrow(my_data)) {
      sector_name <- my_sectors[i]
      if(my_data$sector_category[j]==sector_name) {
        my_data$index_sector[j] <- i
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}

### Reversing method_index
method_collapse_name_fun <- function(my_data) {
  for (i in 1:length(n_method)) {
    for (j in 1:nrow(my_data)) {
      method_name <- n_method[i]
      if(my_data$index_method_collapse[j]==i) {
        my_data$Method_collapse[j] <- method_name
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}

### Reversing country_index
country_name_fun <- function(my_data,my_country) {
  for (i in 1:length(my_country)) {
    for (j in 1:nrow(my_data)) {
      country_name <- my_country[i]
      if(my_data$index_country[j]==i) {
        my_data$Country[j] <- as.character(country_name)
      } 
      else {
        next
      }
    }
  }
  return(my_data)
}


#### Collapsing Some Methods into 'Others' category
collapse_methods_fun <- function(my_data) {
  my_data$Method_collapse <- as.vector(sapply(my_data[,'Method'], FUN = function(x)
    ifelse(x=="Implants","Implants", 
           ifelse(x=="Injectables","Injectables",
                  ifelse(x=="IUD","IUD",
                         ifelse(x %in% c("Female Sterilization","Female sterilization"),"Female sterilization",
                                ifelse(x %in% c("Male Sterilization","Male sterilization"),"Male sterilization",
                                       ifelse(x %in% c("OC Pills", "Pill"),"Pill",
                                              ifelse(x=="Condom","Condom", "Others"
                                              )))))))) )
  
  return(my_data)
}

#### Get splines 

GetSplines <- function( # Get B-splines
  x.i, ##<< Vector of x-values (without NAs) for which splines need to be calculated (determines the number of rows of B.ik)
  x0 = NULL, ##<< x-value which determines knot placement. By default, knot is placed half-interval before last observation
  I , ##<< Interval length between two knots during observation period
  degree = 3 # currently tested only with degree 3
) {
  if (is.null(x0)) {
    x0 <- max(x.i)-0.5*I 
  } 
  # get knots, given that one knot needs to be in year0
  knots <- seq(x0-I, x0+I, I) 
  while (min(x.i) < knots[1]) knots <- c(seq(knots[1]-I, knots[1]-I,I), knots)
  while (max(x.i) > knots[length(knots)]) knots <- c(knots, seq(knots[length(knots)]+I, 
                                                                knots[length(knots)]+I, I)) 
  Btemp.ik <- bs(x.i, knots = knots[-c(1, length(knots))],  degree = degree,
                 Boundary.knots = knots[c(1, length(knots))])
  indicesofcolswithoutzeroes <- which(apply(Btemp.ik, 2, sum) > 0) # only remove columns with zeroes at start and end
  startnonzerocol <- indicesofcolswithoutzeroes[1]
  endnonzerocol <- indicesofcolswithoutzeroes[length(indicesofcolswithoutzeroes)]
  B.ik <- Btemp.ik[,startnonzerocol:endnonzerocol]
  colnames(B.ik) <- paste0("spline", seq(1, dim(B.ik)[2]))
  knots.k <- knots[startnonzerocol:endnonzerocol]
  names(knots.k) <- paste0("spline", seq(1, dim(B.ik)[2]))
  ##value<< List of B-splines containing:
  return(list(B.ik = B.ik, ##<< Matrix, each row is one observation, each column is one B-spline.
              knots.k = knots.k ##<< Vector of knots.
  ))
}

# LOO median and CI p calc for all sectors
loo_p_calc_sectors <- function(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_loo_model) { 
  p_med_loo <- P_median_calc_sectors(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_loo_model) # Median estimate for proportion of modern contraceptives coming from all sector
  P_80quantile <- iq80_P_calc_sectors(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_loo_model) # Mean estimate for proportion of modern contraceptives coming from public sector
  P_95quantile <- iq95_P_calc_sectors(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_loo_model) # Mean estimate for proportion of modern contraceptives coming from public sector
  loo_all_p <- merge(P_80quantile,p_med_loo)
  loo_all_p <- merge(P_95quantile,loo_all_p)
  return(loo_all_p)
}

# Calculation for sector data
P_median_calc_sectors <- function(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_model) { # Median alpha values
  
  years <- unique(year_index_table$floored_year)
  n_years <- years %>% length() # important when using 6-monthly description
  
  #### P median
  
  P_samp <- my_model$BUGSoutput$sims.list$P
  P_dims <- dim(P_samp)
  print(P_dims)

  P_s_med <- array(dim=c(P_dims[4],n_years,P_dims[3],P_dims[2])) # 30x5 matrix for CountryxMethod then 5 into an array for each sector

  # Create a table for storing individual true country public data
  for(k in 1:n_years) { # time loop
    year <- years[k]
    time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
      dplyr::filter(floored_year == year) %>%
      dplyr::select(index_year) %>%
      unlist() %>%
      as.vector()
    for(j in 1:P_dims[4]) { # subnat
      for (r in 1:P_dims[3]) { # method
        for (i in 1:P_dims[2]) { # sector
          P_s_med[j,k,r,i] <- stats::median(P_samp[,i,r,j,time_index])
        }
      }
    }
  }
  
  P_s_med <- plyr::adply(P_s_med, c(1,2,3,4)) 
  colnames(P_s_med) <- c("index_subnat", "index_year", "index_method", "index_sector", "median_p")
  
  P_s_med <- P_s_med %>%
    mutate(index_method = as.numeric(index_method)) %>%
    mutate(index_sector = as.numeric(index_sector)) %>%
    mutate(index_subnat = as.numeric(index_subnat)) %>%
    mutate(index_year = as.numeric(index_year)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(sector_index_table)
  
  return(P_s_med)
}

iq95_P_calc_sectors <- function(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_model) { # Median alpha values
  
  years <- unique(year_index_table$floored_year)
  n_years <- years %>% length() # important when using 6-monthly description
  

  P_samp <- my_model$BUGSoutput$sims.list$P
  P_dims <- dim(P_samp)
  print(P_dims)
  
  upper_P_s_med <- lower_P_s_med <- array(dim=c(P_dims[4],n_years,P_dims[3],P_dims[2])) # 30x5 matrix for CountryxMethod then 5 into an array for each sector
  
  for(t in 1:n_years) { # time loop
    year <- years[t]
    time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
      dplyr::filter(floored_year == year) %>% 
      dplyr::select(index_year) %>% 
      unlist() %>% 
      as.vector()
    for(j in 1:P_dims[4]) {
      for (r in 1:P_dims[3]) { # Create a table for storing individual true country public data
        for (i in 1:P_dims[2]) {
          upper_P_s_med[j,t,r,i] <- stats::quantile(P_samp[,i,r,j,time_index], probs = 0.975, na.rm=TRUE)
        }
      } 
    }
  }
  
  upper_P_s_med <- plyr::adply(upper_P_s_med, c(1,2,3,4)) 
  colnames(upper_P_s_med) <-  c("index_subnat", "index_year", "index_method", "index_sector", "upper_q95_P")
  
  for(t in 1:n_years) { # time loop
    year <- years[t]
    time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
      dplyr::filter(floored_year == year) %>% 
      dplyr::select(index_year) %>% 
      unlist() %>% 
      as.vector()
    for(j in 1:P_dims[4]) {
      for (r in 1:P_dims[3]) { # Create a table for storing individual true country public data
        for (i in 1:P_dims[2]) {
          lower_P_s_med[j,t,r,i] <- stats::quantile(P_samp[,i,r,j,time_index], probs = 0.025, na.rm=TRUE)
        }
      } 
    }
  }
  
  lower_P_s_med <- plyr::adply(lower_P_s_med, c(1,2,3,4)) 
  colnames(lower_P_s_med) <-  c("index_subnat", "index_year", "index_method", "index_sector", "lower_q95_P")
  
  P_s_med <- merge(lower_P_s_med, upper_P_s_med)
  
  P_s_med <- P_s_med %>%
    mutate(index_method = as.numeric(index_method)) %>%
    mutate(index_sector = as.numeric(index_sector)) %>%
    mutate(index_subnat = as.numeric(index_subnat)) %>%
    mutate(index_year = as.numeric(index_year)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(sector_index_table)
  
  return(P_s_med)
}

iq80_P_calc_sectors <- function(subnat_index_table, method_index_table, sector_index_table, year_index_table, my_model) { # Median alpha values
  
  years <- unique(year_index_table$floored_year)
  n_years <- years %>% length() # important when using 6-monthly description
  
  
  P_samp <- my_model$BUGSoutput$sims.list$P
  P_dims <- dim(P_samp)
  print(P_dims)
  
  upper_P_s_med <- lower_P_s_med <- array(dim=c(P_dims[4],n_years,P_dims[3],P_dims[2])) # 30x5 matrix for CountryxMethod then 5 into an array for each sector
  
  for(t in 1:n_years) { # time loop
    year <- years[t]
    time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
      dplyr::filter(floored_year == year) %>% 
      dplyr::select(index_year) %>% 
      unlist() %>% 
      as.vector()
    for(j in 1:P_dims[4]) {
      for (r in 1:P_dims[3]) { # Create a table for storing individual true country public data
        for (i in 1:P_dims[2]) {
          upper_P_s_med[j,t,r,i] <- stats::quantile(P_samp[,i,r,j,time_index], probs = 0.9, na.rm=TRUE)
        }
      } 
    }
  }
  
  upper_P_s_med <- plyr::adply(upper_P_s_med, c(1,2,3,4)) 
  colnames(upper_P_s_med) <-  c("index_subnat", "index_year", "index_method", "index_sector", "upper_q80_P")
  
  for(t in 1:n_years) { # time loop
    year <- years[t]
    time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
      dplyr::filter(floored_year == year) %>% 
      dplyr::select(index_year) %>% 
      unlist() %>% 
      as.vector()
    for(j in 1:P_dims[4]) {
      for (r in 1:P_dims[3]) { # Create a table for storing individual true country public data
        for (i in 1:P_dims[2]) {
          lower_P_s_med[j,t,r,i] <- stats::quantile(P_samp[,i,r,j,time_index], probs = 0.1, na.rm=TRUE)
        }
      } 
    }
  }
  
  lower_P_s_med <- plyr::adply(lower_P_s_med, c(1,2,3,4)) 
  colnames(lower_P_s_med) <-  c("index_subnat", "index_year", "index_method", "index_sector", "lower_q80_P")
  
  P_s_med <- merge(lower_P_s_med, upper_P_s_med)
  
  P_s_med <- P_s_med %>%
    mutate(index_method = as.numeric(index_method)) %>%
    mutate(index_sector = as.numeric(index_sector)) %>%
    mutate(index_subnat = as.numeric(index_subnat)) %>%
    mutate(index_year = as.numeric(index_year)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(sector_index_table)
  
  return(P_s_med)
}

# Calculation for sector data
P_predictive_samp_sectors <- function(match_subnat_test, subnat_index_table,  match_method_test, method_index_table, match_years_test, pred_samp) { # Median alpha values
  
  # #### P median
  pred_samp <- apply(pred_samp, c(2,3,4,5), median, na.rm=TRUE) 
  pred_samp <- plyr::adply(pred_samp, c(2,3,4))
  pred_samp <- as_tibble(pred_samp)
  colnames(pred_samp) <- c("index_method","index_subnat", "index_year", "Public", "Commercial_medical", "Other")
  pred_samp <- pred_samp %>%
    mutate(index_method = as.numeric(index_method),
           index_subnat = as.numeric(index_subnat),
           index_year = as.numeric(index_year)) %>%
    pivot_longer(cols = c(Public, Commercial_medical, Other), names_to = "Sector", values_to = "median_p")
  test_obs <- dplyr::tibble(index_subnat = as.numeric(match_subnat_test), # get all test obs
                  index_year = as.numeric(match_years_test),
                  index_method = as.numeric(match_method_test))

  pred_samp_test <- test_obs %>% left_join(pred_samp) # keep test obs only

  pred_samp_test <- pred_samp_test %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table)
  
  return(pred_samp_test)
}

iq_P_predictive_samp_sectors <- function(match_subnat_test, subnat_index_table,  match_method_test, method_index_table, match_years_test, sector_index_table, pred_samp) { # Median alpha values

  pred_samp_95 <- apply(pred_samp, c(2,3,4,5),  quantile, probs = c(0.025,0.975),  na.rm = TRUE)
  pred_samp_95 <- plyr::adply(pred_samp_95, c(2,3,4,5))
  colnames(pred_samp_95) <- c("index_sector", "index_method", "index_subnat", "index_year", "lower_q95_P", "upper_q95_P") 
  pred_samp_95 <- pred_samp_95 %>%
    mutate(index_method = as.numeric(index_method),
           index_sector = as.numeric(index_sector),
           index_subnat = as.numeric(index_subnat),
           index_year = as.numeric(index_year))
  
  test_obs <- dplyr::tibble(index_sector = rep(1:3, each=length(match_subnat_test)),
                            index_method = rep(match_method_test, 3),
                            index_subnat = rep(match_subnat_test, 3),
                            index_year = rep(match_years_test,3))
  
  pred_samp_95 <- test_obs %>% left_join(pred_samp_95) # keep test obs only
  
  pred_samp_80 <- apply(pred_samp, c(2,3,4,5),  quantile, probs = c(0.1,0.9),  na.rm = TRUE)
  pred_samp_80 <- plyr::adply(pred_samp_80, c(2,3,4,5))
  colnames(pred_samp_80) <- c("index_sector", "index_method", "index_subnat", "index_year", "lower_q80_P", "upper_q80_P") 
  pred_samp_80 <- pred_samp_80 %>%
    mutate(index_method = as.numeric(index_method),
           index_sector = as.numeric(index_sector),
           index_subnat = as.numeric(index_subnat),
           index_year = as.numeric(index_year))
  
  test_obs <- dplyr::tibble(index_sector = rep(1:3, each=length(match_subnat_test)),
                            index_method = rep(match_method_test, 3),
                            index_subnat = rep(match_subnat_test, 3),
                            index_year = rep(match_years_test,3))
  
  pred_samp_80 <- test_obs %>% left_join(pred_samp_80) # keep test obs only
  
  pred_samp_CI <- left_join(pred_samp_95, pred_samp_80)
  
  pred_samp_CI <- pred_samp_CI %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(sector_index_table)
  
  return(pred_samp_CI)
}

iq80_P_predictive_samp_sectors <- function(country_index_table, method_index_table, sector_index_table, my_pred_samps) { # Median alpha values
  
  #### P CI
  
  P_samp <- my_pred_samps
  P_dims <- dim(P_samp)
  print(P_dims)
  
  upper_P_s_med <- lower_P_s_med <- array(dim=c(P_dims[4],P_dims[5],P_dims[3],P_dims[2])) # 30x5 matrix for CountryxMethod then 5 into an array for each sector
  
  for(k in 1:P_dims[5]) {
    for(j in 1:P_dims[4]) {
      for (r in 1:P_dims[3]) { # Create a table for storing individual true country public data
        for (i in 1:P_dims[2]) {
          upper_P_s_med[j,k,r,i] <- stats::quantile(P_samp[,i,r,j,k], probs = 0.9, na.rm = TRUE)
        }
      } 
    }
  }
  
  upper_P_s_med <- plyr::adply(upper_P_s_med, c(1,2,3,4)) 
  colnames(upper_P_s_med) <- c("index_subnat", "index_year", "index_method", "index_sector", "upper_q80_P")
  
  for(k in 1:P_dims[5]) {
    for(j in 1:P_dims[4]) {
      for (r in 1:P_dims[3]) { # Create a table for storing individual true country public data
        for (i in 1:P_dims[2]) {
          lower_P_s_med[j,k,r,i] <- stats::quantile(P_samp[,i,r,j,k], probs = 0.1, na.rm = TRUE)
        }
      } 
    }
  }
  
  lower_P_s_med <- plyr::adply(lower_P_s_med, c(1,2,3,4)) 
  colnames(lower_P_s_med) <- c("index_subnat", "index_year", "index_method", "index_sector", "lower_q80_P")
  
  P_s_med <- merge(lower_P_s_med, upper_P_s_med)
  
  P_s_med <- P_s_med %>%
    mutate(index_method = as.numeric(index_method)) %>%
    mutate(index_sector = as.numeric(index_sector)) %>%
    mutate(index_subnat = as.numeric(index_subnat)) %>%
    mutate(index_year = as.numeric(index_year)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(sector_index_table)
  
  return(P_s_med)
}

# LOO median and CI p calc for all sectors
loo_P_predictive_samp_sectors <- function(match_subnat_test, subnat_index_table,  match_method_test, method_index_table, match_years_test, sector_index_table, my_pred_samps) { 
  
  #### Median P calc
  p_med_loo <- P_predictive_samp_sectors(match_subnat_test, subnat_index_table,  match_method_test, method_index_table, match_years_test, my_pred_samps) # Median estimate for proportion of modern contraceptives coming from all sectors

  P_quantile <- iq_P_predictive_samp_sectors(match_subnat_test, subnat_index_table,  match_method_test, method_index_table, match_years_test, sector_index_table, my_pred_samps) # Mean estimate for proportion of modern contracetievs coming from public sector

  loo_all_p <- merge(P_quantile,p_med_loo)
  loo_all_p <-standard_method_names(loo_all_p)
  
  return(loo_all_p)
  
}

loo_CI_check_sectors <- function(my_CI_calc,my_test_data) { 
  my_test_data <-  loo_test_data %>% 
    dplyr::select(Country, Region, Method, average_year, Commercial_medical, Other, Public) 
  my_test_data_long <- gather(my_test_data, sector_category, observed_proportion, c(Public, Commercial_medical, Other), factor_key=TRUE) %>%
    rename(Sector = sector_category)
  test <- dplyr::left_join(pred_samples, my_test_data_long)  
  test <- test %>% 
    mutate(result95 = ifelse(round(lower_q95_P, 2) <= round(observed_proportion,2) & round(observed_proportion,2)  <= round(upper_q95_P, 2),1,0)) %>%
    mutate(result80 = ifelse(round(lower_q80_P, 2) <= round(observed_proportion,2) & round(observed_proportion,2)  <= round(upper_q80_P, 2),1,0))
    
  print(test)
  
  # Addressing issues with practically 0 estimates and observations
  for(i in 1:nrow(test)) {
    if(test$observed_proportion[i]<=0.005) {
      if(test$lower_q95_P[i] < 0.005) {
        test$result95[i] <- 1
      }
    } else {
      next
    }
    if(test$observed_proportion[i]<=0.005) {
      if(test$lower_q80_P[i] < 0.005) {
        test$result80[i] <- 1
      }
    } else {
      next
    }
  }
  
  print("The average 95% performance of your model is:")
  print(sum(test$result95)/length(test$result95))
  
  print("The average 80% performance of your model is:")
  print(sum(test$result80)/length(test$result80))
  
  return(test)
}

#### Leave one out validation for all sectors
get_leave_one_out_data <- function(my_data) {
  
  suitable_data <-  mytestdata <- NULL
  
  country_subnat_tbl <- my_data %>% ungroup() %>% dplyr::select(Country, Region) %>% distinct()
  my_methods <- as.vector(unique(my_data$Method))
  
  for (i in 1:nrow(country_subnat_tbl)) {
    c <- country_subnat_tbl$Country[i]
    s <- country_subnat_tbl$Region[i]
    for (m in my_methods) {
      sub_data <- my_data %>% dplyr::filter(Country==c & Region==s & Method==m)
      if(nrow(sub_data)>0) {
        sector_props <- sub_data[,c("Public","Commercial_medical","Other")] # Taking proportion data to examine for NAs etc.
        if(nrow(sector_props)==2 & sum(is.na(sector_props))==3) { # skip if there are 2 entries and half missing data
          next
        } else { if(nrow(sector_props)==1 | sum(is.na(sector_props))==ncol(sector_props)) { # Remove combinations of just NA or one liners
          next
        } else {
          suitable_data <- rbind(suitable_data,sub_data)
        }
        }
      } else { next }
    }
  }
  
  new_country_subnat_tbl <- suitable_data %>% ungroup() %>% dplyr::select(Country, Region, Method) %>% distinct()

  for (i in 1:nrow(new_country_subnat_tbl)) {
    c <- new_country_subnat_tbl$Country[i]
    s <- new_country_subnat_tbl$Region[i]
    m <- new_country_subnat_tbl$Method[i]
    my_new_data <- tibble(suitable_data) %>% filter(Country==c & Region==s & Method==m)
    if(nrow(my_new_data)>0) {
      max_year <- max(my_new_data$average_year, na.rm=TRUE)
      max_data <- tibble(my_new_data) %>% dplyr::filter(average_year==max_year)
      mytestdata <- rbind(mytestdata,max_data)
    } else {
      next
    }
  }
  
  print("done test data")
  my_train_data <- data.frame(suitable_data) %>% dplyr::anti_join(data.frame(mytestdata)) # removing test data from training
  print("done train data")
  loo_data_list <- list(train_data = my_train_data, test_data = mytestdata)
  return(loo_data_list)
}

calculate_subnatSE_data <- function(filepath, myresultsfolder, strata_null = FALSE) {
  # Subnational Inputs for SStoEMU
  # Kristin Bietsch, PhD
  # Avenir Health
  # 08/08/19
  
  ##############################################
  # Set Results Location
  ########################################################
  
  folder <- myresultsfolder
  
  ########################################################
  # Analysis of most recent survey 
  ########################################################
  country_codes_DHS <- readxl::read_excel("data/country_codes_DHS.xlsx") # country codes and names
  
  data <- read_dta(file=filepath)
  year <-  min(data$v007)
  survey <- "DHS"
  population <- "AW"
  country_code <- unique(data$v000)
  year
  val_labels(data$v024)
  #print(data$sdevreg)
  
  data$region <- as_factor(data$v024, levels = "labels")
  
  # Percentage currently using a modern method
  data <- data %>% mutate(mcpr= case_when(v313==3 ~ 1, TRUE ~ 0))
  
  # Percentage married
  data <- data %>% mutate(married= case_when(v502==1 ~ 1, TRUE ~ 0))
  
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

calculate_subnatIndiaSE_data <- function(filepath, myresultsfolder, strata_null = FALSE) {
  # Subnational Inputs for SStoEMU
  # Kristin Bietsch, PhD
  # Avenir Health
  # 08/08/19
  
  ##############################################
  # Set Results Location
  ########################################################
  
  folder <- myresultsfolder
  
  ########################################################
  # Analysis of most recent survey 
  ########################################################
  country_codes_DHS <- readxl::read_excel("data/country_codes_DHS.xlsx") # country codes and names
  
  data <- read_dta(file=filepath)
  year <-  min(data$v007)
  survey <- "DHS"
  population <- "AW"
  country_code <- unique(data$v000)
  year
  val_labels(data$v024)
  
  if(country_code=="IA3") {
    state_code <- substr(filepath, 16,17) # region code from file name
    state_name <- as.vector(unlist(country_codes_DHS %>% filter(`India States`==state_code) %>% select(`State Name`)))
  }
  
  data$region <- ifelse(country_code=="IA3",state_name,as_factor(data$v024, levels = "labels"))
  
  # Percentage currently using a modern method
  data <- data %>% mutate(mcpr= case_when(v313==3 ~ 1, TRUE ~ 0))
  
  # Percentage married
  data <- data %>% mutate(married= case_when(v502==1 ~ 1, TRUE ~ 0))
  
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
  
  write.xlsx(as.data.frame(my_SEdf), paste(folder, "/", country_code,"_", year, "_", state_name,"_SEdf.xlsx" , sep=""))
  
  return(head(my_SEdf))
}

bs_bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(splines::bs(x, knots = knots, degree = deg, Boundary.knots = c(knots[1], knots[length(knots)])), nrow = length(x))
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  return(list(B.ik = bs_matrix, ##<< Matrix, each row is one observation, each column is one B-spline.
              knots.k = knots[-c(1,2,length(knots),(length(knots)-1))] ##<< Vector of knots.
  ))
}

# Flatten correlation matrix to vector 
flat_cor_mat <- function(cor_r){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 3 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_r <- cor_r %>% distinct(cor, .keep_all = TRUE)
  cor_r$cor <- round(cor_r$cor,1)
  cor_r <- cor_r %>% filter(cor!=1)
  return(cor_r)
}

bs_bbase_precise <- function(x = x,lastobs = max(x), xl = min(x), xr = max(x), nseg = 10, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Compute position of knot before last observation
  dk <- lastobs
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Find index of closest knot to dk
  dk_index <- which.min(abs(knots-dk))
  # Find transformation to knot placement so that dk is a knot 
  ktrans <- (dk-knots)[dk_index]
  # Add transformation to knots
  knotsnew <- knots + ktrans
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(splines::bs(x, knots = knotsnew, degree = deg, Boundary.knots = c(knotsnew[1], knotsnew[length(knotsnew)])), nrow = length(x))
  
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  used_knots <- knotsnew[-c(1,2,length(knotsnew),(length(knotsnew)-1))]
  Kstar <- which(used_knots==dk)
  
  return(list(B.ik = bs_matrix, ##<< Matrix, each row is one observation, each column is one B-spline.
              knots.k = used_knots, ##<< Vector of transformed knots.
              Kstar = Kstar # Knot point of last observation
  ))
}


# # Calculation for sector data
P_median_global_2files <- function(subnat_index_table,
                            method_index_table,
                            sector_type,
                            year_index_table,
                            file_chain1,
                            file_chain2) {

  mod_P1 <- readRDS(file_chain1)
  mod_P2 <- readRDS(file_chain2)

  P_dims <-  dim(mod_P1)

  time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
    dplyr::filter(average_year>floored_year) %>%
    dplyr::select(index_year) %>%
    unlist() %>%
    as.vector()

  averageyear_index_table <- year_index_table %>%
    dplyr::filter(average_year>floored_year) %>%
    select(average_year, index_year) %>%
    dplyr::mutate(index_year = 1:n())

  P_s_med <- array(dim=c(P_dims[3], length(time_index), P_dims[2])) # subnat, time, method

  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_s_med[j,k,r] <- stats::median(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]))
      }
    }
  }

  # P_dims <- dim(Psamps)
  #
  # #### P median
  # P_s_med <- array(dim=c(length(time_index),P_dims[2],P_dims[3])) # method, year, subnat
  #
  # # Create a table for storing individual true country public data
  # for(k in 1:length(time_index)) { # time loop
  #   for(s in 1:P_dims[3]) { # subnat
  #     for (m in 1:P_dims[2]) { # method
  #       P_s_med[k,m,s] <- stats::median(Psamps[,m,s,time_index[k]])
  #     }
  #   }
  # }

  P_s_med <- plyr::adply(P_s_med, c(1,2,3))
  colnames(P_s_med) <- c("index_subnat", "index_year", "index_method",  "median_p")

  P_s_med <- P_s_med %>%
    group_by(index_method, index_subnat, index_year) %>%
   # summarise(median_p = stats::median(estimate)) %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(averageyear_index_table) %>%
    dplyr::mutate(Sector = sector_type)

  P_Q <- array(dim=c(P_dims[3], length(time_index), P_dims[2], 4)) # subnat, time, method, quantile(95, 80)

  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_Q[j,k,r, 1:2] <- as.vector(unlist(stats::quantile(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]), prob=c(0.025, 0.975))))
        P_Q[j,k,r, 3:4] <- as.vector(unlist(stats::quantile(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]), prob=c(0.1, 0.9))))
      }
    }
  }

  P_Q <- plyr::adply(P_Q, c(1,2,3))
  colnames(P_Q) <- c("index_subnat", "index_year", "index_method",  "lower_95", "uppper_95", "lower_80", "upper_80")

  P_s_med <- P_Q %>%
    group_by(index_method, index_subnat, index_year) %>%
    # summarise(median_p = stats::median(estimate)) %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    left_join(P_s_med)

  return(P_s_med)
}

# # Calculation for sector data
P_sd_global_2files <- function(subnat_index_table,
                                   method_index_table,
                                   sector_type,
                                   year_index_table,
                                   file_chain1,
                                   file_chain2) {
  
  mod_P1 <- readRDS(file_chain1)
  mod_P2 <- readRDS(file_chain2)
  
  P_dims <-  dim(mod_P1)
  
  time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
    dplyr::filter(average_year>floored_year) %>%
    dplyr::select(index_year) %>%
    unlist() %>%
    as.vector()
  
  averageyear_index_table <- year_index_table %>%
    dplyr::filter(average_year>floored_year) %>%
    select(average_year, index_year) %>%
    dplyr::mutate(index_year = 1:n())
  
  P_s_med <- array(dim=c(P_dims[3], length(time_index), P_dims[2])) # subnat, time, method
  
  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_s_med[j,k,r] <- sd(c(mod_P1[,r,j,time_index[k]] , mod_P2[,r,j,time_index[k]]))
      }
    }
  }
  
  P_s_med <- plyr::adply(P_s_med, c(1,2,3))
  colnames(P_s_med) <- c("index_subnat", "index_year", "index_method",  "sd_p")
  
  P_s_med <- P_s_med %>%
    group_by(index_method, index_subnat, index_year) %>%
    # summarise(median_p = stats::median(estimate)) %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(averageyear_index_table) %>%
    dplyr::mutate(Sector = sector_type)
  
  return(P_s_med)
}


# Calculation for sector data
P_median_global <- function(subnat_index_table, 
                            method_index_table, 
                            sector_type, 
                            year_index_table, 
                            file) { 
  Psamps <- readRDS(file)
  
  time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
    dplyr::filter(average_year>floored_year) %>% 
    dplyr::select(index_year) %>% 
    unlist() %>% 
    as.vector()
  
  averageyear_index_table <- year_index_table %>%
    dplyr::filter(average_year>floored_year) %>% 
    select(average_year, index_year) %>%
    dplyr::mutate(index_year = 1:n())
  
  P_dims <- dim(Psamps)

  #### P median
  P_s_med <- array(dim=c(length(time_index),P_dims[2],P_dims[3])) # year, method, subnat

  # Create a table for storing individual true country public data
  for(k in 1:length(time_index)) { # time loop
    for(s in 1:P_dims[3]) { # subnat
      for (m in 1:P_dims[2]) { # method
        P_s_med[k,m,s] <- stats::median(Psamps[,m,s,time_index[k]])
      }
    }
  }
  
  P_s_med <- plyr::adply(P_s_med, c(1,2,3))
  colnames(P_s_med) <- c("index_year", "index_method","index_subnat",  "median_p")
  
  P_s_med <- P_s_med %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(averageyear_index_table) %>%
    dplyr::mutate(Sector = sector_type)
  
  P_Q <- array(dim=c(length(time_index),P_dims[2],P_dims[3], 4)) # subnat, time, method, stats::quantile(95, 80)
  
  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_Q[k,r,j, 1:2] <- as.vector(unlist(stats::quantile(Psamps[,r,j,time_index[k]], prob=c(0.025, 0.975))))
        P_Q[k,r,j, 3:4] <- as.vector(unlist(stats::quantile(Psamps[,r,j,time_index[k]], prob=c(0.1, 0.9))))
      }
    }
  } 
  
  P_Q <- plyr::adply(P_Q, c(1,2,3))
  colnames(P_Q) <- c("index_year", "index_method", "index_subnat",  "lower_95", "uppper_95", "lower_80", "upper_80")
  
  P_s_med <- P_Q %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    dplyr::left_join(P_s_med)
  
  return(P_s_med)
}

# Calculation for sector data
P_quantiles_global <- function(subnat_index_table, 
                               method_index_table, 
                               sector_type, 
                               year_index_table, 
                               file_chain1, 
                               file_chain2) { 
  
  time_index <- year_index_table %>% # match index years to pooled years (pool 6 monthly estimates)
    dplyr::filter(average_year>floored_year) %>% 
    dplyr::select(index_year) %>% 
    unlist() %>% 
    as.vector()
  
  averageyear_index_table <- year_index_table %>%
    dplyr::filter(average_year>floored_year) %>% 
    select(average_year, index_year) %>%
    dplyr::mutate(index_year = 1:n())
  
  mod_P1 <- readRDS(file_chain1) 
  mod_P2 <- readRDS(file_chain2)
  
  P_dims <-  dim(mod_P1)
  # P_dims <- dim(Psamps)
  
  # #### P median
  # P_s_med <- array(dim=c(length(time_index),P_dims[2],P_dims[3],4)) # method, year, subnat, stats::quantile(95, 80)
  # 
  # # Create a table for storing individual true country public data 
  # for(k in 1:length(time_index)) { # time loop
  #   for(s in 1:P_dims[3]) { # subnat
  #     for (m in 1:P_dims[2]) { # method
  #       P_s_med[k,m,s,1:2] <- as.vector(unlist(stats::quantile(Psamps[,m,s,time_index[k]], prob=c(0.025, 0.975))))
  #       P_s_med[k,m,s,3:4] <- as.vector(unlist(stats::quantile(Psamps[,m,s,time_index[k]], prob=c(0.1, 0.9))))
  #     } 
  #   }
  # }
  
  P_s_med <- array(dim=c(P_dims[3], length(time_index), P_dims[2], 4)) # subnat, time, method, stats::quantile(95, 80)
  
  # Create a table for storing individual true country public data
  for(k in  1:length(time_index)) { # time loop
    for(j in 1:P_dims[3]) { # subnat
      for (r in 1:P_dims[2]) { # method
        P_s_med[j,k,r, 1:2] <- as.vector(unlist(stats::quantile(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]), prob=c(0.025, 0.975))))
        P_s_med[j,k,r, 3:4] <- as.vector(unlist(stats::quantile(c(mod_P1[,r,j,time_index[k]], mod_P2[,r,j,time_index[k]]), prob=c(0.1, 0.9))))
      }
    }
  } 
  
  P_s_med <- plyr::adply(P_s_med, c(1,2,3))
  colnames(P_s_med) <- c("index_subnat", "index_year", "index_method",  "lower_95", "uppper_95", "lower_80", "upper_80")
  
  P_s_med <- P_s_med %>%
    group_by(index_method, index_subnat, index_year) %>%
    # summarise(median_p = stats::median(estimate)) %>%
    dplyr::mutate(index_year = as.numeric(index_year)) %>%
    dplyr::mutate(index_method = as.numeric(index_method)) %>%
    dplyr::mutate(index_subnat = as.numeric(index_subnat)) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(averageyear_index_table) %>%
    dplyr::mutate(Sector = sector_type)
  
  return(P_s_med)
}

# Connect islands for adjancy matrix GIS:  https://gis.stackexchange.com/questions/413159/how-to-assign-a-neighbour-status-to-unlinked-polygons 
mstconnect <- function(polys, nb, distance="centroid", style="M"){
  
  if(distance == "centroid"){
    coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(polys)))
    dmat = as.matrix(dist(coords))
  }else if(distance == "polygon"){
    dmat = sf::st_distance(polys) + 1 # offset for adjacencies
    diag(dmat) = 0 # no self-intersections
  }else{
    stop("Unknown distance method")
  }
  
  gfull = igraph::graph.adjacency(dmat, weighted=TRUE, mode="undirected")
  gmst = igraph::mst(gfull)
  edgemat = as.matrix(igraph::as_adj(gmst))
  edgelistw = spdep::mat2listw(edgemat, style=style)
  edgenb = edgelistw$neighbour
  attr(edgenb,"region.id") = attr(nb, "region.id")
  allnb = spdep::union.nb(nb, edgenb)
  allnb
}

# Get alpha intercepts from JAGS file 
get_alpha <- function(mymod, modtype, subnat_index_table, method_index_table, sector_index_table) {
  subnatalpha_samps <-  mymod$BUGSoutput$sims.matrix[, which(grepl("alpha_snms", colnames(mymod$BUGSoutput$sims.matrix)))]
  subnatalpha_samps <- as_tibble(subnatalpha_samps) %>% pivot_longer(cols = everything(), names_to = "col_name", values_to = "estimate")
  subnatalpha_samps$col_details <- gsub(".*[[]([^.]+)[]].*", "\\1", subnatalpha_samps$col_name) # extract indexing
  
  subnatalpha_samps <- subnatalpha_samps %>%
    dplyr::rowwise() %>%
    dplyr::mutate(index_sector = as.numeric(unlist(strsplit(col_details, ","))[1]),
           index_method = as.numeric(unlist(strsplit(col_details, ","))[2]),
           index_subnat = as.numeric(unlist(strsplit(col_details, ","))[3])) %>%
    dplyr::left_join(subnat_index_table) %>%
    dplyr::left_join(method_index_table) %>%
    dplyr::left_join(sector_index_table) %>%
    dplyr::select(Country, Region, Method, Sector, estimate) %>%
    dplyr::mutate(Model = modtype, component = "Posterior")
  
  subnatalpha_samps_median <- subnatalpha_samps %>%
    dplyr::group_by(Sector, Method, Region) %>%
    dplyr::summarise(median_alpha = stats::median(estimate, na.rm = TRUE)) 
  
  subnatalpha_samps_sd <- subnatalpha_samps %>%
    dplyr::group_by(Sector, Method, Region) %>%
    dplyr::summarise(sd_alpha = sd(estimate, na.rm = TRUE)) 
  
  subnatalpha_samps_95CI <- left_join(subnatalpha_samps_median, subnatalpha_samps_sd)  %>%
    dplyr::rowwise() %>%
    dplyr::mutate(upper_95 = median_alpha+2*sd_alpha, 
           lower_95 =  median_alpha-2*sd_alpha) %>%
    dplyr::mutate(Model = modtype)
  
  return(list( alpha_samps = subnatalpha_samps, alpha_summary = subnatalpha_samps_95CI))
}
