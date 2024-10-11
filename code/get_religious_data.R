# Source code --------------------------------------
library(tidyverse)
library(tidybayes)
source('code/load_functions.R')
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")


set.seed(1209)

country_regions_combos <- FP_source_data_wide %>%
  select(Country, Region) %>%
  distinct()

# read in data 
religion <- readRDS('data/IPUMS/SEdf/religion_ipums.RDS') %>%
  mutate(religion_category = case_when(Religion=="MUSLIM" ~ "Muslim",
                                       Religion=="NO RELIGION" ~ "None",
                                       Religion=="CHRISTIAN" ~ "Christian",
                                       Religion=="Catholic" ~ "Christian",
                                       Religion=="Orthodox" ~ "Christian",
                                       Religion=="Protestant" ~ "Christian",
                                       Religion=="Anglican" ~ "Christian",
                                       Religion=="Presbyterian" ~ "Christian",
                                       Religion=="Baptist/Seventh-day Adventist" ~ "Christian",
                                       Religion=="Baptist" ~ "Christian",
                                       Religion=="Seventh-day Adventist" ~ "Christian",
                                       Religion=="Apostolic" ~ "Christian",
                                       Religion=="Salvation Army" ~ "Christian",
                                       Religion=="Methodist" ~ "Christian",
                                       Religion=="Pentecostal-based" ~ "Christian",
                                       Religion=="Pentecostal" ~ "Christian",
                                       Religion=="Celestial Church of Christ" ~ "Christian",
                                       Religion=="Assemblies of God" ~ "Christian",
                                       Religion=="Charismatic" ~ "Christian",
                                       Religion=="Other Protestant" ~ "Christian",
                                       Religion=="Evangelical" ~ "Christian",
                                       Religion=="Kimbanguist (Congo, Democratic Republic and Republic)" ~ "Christian",
                                       Religion=="Christian Zionist" ~ "Christian",
                                       Religion=="Jehovah's Witness" ~ "Christian",
                                       Religion=="Other Christian, country-specific" ~ "Christian",
                                       Religion=="Mammon (Uganda)" ~ "Christian",
                                       Religion=="BUDDHIST/NEO-BUDDHIST" ~ "Buddhist",
                                       Religion=="Buddhist" ~ "Buddhist",
                                       Religion=="HINDU" ~ "Hindu",
                                       Religion=="JEWISH" ~ "Jewish",
                                       Religion=="TRADITIONAL/SPIRITUAL/ANIMIST" ~ "Other",
                                       Religion=="Traditional" ~ "Other",
                                       Religion=="Spiritual" ~ "Other",
                                       Religion=="Animist" ~ "Other",
                                       Religion=="Donyi-Polo" ~ "Other",
                                       Religion=="Sanamahi" ~ "Other",
                                       Religion=="Vodun" ~ "Other",
                                       Religion=="Baha'i" ~ "Other",
                                       Religion=="Sikh" ~ "Other",
                                       Religion=="Zoroastrian" ~ "Other",
                                       Religion=="Jain" ~ "Other",
                                       Religion=="Bundu dia Kongo (Congo, Democratic Republic)" ~ "Other",
                                       Religion=="Vuvamu (Congo, Democratic Republic)" ~ "Other",
                                       Religion=="Kirat Mundhum (Nepal)" ~ "Other",
                                       Religion=="OTHER" ~ "Other",
                                       Religion=="Missing" ~ 'Missing')) %>%
  group_by(Country, Region, Year) %>%
  rename(proportion = Subnat_Freq) %>%
  distinct() %>%
  filter(proportion >0)

# clean data
religion_1 <- religion %>%
  filter(Religion %in% c('CHRISTIAN', 'MUSLIM', "NO RELIGION", "BUDDHIST/NEO-BUDDHIST", "HINDU", "JEWISH", "TRADITIONAL/SPIRITUAL/ANIMIST", "OTHER")) %>%
  group_by(Country, Region, Year) %>%
  mutate(check_sum = sum(proportion, na.rm=TRUE)) %>%
  filter(check_sum==1)

religion_2 <- religion %>%
  group_by(Country, Region, Year, religion_category) %>%
  anti_join(religion_1) %>%
  mutate(proportion_category = sum(proportion, na.rm=TRUE)) %>%
  select(Country, Region, Year, religion_category, proportion_category) %>%
  distinct() %>%
  ungroup() %>%
  group_by(Country, Region, Year) %>%
  mutate(check_sum = sum(proportion_category))

R1 <- religion_1 %>%
  rename(proportion_category = proportion) %>%
  select(Country, Region, Year, religion_category, proportion_category)
R2 <- religion_2 %>%
  select(Country, Region, Year, religion_category, proportion_category)

religion_data <- rbind(R1, R2) %>%
  mutate(average_year = as.numeric(Year)+0.5) %>%
  ungroup() %>%
  select(!Year) %>%
  rowwise() %>%
  filter(is.na(religion_category)==FALSE) %>%
  pivot_wider(names_from = religion_category, values_from=proportion_category) %>%
  mutate(across(None:Missing, .fns = ~replace_na(.,0))) %>%
  rename(Other_religion = Other)

# Set up rates of change -------------
lagged_data <- FP_source_data_wide %>%
  left_join(area_classification) %>%
  group_by(Country, Region, Method) %>%
  arrange(average_year) %>%
  mutate(logit.public = log(Public/(1-Public)),
         lag.logit.public = lag(logit.public),
         rate_of_change.public = logit.public - lag.logit.public,
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         lag.logit.CM = lag(logit.CM),
         rate_of_change.CM = logit.CM - lag.logit.CM,
         logit.other = log(Other/(1-Other)),
         lag.logit.other = lag(logit.other),
         rate_of_change.other = logit.other - lag.logit.other) %>%
  select(!c(Commercial_medical.SE, Public.SE, Other.SE, n_Other, n_Commercial_medical, n_Public,
            DEFT, index_subnat, index_country, index_method, index_superregion, index_year,
            check_total, count_SE.NA))

# Create dataframe for clustering
cluster_df <- lagged_data %>% 
  select(!c(Commercial_medical, Other, Public, logit.public, logit.CM, logit.other, lag.logit.public, lag.logit.CM, lag.logit.other)) %>% 
  group_by(Country, Region, Method) %>%
  filter(average_year==max(average_year)) %>% 
  filter(is.na(rate_of_change.public)==FALSE | is.na(rate_of_change.CM)==FALSE | is.na(rate_of_change.other)==FALSE) %>%
  left_join(religion_data) 

cluster_df <- cluster_df[which(complete.cases(cluster_df)==TRUE),]

X <- cluster_df %>% # add in in wealth data
  ungroup() %>%
  select(rate_of_change.public:Missing) %>% # Poorest) %>%
  scale() 

# Proceed with clusters
km.res <- kmeans(X,6, nstart = 10000)
print(km.res)

aggregate(X, by=list(cluster=km.res$cluster), mean)

# Reassign clusters to investigate trends 
info_df <- as.data.frame(cluster_df) %>%
  mutate(cluster = km.res$cluster)  %>%
  select(Super_region, Country, Region, Method, average_year, cluster)

data <- FP_source_data_wide %>%
  group_by(Country, Region, Method) %>%
  arrange(average_year) %>%
  select(!c(Commercial_medical.SE, Public.SE, Other.SE, n_Other, n_Commercial_medical, n_Public,
            DEFT, check_total, count_SE.NA)) 

res_df <- info_df %>% 
  left_join(data, relationship = "many-to-many") %>%
  left_join(updated_religion, relationship = "many-to-many") %>%
  filter(is.na(Commercial_medical)==FALSE & is.na(Public)==FALSE & is.na(Other)==FALSE)




