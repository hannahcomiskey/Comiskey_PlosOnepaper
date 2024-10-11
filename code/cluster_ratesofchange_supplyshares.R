set.seed(123)
library(mclust)
library(factoextra)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

# Source code --------------------------------------
source('code/load_functions.R')
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

# get country-region combos ---------------
country_regions_combos <- FP_source_data_wide %>%
  select(Country, Region) %>%
  distinct() %>%
  mutate(Region = str_to_sentence(Region))

# Join the subcontinental data ------------

area_classification <- mcmsupply::Country_and_area_classification %>% 
  dplyr::select(`Country or area`, Region) %>% 
  dplyr::rename(Country = `Country or area`) %>%
  dplyr::rename(Super_region = Region) %>%
  dplyr::mutate(Country = case_when(Country== "Bolivia (Plurinational State of)" ~ "Bolivia",
                                    Country== "Republic of Moldova" ~ "Moldova",
                                    Country== "Viet Nam" ~ "Vietnam",
                                    Country== "Kyrgyzstan" ~ "Kyrgyz Republic",
                                    .default = as.character(Country)))

head(area_classification)

# Read in wealth quintiles 
wq_data <- readRDS('data/IPUMS/SEdf/wealth_quintiles_ipums.RDS') %>%
  rename(wealth_proportion = Subnat_Freq) %>%
  pivot_wider(names_from = Wealth_quintiles, values_from = wealth_proportion) %>%
  mutate_all(~ifelse(is.nan(.), 0, .)) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>% 
  rowwise() %>% 
  mutate(check_sum = sum(Poorest, Poorer, Middle, Richer, Richest, Missing)) %>%
  mutate(Missing = ifelse(check_sum==0, 1, Missing)) %>%
  mutate(average_year = as.numeric(Year) + 0.5) %>%
  select(!Year)

religion <- readRDS('data/IPUMS/SEdf/religion_ipums.RDS') %>%
  pivot_wider(names_from = 'Religion', values_from = 'Subnat_Freq') %>%
  rowwise() %>%
  mutate(Christianity = ifelse(CHRISTIAN==0 | is.na(CHRISTIAN)==TRUE, sum(c_across(7:29), na.rm=TRUE), CHRISTIAN)) %>%
  mutate(Christianity = case_when(Country=='Madagascar'~ sum(c_across(7:29), na.rm=TRUE), 
                                  .default = Christianity)) %>%
  mutate(Buddhism = ifelse(`BUDDHIST/NEO-BUDDHIST`==0 | is.na(`BUDDHIST/NEO-BUDDHIST`)==TRUE, sum(c_across(31)), `BUDDHIST/NEO-BUDDHIST`)) %>%
  mutate(Buddhism = case_when(Country=='Madagascar'~ sum(c_across(30:31), na.rm=TRUE),
                                 .default = Buddhism)) %>%
  rename(Hinduism = HINDU) %>%
  rename(Islam = MUSLIM) %>%
  rename(Judaism = JEWISH) %>%
  mutate(Traditional = ifelse(`TRADITIONAL/SPIRITUAL/ANIMIST`==0 | is.na(`TRADITIONAL/SPIRITUAL/ANIMIST`)==TRUE, sum(c_across(35:47)), `TRADITIONAL/SPIRITUAL/ANIMIST`)) %>%
  mutate(Traditional = case_when(Country=='Madagascar'~ sum(c_across(34:47), na.rm=TRUE), 
                                  .default = Traditional)) %>%
  rename(None = `NO RELIGION`) %>%
  rename(Other = OTHER) %>%
  select(Country, Region, Year, Christianity, Buddhism, Hinduism, Islam, Judaism, Traditional, None, Other, Missing) %>%
  rowwise() %>% 
  drop_na(Christianity:Missing) %>%
  rowwise() %>%
  mutate(check_sum = sum(c_across(Christianity:Missing), na.rm = T)) %>%
  mutate(Missing = case_when(check_sum<0.99 ~ sum(c(Missing,(1-check_sum)), na.rm=TRUE),
                             .default = Missing)) %>%
  mutate(check_sum = sum(c_across(Christianity:Missing), na.rm = T)) %>%
  group_by(Country, Region) %>%
  filter(Year==max(Year)) %>%
  select(!c(check_sum,Year))

set1_countries <- religion %>% 
  ungroup() %>% 
  select(Country, Region) %>% 
  distinct() %>%
  mutate(Region = str_to_sentence(Region))

missingcountries <- setdiff(country_regions_combos , set1_countries)

# Read in religious identities for missing countries (using country-level identities in absence of regional info)
religious_identities <- readxl::read_excel('data/religious_identities_data.xlsx')

religious_identities <- missingcountries %>%
  left_join(religious_identities)

updated_religion <- rbind(religion, religious_identities) %>%
  mutate(Region = stringr::str_to_title(Region))

# # Read in rural % data
# rural_data <- read.csv('data/P_Data_Extract_From_World_Development_Indicators/rural_data.csv')
# rural_data <- rural_data[,c(3, 5:36)]
# colnames(rural_data) <- c('Country', 1990:2021)
# rural_data <- rural_data %>%
#   mutate(Country = case_when( Country=='Congo, Dem. Rep.' ~ "Congo Democratic Republic",
#                               TRUE ~ as.character(Country))) %>%
#   pivot_longer(cols = `1990`:`2021`, names_to='average_year', values_to = 'rural_proportion') %>%
#   mutate(average_year = as.numeric(average_year)+0.5,
#          rural_proportion = as.numeric(rural_proportion)/100)
#
# head(rural_data)

# Set up rates of change -------------
lagged_data <- FP_source_data_wide %>%
  left_join(area_classification) %>%
  #left_join(rural_data) %>%
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

# check data
#lagged_data %>% filter(Country=='Kenya' & Method=='Implants') %>% View()

richest_data <- wq_data %>% select(Country, Region, average_year, Richest, Poorest)

# Create dataframe for clustering
cluster_df <- lagged_data %>% 
  select(!c(Commercial_medical, Other, Public, logit.public, logit.CM, logit.other, lag.logit.public, lag.logit.CM, lag.logit.other)) %>% 
  group_by(Country, Region, Method) %>%
  filter(average_year==max(average_year)) %>% 
  filter(is.na(rate_of_change.public)==FALSE | is.na(rate_of_change.CM)==FALSE | is.na(rate_of_change.other)==FALSE) #%>%
  #left_join(richest_data, relationship = "many-to-many") %>%
  #left_join(updated_religion, relationship = "many-to-many") 

glimpse(cluster_df)

plotting_data <- cluster_df %>%
  pivot_longer(cols = c(rate_of_change.public, rate_of_change.CM, rate_of_change.other), names_to = 'Rate_type', values_to = 'value')


ggplot(plotting_data) +
  geom_density(aes(x=value, fill=Rate_type), alpha=0.5) +
  facet_wrap(~Method, scales='free')

glimpse(plotting_data)

# plot density of rates of change 
ggplot(cluster_df) +
  geom_density(aes(x=rate_of_change.public, fill=Method), alpha=0.5)

ggplot(cluster_df) +
  geom_density(aes(x=Christianity, fill=Method), alpha=0.5) +
  facet_wrap(~Method)

ggplot(cluster_df) +
  geom_density(aes(x=rate_of_change.CM, fill=Method), alpha=0.5)

ggplot(cluster_df) +
  geom_density(aes(x=rate_of_change.other, fill=Method), alpha=0.5)


ggplot(test) +
  geom_density(aes(x=Proportion, fill = Religion), alpha=0.5) +
  facet_wrap(~interaction(cluster, Method), scales='free')


# Clustering -------------------------------------------------------------------
X <- cluster_df %>% # add in in wealth data
  ungroup() %>%
  select(rate_of_change.public:rate_of_change.other)

# cluster on subcontinental regions ------------
ESS <- matrix(data=NA, nrow=2, ncol=9)
for(i in 2:10) {
  km.res <- kmeans(X,i, nstart = 25)
  print(km.res)
  ESS[1,(i-1)] <- i
  ESS[2,(i-1)] <- km.res$betweenss/km.res$totss
  
  aggregate(X, by=list(cluster=km.res$cluster), mean)
  rownames(X) <- paste0(cluster_df$religion_category,'_',1:nrow(cluster_df))
  g = fviz_cluster(km.res, X, 
                   repel=TRUE,
                   labelsize = NULL,
                   ellipse.type = "norm",  
                   geom="point",
                   palette = "Set2", 
                   ggtheme = theme_minimal())
  ggsave(g, filename = paste0('cluster_plot_all_methods_',i,'.pdf'), path = 'visualisations/global_subnat/', height=12, width=15) 
  
}

plot(ESS[1,], ESS[2,],  type="b", col="red", lwd=5, pch=15, xlab="# clusters", ylab="BSS/TSS ratio")

# Proceed with clusters
km.res <- kmeans(X,6, nstart = 10000)
print(km.res)

aggregate(X, by=list(cluster=km.res$cluster), mean)

# Reassign clusters to investigate trends 
cluster_df <- as.data.frame(cluster_df) %>%
  mutate(cluster = km.res$cluster)

raw_data <- FP_source_data_wide %>%
  group_by(Country, Region, Method) %>%
  arrange(average_year) %>%
  mutate(ROC.prop.public = Public - lag(Public),
         ROC.prop.CM = Commercial_medical - lag(Commercial_medical),
         ROC.prop.other = Other - lag(Other)) %>%
  filter(average_year == max(average_year)) %>%
  select(!c(Commercial_medical.SE, Public.SE, Other.SE, n_Other, n_Commercial_medical, n_Public,
            DEFT, check_total, count_SE.NA,index_subnat, index_country, index_method, index_superregion, index_year))

res_df <- cluster_df %>% 
  left_join(raw_data, relationship = "many-to-many") %>%
  left_join(updated_religion, relationship = "many-to-many")

res_df %>% group_by(cluster) %>% 
  summarise(across(c(rate_of_change.public:Missing, ROC.prop.public:ROC.prop.other), \(x) mean(x, na.rm=TRUE))) %>% 
  #mutate(across(c(rate_of_change.public:Missing, ROC.prop.public:ROC.prop.other), \(x) round(x, digits=3))) %>%
  View()

fviz_cluster(km.res, X, 
             repel=TRUE,
             labelsize = NULL,
             ellipse.type = "norm",  
             geom="point",
             palette = "Set2", 
             ggtheme = theme_minimal())


test <- res_df %>% pivot_longer(cols=c(Christianity:Missing), names_to = 'Religion', values_to = 'Proportion') %>%
  filter(Proportion>0) %>% 
  select(Country, Region, Method, Religion, Proportion, cluster) %>% 
  filter(Religion %in% c('Christianity', 'Islam', 'Hinduism', 'Traditional')) %>% 
  distinct()
  

# plot density of rates of change by religion 
ggplot(test) +
  geom_density(aes(x=Proportion, fill = Religion), alpha=0.5) +
  facet_wrap(~cluster, scales='free')

ggplot(test) +
  geom_density(aes(x=Proportion, fill = Religion), alpha=0.5) +
  facet_wrap(~interaction(cluster, Method), scales='free')


