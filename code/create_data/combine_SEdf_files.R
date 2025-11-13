# Source files
library(tidyverse)

# Read in results
my_SEestimate_files <- list.files("data/IPUMS/SEdf")

# remove dudlist
my_SEestimate_files <- my_SEestimate_files[my_SEestimate_files != "dud_SElist_jobrun.RDS"]
my_SEestimate_files <- my_SEestimate_files[my_SEestimate_files != "wealth_quintiles_ipums.RDS"]
my_SEestimate_files <- my_SEestimate_files[my_SEestimate_files != "religion_ipums.RDS"]

# Create list of file names
tmp <- strsplit(sub("(_)(?=[^_]+$)", " ", my_SEestimate_files, perl=T), " ")
file_names <- unlist(lapply(tmp, function(l) l[[1]]))

# Read in all data into list 
df.list <- lapply(paste0("data/IPUMS/SEdf/",my_SEestimate_files), function(i){
  x = readRDS(i)
  x
})

names(df.list) <- file_names

# Match country codes and country to create new column in each file
country_codes_DHS <- readxl::read_excel("data/country_codes_DHS.xlsx") # country codes and names

# Run function to combine dfs into one
country_SEestimates <- tibble()
for(i in 1:length(file_names)) {
  mydata <- tibble(df.list[[i]])
  mydata <- mydata %>%
    rename(proportion = Subnat_Freq) %>%
    rename(SE.proportion = Subnat_Freq_SE)
  country_SEestimates <- bind_rows(mydata, country_SEestimates)
}

# Checking original data to see if sums to 1
country_SEestimates <- country_SEestimates %>%
  group_by(Year, Method, Country, Region) %>%
  mutate(check_sum = sum(proportion, na.rm=TRUE))

# Checking what surveys are present 
View(country_SEestimates %>% select(Country, Year) %>% distinct())

# Extracting the method we are interested in for now (removing condoms, other, emergency)
country_SEestimates <- country_SEestimates %>%
  filter(Method %in% c("Norplant/Implants", "Injections", "IUD", "Female Sterilization", "Pill")) %>%
  arrange(Country, Year)

# Recoding methods to match original data
cleaned_SE_source <- country_SEestimates %>%
  mutate(Method = recode(Method, "Female Sterilization" = "Female Sterilization",
                         "Norplant/Implants" = "Implants",
                         "Pill" = "OC Pills",
                         "IUD" = "IUD",
                         "Injections" = "Injectables")) %>%
  rename(sector_category = sector_categories) %>%
  mutate(average_year = Year + 0.5) %>%
  mutate(Region = str_to_title(Region)) %>%
  mutate(Region = case_when(Country=='Burkina Faso' & Region =='Est' ~ 'East', 
                            Country=='Burkina Faso' & Region =='Sud' ~ 'South',
                            Country=='Burkina Faso' & Region =='Nord' ~ 'North',
                            Country=='Burkina Faso' & Region =='Centre-Ouest' ~ 'Central/West',
                            Country=='Burkina Faso' & Region =='Centre-Sud' ~ 'Central/South',
                            Country=='Burkina Faso' & Region =='Centre-Nord' ~ 'Central/North',
                            Country=='Burkina Faso' & Region =='Centre-Est' ~ 'Central/East',
                            Country=='Burkina Faso' & Region =='Sud-Ouest' ~ 'South/West',
                            Country=="Cote d'Ivoire" & Region =='South Without Abidjan' ~ 'South',
                            Country=="Cote d'Ivoire" & Region =='Center' ~ 'Centre',
                            Country=="Cote d'Ivoire" & Region =='Capital (Abidjan)' ~ 'Abidjan',
                            Country=="Cote d'Ivoire" & Region =='City Of Abidjan' ~ 'Abidjan',
                            Country=="Cote d'Ivoire" & Region =='Center East' ~ 'Centre-East',
                            Country=="Cote d'Ivoire" & Region =='Center North' ~ 'Centre-North',
                            Country=="Cote d'Ivoire" & Region =='Center West' ~ 'Centre-West',
                            Country=="Cote d'Ivoire" & Region =='North East' ~ 'North-East',
                            Country=="Cote d'Ivoire" & Region =='North West' ~ 'North-West',
                            Country=="Cote d'Ivoire" & Region =='South West' ~ 'South-West',
                            Country=="Nigeria" & Region =='Southeast' ~ 'South-East',
                            Country=="Nigeria" & Region =='South East' ~ 'South-East',
                            Country=="Nigeria" & Region =='Southwest' ~ 'South-West',
                            Country=="Nigeria" & Region =='South West' ~ 'South-West',
                            Country=="Nigeria" & Region =='Northwest' ~ 'North-West',
                            Country=="Nigeria" & Region =='North West' ~ 'North-West',
                            Country=="Nigeria" & Region =='Northeast' ~ 'North-East',
                            Country=="Nigeria" & Region =='North East' ~ 'North-East',
                            .default = as.character(Region))) %>%
  filter(Region!='Na' & Region != "Countryside" & Region != "Small City") 

# relgion <- readRDS('data/IPUMS/SEdf/religion_ipums.RDS') %>% 
#   rename(religion_proportion = Subnat_Freq) %>%
#   mutate(average_year = as.numeric(Year)+0.5) %>%
#   select(!c(Year)) %>%
#   filter(is.na(religion_proportion)==FALSE)
# 
# wealth <- readRDS('data/IPUMS/SEdf/wealth_quintiles_ipums.RDS') %>% 
#   rename(wealth_quintile_proportion = Subnat_Freq)  %>%
#   mutate(average_year = as.numeric(Year)+0.5) %>%
#   select(!c(Year)) %>%
#   pivot_wider(names_from = Wealth_quintiles, values_from = wealth_quintile_proportion)
# cleaned_SE_source <- cleaned_SE_source %>% left_join(relgion, relationship = "many-to-many") %>% left_join(wealth, relationship = "many-to-many")

saveRDS(cleaned_SE_source, file="data/subnat_SE_source_data.RDS")

# Removing observations with less than 15 sampling units
SE15_source <- cleaned_SE_source %>%
  filter(n >=15)
included_groups <- SE15_source %>% select(Country, Method, average_year) %>% distinct()
SE15_source1 <- left_join(included_groups, cleaned_SE_source)
saveRDS(SE15_source1, file="data/subnat_SE_source_data_15.RDS")

# Removing observations with less than 20 sampling units and keeping 3 sectors if at least one is >= 20.
SE20_source <- cleaned_SE_source %>%
  filter(n >=20 | is.na(n))
included_groups <- SE20_source %>% select(Country, Method, average_year) %>% distinct()
SE20_source1 <- left_join(included_groups, cleaned_SE_source)
saveRDS(SE20_source1, file="data/subnat_SE_source_data_20.RDS")
