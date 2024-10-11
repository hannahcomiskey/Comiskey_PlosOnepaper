# Install package from CRAN
library(scales)
library(countries)
library(ggmap)
library(geosphere)
source("code/read_in_subnational_SEdata.R")

n_country_names <- as_tibble(Country = n_country) %>%
  mutate(Country = case_when(Country=='Rwanda' ~ 'Rawanda',
                             Country=='Congo Democratic Republic' ~ 'Democratic Republic of the Congo',
                             Country=='Tanzania' ~ 'United  Republic of Tanzania',
                             .default = as.character(Country)))

capitals_data <- mapdeck::capitals %>% filter(country %in% n_country_names$Country) %>%
  mutate(Country = case_when(country=='Rawanda' ~ 'Rwanda',
                             country=='Democratic Republic of the Congo' ~ 'Congo Democratic Republic',
                             country=='United  Republic of Tanzania' ~ 'Tanzania',
                             .default = as.character(country))) %>%
  select(!country) 

distmat <- matrix(NA, nrow=nrow(capitals_data), ncol=nrow(capitals_data))
for(i in 1:nrow(distmat)) {
  d1 <- as.vector(unlist(capitals_data[i, c('lon', 'lat')]))
  for(j in 1:nrow(distmat)) {
    d2 <-  as.vector(unlist(capitals_data[j, c('lon', 'lat')]))
    distmat[i,j] <- distm(d1,d2, fun = distGeo)
  }
}

colnames(distmat) <- capitals_data$Country

distmat_df <- as_tibble(bind_cols(Country = capitals_data$Country, distmat))

norm_distmat <- apply(distmat, 1, rescale)