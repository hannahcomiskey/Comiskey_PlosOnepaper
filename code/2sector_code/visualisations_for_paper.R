library(tidyverse)
library(tidybayes)
library(bayesplot)
library(coda)
library(shinystan)
library(rjags)
library(sf)

# Source code --------------------------------------
source('code/stan_utility.R')
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

vispath = 'visualisations/MVN_alpha/'

# fit_sim model ----------------------------------------------------------------------------------------
fit_sim <- readRDS('results/JAGS/JAGS_model_MVN_NCP_kstar_SE.RDS')

P_samps_df <- readRDS('results/JAGS/P_samps/MVN_alpha_chains/Nov2025/JAGS_MVN_alpha_Psamps_df.RDS')

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))


# Get observed data ------------------------------------------------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(y.cap = pmin(pmax(Public, 0.005), 0.995),
         logit.Public = log(y.cap/(1-y.cap)),
         logit.Public.Var = ((1/(y.cap*(1-y.cap)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))
M =  length(n_method)

P_df<- FP_source_data_wide %>%
  select(Country, Region, Method, index_year, average_year, Public, Private) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

P_SEdf<- FP_source_data_wide %>%
  select(Country, Region, Method, index_year, average_year, Public.SE, Private.SE) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public.SE, Private.SE), names_to = 'Sector', values_to = 'SE') %>%
  mutate(Sector = str_replace(Sector, '.SE', ''))

P_df <- left_join(P_df, P_SEdf) %>%
  rowwise() %>%
  mutate(lower_95 = Observed - 2*SE,
         upper_95 = Observed + 2*SE)

# Get correlations -------------------------------------------------------------

corr_country <- cov2cor(solve(fit_sim$BUGSoutput$mean$inv.Sigma.alpha_cms))
rownames(corr_country) <- colnames(corr_country) <- n_method

corr_prov <- cov2cor(solve(fit_sim$BUGSoutput$mean$inv.Sigma.alpha_pms))
rownames(corr_prov) <- colnames(corr_prov) <- n_method


country_plot <- ggcorrplot::ggcorrplot(corr_country, 
                                       type = "upper",
                                       tl.cex = 18,
                                       lab_size = 14,
                                       lab = TRUE)  + theme(legend.position = 'none')
prov_plot <- ggcorrplot::ggcorrplot(corr_prov,                                        
                                    type = "upper",
                                    tl.cex = 18,
                                    lab_size = 14,
                                    lab = TRUE)  + theme(legend.position = 'none')

corr_plots <- ggpubr::ggarrange(country_plot, prov_plot, labels = "AUTO")
ggsave(corr_plots, filename = paste0(vispath,'/parameter_plots/corr_plots.pdf'), height=12, width=15) 

# Coreelations with more conservative prior 

fit_cons <- readRDS('results/JAGS/JAGS_model_MVN_NCP_kstar_SE_higherDF_Wishart.RDS')

corr_country <- cov2cor(solve(fit_cons$BUGSoutput$mean$inv.Sigma.alpha_cms))
rownames(corr_country) <- colnames(corr_country) <- n_method

corr_prov <- cov2cor(solve(fit_cons$BUGSoutput$mean$inv.Sigma.alpha_pms))
rownames(corr_prov) <- colnames(corr_prov) <- n_method


country_plot <- ggcorrplot::ggcorrplot(corr_country, 
                                       type = "upper",
                                       tl.cex = 18,
                                       lab_size = 14,
                                       lab = TRUE)  + theme(legend.position = 'none')
prov_plot <- ggcorrplot::ggcorrplot(corr_prov,                                        
                                    type = "upper",
                                    tl.cex = 18,
                                    lab_size = 14,
                                    lab = TRUE)  + theme(legend.position = 'none')

corr_plots <- ggpubr::ggarrange(country_plot, prov_plot, labels = "AUTO")
ggsave(corr_plots, filename = paste0(vispath,'/parameter_plots/corr_plots_conservative.pdf'), height=12, width=15) 


# Subnational maps -------------------------------------------------------------
# load shapefiles

nigeria_preds <- P_samps_df %>% filter(Country == "Nigeria" & average_year %in% c(1990, 2000, 2010, 2020, 2025))
nigeria_shp <- read_sf("data/shape/dhs_ipumsi_ng/dhs_ipumsi_ng.shp") %>%
  mutate(ADMIN_NAME = str_trim(str_replace_all(ADMIN_NAME, "\r|\n", ""))) %>%
  filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON"))

nigeria_shp$Region <- dplyr::case_when(
  nigeria_shp$ADMIN_NAME %in% c("Federal Capital Territory Abuja", "Lagos", "Ogun", "Oyo", "Osun", "Ondo", "Ekiti") ~ "South-West",
  nigeria_shp$ADMIN_NAME %in% c("Edo", "Delta", "Bayelsa", "Rivers", "Akwa Ibom", "Cross River") ~ "South South",
  nigeria_shp$ADMIN_NAME %in% c("Abia", "Anambra", "Ebonyi", "Enugu", "Imo") ~ "South-East",
  nigeria_shp$ADMIN_NAME %in% c("Benue", "Kogi", "Kwara", "Nasarawa", "Niger", "Plateau", "FCT") ~ "North Central",
  nigeria_shp$ADMIN_NAME %in% c("Adamawa", "Bauchi", "Borno", "Gombe", "Taraba", "Yobe") ~ "North-East",
  nigeria_shp$ADMIN_NAME %in% c("Jigawa", "Kaduna", "Kano", "Katsina", "Kebbi", "Sokoto", "Zamfara") ~ "North-West"
)

table(nigeria_shp$ADMIN_NAME)
table(nigeria_shp$Region)

nga_zones <- nigeria_shp %>%
  st_make_valid() %>%          # ensure geometries are valid
  group_by(Region) %>%
  reframe(geometry = st_union(geometry)) %>%  
  st_as_sf() %>%  
  st_cast("MULTIPOLYGON") 

unique(nga_zones$Region)
st_geometry_type(nga_zones)  # should all be POLYGON or MULTIPOLYGON
nrow(nga_zones) 

nigeria_preds <- nigeria_preds %>% left_join(nga_zones)
table(nigeria_preds$Region)


nigeria_preds$Region <- dplyr::case_when(
  nigeria_preds$Region %in% c("Central", "North Central") ~ "N.Central",
  nigeria_preds$Region %in% c("North-East") ~ "N.East",
  nigeria_preds$Region %in% c("North-West") ~ "N.West",
  nigeria_preds$Region %in% c("South South") ~ "S.South",
  nigeria_preds$Region %in% c("South-East") ~ "S.East",
  nigeria_preds$Region %in% c("South-West") ~ "S.West"
)


ggplot() +
  geom_sf(data = nigeria_preds %>%   filter(Sector == "Private", Method %in% c("OC Pills", "IUD")), aes(geometry = geometry, fill = Mean)) +
  scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 1)) +
  geom_sf_text(
    data = nigeria_preds %>%   filter(Sector == "Private", Method %in% c("OC Pills", "IUD")),
    aes(geometry = geometry, label = Region),
    size = 3, color = "black", position = position_jitter(width = 0.3, height = 0.3)) +
  theme_minimal() +
  labs(x="", y="")+
  theme(legend.position = 'bottom') +
  facet_grid(Method~average_year)
ggsave(filename = paste0(vispath,'nigeria_maps.pdf'), height=12, width=15) 



rwanda_preds <- P_samps_df %>% filter(Country == "Rwanda" & average_year %in% c(1990, 2000, 2010, 2020, 2025))
rwanda_shp <- read_sf("data/shape/dhs_ipumsi_rw/dhs_ipumsi_rw.shp") %>%
  rename(Region = ADMIN_NAME)
  
rwanda_preds <- rwanda_preds %>% left_join(rwanda_shp)

ggplot() +
  geom_sf(data = rwanda_preds %>% filter(Sector=="Private"), aes(geometry = geometry, fill = Mean)) +
  scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 1)) +
  geom_sf_text(data = rwanda_preds %>% filter(Sector=="Private"), aes(geometry = geometry, label = Region), size = 3, color = "black") +
  theme_minimal() +
  labs(x="", y="")+
  theme(legend.position = 'bottom') +
  facet_grid(Method~average_year)
ggsave(filename = paste0(vispath,'rwanda_maps.pdf'), height=12, width=15) 



# Plot predictions ----------------------------------------------------------------
rwanda_obs <- P_df %>% filter(Country=="Rwanda" & Method %in% c("Implants", "Injectables", "IUD", "OC Pills"))
rwanda_preds <- P_samps_df %>% filter(Country=="Rwanda" & Method %in% c("Implants", "Injectables", "IUD", "OC Pills"))

ggplot() +
  geom_point(data = rwanda_obs, aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_errorbar(data = rwanda_obs, aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
  geom_line(data = rwanda_preds, aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = rwanda_preds, aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=8), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
  facet_grid(Region ~ Method)
ggsave(filename = paste0("rwanda_supplyshares_plot.pdf"), path = paste0(vispath, '/figs'), height=12, width=15) 

nga_obs <- P_df %>% filter(Country=="Nigeria" & Method %in% c("Female Sterilization", "Injectables", "IUD"))
nga_preds <- P_samps_df %>% filter(Country=="Nigeria" & Method %in% c("Female Sterilization", "Injectables", "IUD"))

ggplot() +
  geom_point(data = nga_obs, aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_errorbar(data = nga_obs, aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
  geom_line(data = nga_preds, aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = nga_preds, aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=8), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
  facet_grid(Region ~ Method)
ggsave(filename = paste0("nigeria_supplyshares_plot.pdf"), path = paste0(vispath, '/figs'), height=12, width=15) 


# Get survey data summary ----------------------------
n_province <- FP_source_data_long %>% 
  group_by(Country) %>%
  distinct(Region) %>%
  summarise(n_province = n())

n_survey <- FP_source_data_long %>%  
  group_by(Country) %>%
  distinct(average_year) %>%
  summarise(n_survey = n())
