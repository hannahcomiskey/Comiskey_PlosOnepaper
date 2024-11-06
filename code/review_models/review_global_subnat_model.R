# Source code --------------------------------------
library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
#library(sf)
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

respath = "results/global_subnat/" 
filepath = "visualisation/global_subnat/"

mod <- readRDS(paste0(respath, "test_corr/sdZ/subnational_model_all_mvn_alpha_4country_sdZ.RDS"))


# Summary of data ----------------------------------

n_provinces <- FP_source_data_wide %>% 
  ungroup() %>% 
  select(Country, Region) %>%
  distinct() %>%
  group_by(Country) %>%
  summarise(n_province = n())


n_survey <- FP_source_data_wide %>% 
  ungroup() %>% 
  select(Country, average_year) %>%
  distinct() %>%
  group_by(Country) %>%
  summarise(n_survey = n(),
            max_survey = max(average_year))

summary_info <- merge(n_provinces, n_survey)

summary_info

# Get standard error summary ---------------------------
SEdat <- FP_source_data_wide %>%
  select(Country, Region, Method, average_year, Commercial_medical.SE, Public.SE, Other.SE) %>%
  pivot_longer(cols = c(Commercial_medical.SE, Public.SE, Other.SE),  names_to = "Sector.SE", values_to =  "SE")

range(SEdat$SE)
median(SEdat$SE)

SEdat %>%
  ungroup() %>%
  group_by(Method) %>%
  summarise(mean.SE = mean(SE))

SEdat %>%
  ungroup() %>%
  group_by(Method) %>%
  summarise(range.SE = range(SE))

# Get params for single-pop model -------------------------------------
trans_rho <- mod$BUGSoutput$median$rho.trans
rho_alpha_cms <- (exp(2*trans_rho)-1)/(exp(2*trans_rho)+1)

Sigma_alpha_cms_hat <- mod$BUGSoutput$median$Sigma.alpha_cms

alpha_cms_hat <- mod$BUGSoutput$median$alpha_cms
#saveRDS(alpha_cms_hat, file=paste0(respath,"alpha_cms_hat.RDS"))

tau_alpha_pms_hat <- mod$BUGSoutput$median$tau_alpha_pms
#saveRDS(tau_alpha_pms_hat, file=paste0(respath,"tau_alpha_pms_hat.RDS"))

sigma_delta_hat <- mod$BUGSoutput$median$inv.sigma_delta
#saveRDS(sigma_delta_hat, file=paste0(respath,"inv.sigma_delta_hat.RDS"))

# Creating index tables for reference ------------------------------------------------------------
subnat_index_table <- FP_source_data_wide %>% select(Country, Region, index_subnat) %>% ungroup() %>% distinct()
glimpse(subnat_index_table)

country_index_table <- tibble(Country = unique(FP_source_data_wide$Country), index_country = unique(FP_source_data_wide$index_country))
glimpse(country_index_table)

method_index_table <- tibble(Method = n_method, index_method = 1:5)
glimpse(method_index_table)

sector_index_table <- tibble(Sector = c("Public", "Commercial_medical", "Other"), index_sector = 1:3)

year_index_table <- tibble(average_year = all_years, 
                           index_year = 1:n_all_years, 
                           floored_year = floor(all_years))

# Load calculations ----------------------------------------------------------
p_df <- readRDS(paste0(respath,"global_subnational_estimates.RDS"))

# Recreating data used in model for plotting ----------------------------------------------------------
FP_source_data_long_SE <- FP_source_data_wide %>% 
  select(Country, Region, Method, average_year, Commercial_medical.SE, Public.SE, Other.SE) %>%
  rename(Commercial_medical = Commercial_medical.SE , Public = Public.SE, Other = Other.SE) %>%
  gather(sector_category, SE.proportion, Commercial_medical:Other, factor_key=TRUE) 

FP_source_data_long <- FP_source_data_wide %>% 
  select(Country, Region, Method, average_year, Commercial_medical, Public, Other) %>%
  gather(sector_category, proportion, Commercial_medical:Other, factor_key=TRUE)

FP_source_data_long <- merge(FP_source_data_long, FP_source_data_long_SE)

# Adding error limits to data using SE ---------------------------------------------------------------
FP_source_data_long <- FP_source_data_long %>%
  rename(Sector = sector_category) %>%
  mutate( prop_min = proportion - 2*SE.proportion,
          prop_max = proportion + 2*SE.proportion ) %>%
  mutate(prop_max = ifelse(prop_max > 1, 1, prop_max)) %>%
  mutate(prop_min = ifelse(prop_min < 0, 0, prop_min))

# Getting TPR data ---------------------------------------------------------------
TPR_subset <- FP_source_data_wide %>%
  rowwise() %>%
  ungroup() %>%
  group_by(Country, Method, Region) %>%
  mutate(Total_private = Commercial_medical+Other) %>%
  mutate(CM_TP.ratio = Commercial_medical/Total_private)

TPR_long <- TPR_subset %>% 
  select(Country, Region, Method, average_year, Commercial_medical, Public, Other, CM_TP.ratio) %>%
  gather(sector_category, proportion, c(Commercial_medical, Public, Other, CM_TP.ratio), factor_key=TRUE)

TPR_long <- left_join(TPR_long, FP_source_data_long_SE)

TPR_long <- TPR_long %>%
  mutate( prop_min = ifelse(is.na(SE.proportion)==FALSE,proportion - 2*SE.proportion,NA),
          prop_max = ifelse(is.na(SE.proportion)==FALSE,proportion + 2*SE.proportion,NA)) %>%
  mutate(prop_max = ifelse(prop_max > 1, 1, prop_max)) %>%
  mutate(prop_min = ifelse(prop_min < 0, 0, prop_min)) 

# Country estimates -----------------------------------------------------------------------------------
safe_colorblind_palette <- c("#56B4E9", "#E69F00", "#999999",  "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
for(i in 1:nrow(subnat_index_table)) {
  country_data <- FP_source_data_long %>% 
    filter(Country==subnat_index_table$Country[i] & Region==subnat_index_table$Region[i])
    
  country_calc <- p_df %>% filter(Country==subnat_index_table$Country[i] & Region==subnat_index_table$Region[i]) %>%
    rename(upper_95 = uppper_95)
  
  ci_plot = ggplot() +  # plot of true p value vs time using facet wrap
    geom_line(data=country_calc, aes(x=average_year, y=median_p, color=Sector)) +
    geom_point(data=country_data, aes(x=average_year, y=proportion, colour=Sector))+
    geom_errorbar(data=country_data, aes(ymin = prop_min, ymax = prop_max, x=average_year, colour=Sector), width = 1.5) +
    geom_ribbon(data=country_calc, aes(ymin = lower_95, ymax = upper_95, x=average_year, fill=Sector), alpha=0.2) +
    geom_ribbon(data=country_calc, aes(ymin = lower_80, ymax = upper_80, x=average_year, fill=Sector), alpha=0.26) +
    labs(y="Proportion of contraceptives supplied", x = "Year", title = paste0(unique(country_calc$Region),", ",unique(country_calc$Country))) +
    scale_y_continuous(limits=c(0,1))+
    theme_bw() +
    theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
    scale_colour_manual(breaks = c("Public", "Commercial_medical", "Other"), values=safe_colorblind_palette) +
    scale_fill_manual(breaks = c("Public", "Commercial_medical", "Other"), values=safe_colorblind_palette) +
    labs(fill = "Sector") +
    guides(color="none") + 
    facet_wrap(~Method)
  
  subset_C <- subnat_index_table$Country[i]
  region_name <- str_replace_all(subnat_index_table$Region[i], "[[:punct:]]", "_") # remove special characters
  region_name <- str_replace_all(region_name, " ", "") # remove spaces
  
  ggsave(ci_plot, filename = paste0(i,"_p_estimates.pdf"), path = paste0(filepath,"overleaf"), height=12, width=15) 
}

# Get alpha intercepts  ----------------------------------------------------------------------------------------------
alpha_pms <- mod$BUGSoutput$median$alpha_pms 
invlogit.alpha_pms <- exp(alpha_pms)/(1+exp(alpha_pms))
invlogit.alpha_pms <- plyr::adply(invlogit.alpha_pms, c(2,3))
colnames(invlogit.alpha_pms) <- c("index_method","index_subnat","Public", "CM_TP.ratio")
invlogit.alpha_pms <- as_tibble(invlogit.alpha_pms) %>%
  pivot_longer(names_to="Sector", values_to = "alpha_intercept", cols = c(Public,CM_TP.ratio)) %>%
  mutate(index_method = as.numeric(index_method)) %>%
  mutate(index_subnat = as.numeric(index_subnat)) %>%
  left_join(subnat_index_table) %>%
  left_join(method_index_table)

# Plotting estimates --------------------------------
med_P <- p_df %>%
  select(Country, Region, Method, Sector, median_p, average_year) %>%
  tidyr::pivot_wider(
    names_from  = c(Sector), # Can accommodate more variables, if needed.
    values_from = c(median_p)) %>%
  rowwise() %>%
  mutate(Total_private = sum(Commercial_medical, Other, na.rm = TRUE)) %>%
  mutate(CM_TP.ratio = Commercial_medical/Total_private)
med_P <- med_P %>% pivot_longer(!c(Country, Region, Method, average_year),names_to = "Sector", values_to = "median_p")

lower_P <- p_df %>%
  select(Country, Region, Method, Sector, lower_95, average_year) %>%
  tidyr::pivot_wider(
    names_from  = c(Sector), # Can accommodate more variables, if needed.
    values_from = c(lower_95)) %>%
  rowwise() %>%
  mutate(Total_private = sum(Commercial_medical, Other, na.rm = TRUE)) %>%
  mutate(CM_TP.ratio = Commercial_medical/Total_private)
lower_P <- lower_P %>% pivot_longer(!c(Country, Region, Method, average_year),names_to = "Sector", values_to = "lower_95")

upper_P <- p_df %>%
  select(Country, Region, Method, Sector, upper_95, average_year) %>%
  tidyr::pivot_wider(
    names_from  = c(Sector), # Can accommodate more variables, if needed.
    values_from = c(upper_95)) %>%
  rowwise() %>%
  mutate(Total_private = sum(Commercial_medical, Other, na.rm = TRUE)) %>%
  mutate(CM_TP.ratio = Commercial_medical/Total_private)
upper_P <- upper_P %>% pivot_longer(!c(Country, Region, Method, average_year),names_to = "Sector", values_to = "upper_95")

TP_ratio_estimates <- left_join(med_P, lower_P)
TP_ratio_estimates <- left_join(TP_ratio_estimates, upper_P)

TP_ratio_estimates <- TP_ratio_estimates %>% filter(Sector %in% c("Public", "Commercial_medical", "Other", "CM_TP.ratio"))
TP_ratio_estimates <- standard_method_names(TP_ratio_estimates)

# Plot intercepts ------------------------
for(i in 1:nrow(subnat_index_table)) {
  print(i)
  country_data <- TPR_long %>% 
    filter(Country==subnat_index_table$Country[i] & Region==subnat_index_table$Region[i] & sector_category %in% c("Public","CM_TP.ratio")) %>%
    rename(Sector = sector_category)
  TP_ratio_df <- TP_ratio_estimates %>% 
    filter(Country==subnat_index_table$Country[i] & Region==subnat_index_table$Region[i] & Sector %in% c("Public","CM_TP.ratio"))
  alpha_data <- invlogit.alpha_pms %>% filter(Country==subnat_index_table$Country[i] & Region==subnat_index_table$Region[i] & Sector %in% c("Public","CM_TP.ratio")) 
  
  region_name <- str_replace_all(subnat_index_table$Region[i], "[[:punct:]]", "_") # remove special characters
  region_name <- str_replace_all(region_name, " ", "") # remove spaces
  
  ci_plot = 
    ggplot() + # HC: Using my all_p_med table calculated above. Plot looks the same as yours calculated above.
    ggtitle(paste0("Contraceptive Supply Intercepts for ",region_name,",",subnat_index_table$Country[i])) +
    geom_point(data=country_data, aes(x=average_year, y=proportion, colour=Sector))+
    geom_hline(data = alpha_data, aes(yintercept = alpha_intercept, colour=Sector), linetype = "longdash") +
    labs(y="Proportion of contraceptives supplied", x = "Year") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(strip.text.x = element_text(size = 9)) +
    theme(legend.position = "bottom") +
    scale_colour_manual(values=safe_colorblind_palette) +
    scale_fill_manual(values=safe_colorblind_palette) +
    labs(fill = "Sector") +
    ylim(0,1) +
    facet_wrap(~Method)
  
  ggsave(ci_plot, filename = paste0("plotting_intercepts_",subnat_index_table$Country[i],"_",region_name,".pdf"), path = paste0(filepath,"/intercepts"), height=12, width=15) 
}

# Review implants ----------------------------------

implants_calc <- p_df %>% 
  filter(Sector=="Commercial_medical" & Method=="Implants") %>%
  mutate(Country_region = paste0(Country, "_", Region))
implants_data <- FP_source_data_long  %>% 
  filter(Sector=="Commercial_medical" & Method=="Implants") %>%
  mutate(Country_region = paste0(Country, "_", Region))

ggplot() +
  geom_line(data=implants_calc, aes(x=average_year, y=median_p, group=Country_region, colour=Country_region))+
  geom_point(data=implants_data, aes(x=average_year, y=proportion, group=Country_region, colour=Country_region)) +
  geom_errorbar(data=implants_data, aes(x=average_year, ymin=prop_min, ymax=prop_max, group=Country_region, colour=Country_region)) +
  theme(strip.text.x = element_text(size = 9), axis.text.x = element_text(angle = 90), legend.position = "none") +
  facet_wrap(~Country)
#ggsave(filename = "CM_sector_implants_global.pdf", path = filepath, height=12, width=15) 

## Look at sample size vs. direct:estimates ---------------------------------

preds <- p_df %>% select(Country, Region, Method, Sector, average_year, median_p)

FP_source_data_long_N <- FP_source_data_wide %>% 
  select(Country, Region, Method,  average_year, n_Commercial_medical, n_Public, n_Other) %>%
  rename(Commercial_medical = n_Commercial_medical , Public = n_Public, Other = n_Other) %>%
  gather(Sector, N, Commercial_medical:Other, factor_key=TRUE) 

FP_source_data_long <- merge(FP_source_data_long_N, FP_source_data_long)

data <- FP_source_data_long %>% 
  left_join(preds) %>%
  rowwise() %>%
  mutate(obs_pred = proportion/median_p)

ggplot() +
  geom_point(data = data, aes(x=N, y=obs_pred)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.text.x = element_text(size = 11), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
  labs(y="observed/modelled")
ggsave(filename = "sample_size_vs_estimate_ratio.pdf", path = filepath, width=15, height=12)

ggplot() +
  geom_point(data = data, aes(x=log(N), y=log(obs_pred))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.text.x = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
  labs(y="log(observed/modelled)")
ggsave(filename = "log_sample_size_vs_log_estimate_ratio.pdf", path = filepath, width=15, height=12)

# Compare variance of estimates to standard error of data

sd_p_samps <- readRDS(paste0(respath,"global_subnational_P_sd.RDS"))
sd_obs <- left_join(FP_source_data_long, sd_p_samps) %>%
  rowwise() %>%
  mutate(sd_diff = SE.proportion - sd_p)

ggplot() +
  geom_boxplot(data = sd_obs %>% filter(Country %in% c("Niger", "Rwanda")), aes(x=Method, y=sd_diff, colour=Country)) +
  theme_bw() +
  theme(strip.text.x = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20)) +
  theme(legend.position = "bottom", axis.text.x = element_text(size =20, angle = 90), strip.text.x = element_text(size=20)) +
  facet_wrap(~Sector) 

ggplot() +
  geom_point(data = sd_obs, aes(x=SE.proportion, y=sd_p, colour=Method)) +
  geom_abline(intercept=0, slope=1) +
  labs(x= "Observed SE", y="Estimated sd") +
  theme_bw() +
  theme(strip.text.x = element_text(size=20),  axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20)) +
  theme(legend.position = "bottom", axis.text.x = element_text(size =20, angle = 90), strip.text.x = element_text(size=20)) +
  facet_wrap(~Sector) 
ggsave(filename = "obsevedSE_estimateSD.pdf", path = filepath, width=15, height=12)


