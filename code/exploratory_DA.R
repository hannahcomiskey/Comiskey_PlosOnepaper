library(mcmsupply)
library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
library(ggcorrplot)

set.seed(1209)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

main_path = "results/global_subnat/"

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE", "Commercial_medical", "Commercial_medical.SE")]

most_recent_data <- FP_source_data_wide %>%
  rowwise() %>%
  mutate(CM_ratio = Commercial_medical/(1-Public),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.CM_ratio = log(CM_ratio/(1-CM_ratio))) %>%
  group_by(Country, Method, Region) %>%
  filter(average_year==max(average_year))

# Make the histogram
most_recent_data %>%
  ggplot(aes(x= logit.CM_ratio)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
  theme_bw()+
  ggtitle('Density of the most recently observed logit private ratio observations') +
  facet_wrap(~ Method) 
#ggsave('visualisations/global_subnat/density_rawdata_CMratio_bymethod.pdf')

all_data <- FP_source_data_wide %>%
  rowwise() %>%
  mutate(CM_ratio = Commercial_medical/(1-Public),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.CM_ratio = log(CM_ratio/(1-CM_ratio))) %>%
  group_by(Country, Method, Region) %>%
  select(Country, Region, Method, average_year, logit.CM_ratio)

n_method <- unique(most_recent_data$Method)

# Make the histogram
for(m in n_method) {
  ggplot(data=all_data %>% filter(Method==m), aes(x= logit.CM_ratio)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
    theme_bw()+
    ggtitle(paste0('Density of the observed logit private ratio ', m,' observations')) +
    facet_wrap(~average_year, scales='free_y')
  ggsave(paste0('visualisations/global_subnat/density_rawdata_CMratio_',m,'.pdf'))
}

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var),
         logit.CM.Var = ((1/(Commercial_medical*(1-Commercial_medical)))^2)*Commercial_medical.SE^2,
         logit.CM.SE = sqrt(logit.CM.Var))


# exploring the correlations for the logit.CM ratio in the most recent year

wide_test <- most_recent_data  %>% 
  select(Country, Region, Method, logit.CM_ratio) %>%
  distinct() %>%
  pivot_wider(names_from=Method, values_from = logit.CM_ratio) 

# Implants vs. Injectable
wide_test %>%
  ggplot(aes(x= Implants, y=Injectables)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_abline(slope=1, col='red')

wide_test %>%
  ungroup() %>%
  select(Implants, Injectables) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Implants, Injectables, method = "sp"))

# Implants vs IUD
wide_test %>%
  ggplot(aes(x= Implants, y=IUD)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_abline(slope=1, col='red')

wide_test %>%
  ungroup() %>%
  select(Implants, IUD) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Implants, IUD, method = "sp"))

# Implants vs. OC pills
wide_test %>%
  ggplot(aes(x= Implants, y=`OC Pills`)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_smooth()

wide_test %>%
  ungroup() %>%
  select(Implants, `OC Pills`) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Implants, `OC Pills`, method = "sp"))

# Implants vs. Female sterilization
wide_test %>%
  ggplot(aes(x= Implants, y=`Female Sterilization`)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_smooth() +
  xlim(5,7)

wide_test %>%
  ungroup() %>%
  select(Implants, `Female Sterilization`) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Implants, `Female Sterilization`, method = "sp"))

# Injectables vs OC pills
wide_test %>%
  ggplot(aes(x= Injectables, y=`OC Pills`)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_smooth()

wide_test %>%
  ungroup() %>%
  select(Injectables, `OC Pills`) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Injectables, `OC Pills`, method = "sp"))

# Injectables vs. IUD
wide_test %>%
  ggplot(aes(x= Injectables, y=IUD)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_smooth()

wide_test %>%
  ungroup() %>%
  select(Injectables, IUD) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Injectables, IUD, method = "sp"))

# Injectables vs Female sterilization
wide_test %>%
  ggplot(aes(x= Injectables, y=`Female Sterilization`)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_smooth()

wide_test %>%
  ungroup() %>%
  select(Injectables, `Female Sterilization`) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(Injectables, `Female Sterilization`, method = "sp"))

# OC pills vs Female sterilization
wide_test %>%
  ggplot(aes(x= `OC Pills`, y=`Female Sterilization`)) +
  geom_point() + 
  theme_bw()+
  ggtitle('Scatterplot of logit(private ratio) most recent observations') +
  geom_smooth()

wide_test %>%
  ungroup() %>%
  select(`OC Pills`, `Female Sterilization`) %>%
  filter(if_any(everything(), is.na)==FALSE) %>%
  summarize(corr = cor(`OC Pills`, `Female Sterilization`, method = "sp"))

# Correlation plot over all combinations
model.matrix(~0+., data=wide_test[,3:7]) %>% 
  cor(use="pairwise.complete.obs", method='kendall') %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2)


# Fishers transformation of correlation 
n_complete=8
sd.Z = (1/sqrt(n_complete-3))
z <- rnorm(10000,0, sd.Z)
plot(density(z))

rho.z <- ((exp(z)-1)/(exp(z)+1)) 
range(rho.z)
plot(density(rho.z))

# Rolling window of covariance -------------------------------------------------
# From the data 
public_covar_wide <- FP_source_data_wide %>%
  rowwise() %>%
  mutate(CM_ratio = Commercial_medical/(1-Public),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.CM_ratio = log(CM_ratio/(1-CM_ratio)),
         logit.public = log(Public/(1-Public))) %>%
  group_by(Country, Method, Region) %>%
  arrange(average_year) %>%
  select(Method, Country, Region, average_year, logit.public) %>%
  distinct() %>%
  pivot_wider(names_from=average_year, values_from= logit.public) 
public_covar_wide <- public_covar_wide[,4:31]

cov_public <- matrix(nrow=14, ncol=14)
for(i in 1:13) {
  for(j in (i+1):14) {
    cov_public[i,j] <- cov(public_covar_wide[,i],public_covar_wide[,j], use="na.or.complete")
  }
}

CM_covar_wide <- FP_source_data_wide %>%
  rowwise() %>%
  mutate(CM_ratio = Commercial_medical/(1-Public),
         logit.CM = log(Commercial_medical/(1-Commercial_medical)),
         logit.CM_ratio = log(CM_ratio/(1-CM_ratio)),
         logit.public = log(Public/(1-Public))) %>%
  group_by(Country, Method, Region) %>%
  arrange(average_year) %>%
  select(Method, Country, Region, average_year, logit.CM_ratio) %>%
  distinct() %>%
  pivot_wider(names_from=average_year, values_from= logit.CM_ratio)
CM_covar_wide <- CM_covar_wide[,4:17]

cov_CM <- matrix(nrow=14, ncol=14)
for(i in 1:13) {
  for(j in (i+1):14) {
    cov_CM[i,j] <- cov(CM_covar_wide[,i],CM_covar_wide[,j], use="na.or.complete")
  }
}


# From the model - The covariances are practically zero. This is due to the near-0 correlations used in the calculation of this parameter.
global_subnat_betas <- readRDS("~/Downloads/global_subnat_betas.RDS")
delta.k <- array(dim=c(2,5,173,12))
for(s in 1:2) {
  for(m in 1:5) {
    for(p in 1:173) { 
      for(h in 1:12) {
        delta.k[s,m,p,h] <- global_subnat_betas[s,m,p,h+1] - global_subnat_betas[s,m,p,h]
      }
    }
  }
}

delta.k_df <- plyr::adply(delta.k, c(1,2,3))
colnames(delta.k_df) <- c('Sector', 'Method', 'index_region', 1:12)

public_delta.k <- delta.k_df %>% filter(Sector==1)

# Knot 1
knot_data <- public_delta.k %>%
  select(!c(`2`:`12`)) %>%
  pivot_wider(names_from='Method', values_from = `1`) %>%
  select(!c(Sector, index_region))

cov <- matrix(nrow=5, ncol=5)
for(i in 1:4) {
  for(j in (i+1):5) {
    cov[i,j] <- cov(knot_data[,i], knot_data[,j])
  }
}

# EDA 
test <- FP_source_data_wide %>% 
  ungroup() %>% 
  select(Country, average_year) %>% 
  group_by(Country) %>% 
  filter(average_year==max(average_year)) %>% 
  distinct() %>% 
  ungroup() %>% 
  count(average_year)



