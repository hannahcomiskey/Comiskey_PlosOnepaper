library(tidyverse)

# Load data  -------------------------------------------------------------------
source('code/stan_utility.R')
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# Load results  ----------------------------------------------------------------
mod <- readRDS('results/JAGS/JAGS_model_MVN_NCP_kstar_SE.RDS') # no India, No AFG.

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))
subnat_index_table <- FP_source_data_wide %>% select(Country, Region, index_country, index_subnat) %>% distinct()

vars <- as.vector(unlist(dimnames(mod$BUGSoutput$sims.array)[3]))

subnat_index_table %>% tail()

grep("alpha_pms\\[5,165\\]", vars)
grep("beta.k\\[5,165,13\\]", vars)

# Get alpha_pms 
alpha_pms <- mod$BUGSoutput$sims.array[,,grep("alpha_pms\\[1,1\\]", vars)[1]:grep("alpha_pms\\[5,165\\]", vars)[1]]
alpha_pms <- rbind(alpha_pms[,1,], alpha_pms[,2,])
dim(alpha_pms)
beta.k <- mod$BUGSoutput$sims.array[,,grep("beta.k\\[1,1,1\\]", vars):grep("beta.k\\[5,165,13\\]", vars)]
beta.k <- rbind(beta.k[,1,], beta.k[,2,])
dim(beta.k)

# Calculate logit proportions --------------------------------------------------
n_samps=4000
z <- array(NA, dim=c(n_samps, length(n_method), length(n_subnat), n_years))
P <- array(NA, dim=c(n_samps, 2, length(n_method), length(n_subnat), n_years))
for(p in 1:length(n_subnat)) { # province loop matched to C
  for(m in 1:length(n_method)){ # method loop
    alpha_pms_samp <- alpha_pms[,grep(paste0('alpha_pms\\[',m,',',p, '\\]'), colnames(alpha_pms))]
    beta.k_samp <- beta.k[,grep(paste0('beta.k\\[',m,',',p,','), colnames(beta.k))]
    for (t in 1:n_years) { 
      for(s in 1:n_samps) {
        z[s,m,p,t] <- alpha_pms_samp[s] + B.ik[p,t,]%*%beta.k_samp[s,1:13]
        P[s,1,m,p,t] <- 1/(1+exp(-(z[s,m,p,t]))) 
        P[s,2,m,p,t] <- 1-P[s,1,m,p,t] 
      }
    } # end time loop
  } # end M loop 
} # end P loop

saveRDS(z, 'results/JAGS/P_samps/MVN_delta_chains/Nov2025/z_samps.RDS')
saveRDS(P, 'results/JAGS/P_samps/MVN_delta_chains/Nov2025/P_samps.RDS')

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))


# Get P samples for original model
P_samps <- readRDS('results/JAGS/P_samps/MVN_delta_chains/Nov2025/P_samps.RDS')
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(2,3,4))
colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
dim(P_samps.quantile)
colnames(P_samps.quantile) <- c('index_sector', 'index_method', 'index_subnat', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) %>%
  left_join(index_country_subnat_tbl) %>%
  left_join(year_index_table)

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
saveRDS(P_samps_df, file='results/JAGS/P_samps/MVN_delta_chains/Nov2025/JAGS_MVN_delta_Psamps_df.RDS')

