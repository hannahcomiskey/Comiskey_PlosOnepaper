library(tidyverse)

# Load data  -------------------------------------------------------------------
source('code/stan_utility.R')
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# Load results  ----------------------------------------------------------------
mod <- readRDS('results/JAGS/JAGS_mod_Bspline_N_TG_all.RDS')

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))
subnat_index_table <- FP_source_data_wide %>% select(Country, Region, index_country, index_subnat) %>% distinct()

vars <- as.vector(unlist(dimnames(mod$BUGSoutput$sims.array)[3]))

subnat_index_table %>% tail()

vars[grep("a0\\[1,1\\]", vars)]

grep("a\\[5,195,13\\]", vars)


# a0[m, p] + inprod(B[p, t, 1:num_knots], a[m, p, 1:num_knots])

# Get alpha_pms
alpha_pms <- mod$BUGSoutput$sims.array[,,grep("a0\\[1,1\\]", vars)[1]:grep("a0\\[5,195\\]", vars)[1]]
alpha_pms <- rbind(alpha_pms[,1,], alpha_pms[,2,], alpha_pms[,3,])
dim(alpha_pms)
beta.k <- mod$BUGSoutput$sims.array[,,grep("a\\[1,1,1\\]", vars)[1]:grep("a\\[5,195,13\\]", vars)[1]]
beta.k <- rbind(beta.k[,1,], beta.k[,2,], beta.k[,3,])
dim(beta.k)

# Calculate logit proportions --------------------------------------------------
n_samps=2000
n_chains=3
for(c in 1:n_chains){
  print(c)
  z <- array(NA, dim=c(n_samps, length(n_method), length(n_subnat), n_years))
  P <- array(NA, dim=c(n_samps, 2, length(n_method), length(n_subnat), n_years))
  for(p in 1:length(n_subnat)){ # province loop matched to C
    for(m in 1:length(n_method)){ # method loop
      alpha_pms_samp <- alpha_pms[,grep(paste0('a0\\[',m,',',p,'\\]'), colnames(alpha_pms))]
      beta.k_samp <- beta.k[,grep(paste0('a\\[',m,',',p,','), colnames(beta.k))]
      for(t in 1:n_years){
        x=n_samps*c
        i=0
        for(s in c(x-1999):x){
          i=i+1
          z[i,m,p,t] <- alpha_pms_samp[s] + B.ik[p,t,]%*%beta.k_samp[s,1:13]
        }
        P[,1,m,p,t] <- 1/(1+exp(-(z[,m,p,t])))
        P[,2,m,p,t] <- 1-P[,1,m,p,t]
      } # end t loop
    } # end M loop
  } # end P loop
  saveRDS(z, paste0('results/JAGS/JAGS_mod_MVN_Bspline_z_samps_chain',c,'.RDS'))
  saveRDS(P, paste0('results/JAGS/JAGS_mod_MVN_Bspline_P_samps_chain',c,'.RDS'))
} # end C loop


P1 <- readRDS('results/JAGS/JAGS_mod_MVN_Bspline_P_samps_chain1.RDS')
P2 <- readRDS('results/JAGS/JAGS_mod_MVN_Bspline_P_samps_chain2.RDS')
P3 <- readRDS('results/JAGS/JAGS_mod_MVN_Bspline_P_samps_chain3.RDS')

P_array <- P1 %>% abind::abind(P2, P3, along = 1)
saveRDS(P_array, 'results/JAGS/P_array_MVN_bspline.RDS')
