# set seed
set.seed(1209)
library(tidyverse)

M = 5
P = 6
C = 1
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix
method_index_table <- tibble(Method = n_method, index_method = 1:length(n_method))
index_sector_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)

# Set up the parameters to construct the proportions
# beta_c_sim <- matrix(NA, nrow = M, ncol = C)
# 
# sd_beta <- c(1, 0.5, 0.3 , 1.5, 1.1)
# 
# # Simulate the country-level parameters 
# for(c in 1:ncol(beta_c_sim)) {
#   for(m in 1:M) {
#     beta_c_sim[m,c] <-  rnorm(1, mean = 0, sd=sd_beta[m])
#   }
# }

alpha_sim <- matrix(NA, nrow = M, ncol = P)
beta_c_mean <- readRDS('data/simulated_data/mle_mean_beta.RDS') %>%
  filter(Country=='Kenya') %>%
  ungroup() %>%
  select(!Country)

beta_c_new <- c(0.9, 0.8, 0.5, 0.6, -0.2)

# simulate the province-level parameters 
matchcountry = rep(1,P)
sigma_alpha_sim <- rep(0.2, 6)
for(m in 1:M) {
  for(p in 1:P) {
    alpha_sim[m,p] <- rnorm(1, mean = beta_c_new[m], sd=sigma_alpha_sim[m]) # as.vector(unlist(beta_c_mean[, m])), sd=sigma_alpha_sim[m])
  }
}

all_years = 1:10 #1:50 #1:20
n_years = length(all_years)
# # testing splines ---------------------------------

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


B <- bs_bbase_precise(all_years)
Bik <- B$B.ik
K <-dim(Bik)[2]
H = K-1
sigma_y = 0.5
kstar = B$Kstar

P_sim <- array(NA, dim=c(M,P,n_years,2))
Z_tmp <- Z <- array(NA, dim=c(M,P,n_years))
delta_sim <- array(NA, dim=c(M, P, H))
beta_k_sim <- array(NA, dim=c(M, P, K))

# Simulate FOD spline coefficients 
var_delta <- readRDS('data/simulated_data/mle_var_delta.RDS')

for(m in 1:M){
  for(p in 1:P) {
    for(h in 1:H) {
      delta_sim[m,p,h] <- rnorm(1, mean=0, sd= 0.4) #sqrt(var_delta[m]))
    }
  }
}


for(m in 1:M){
  for(p in 1:P){
    # Spline coefficients
    beta_k_sim[m,p,kstar] = 0 
    for(j in (kstar+1):K) {
      beta_k_sim[m,p,j] = beta_k_sim[m,p,j-1] + delta_sim[m, p, j-1] } 
      for(j in 1:(kstar-1)) {
        t = kstar - j
        beta_k_sim[m,p,t] = beta_k_sim[m,p,t+1] - delta_sim[m,p, t]
      } # before kstar
    
    # Latent variable
    for(t in 1:n_years) {
      Z_tmp[m,p,t] = as.numeric(alpha_sim[m,p] + Bik[t,1:K]%*%beta_k_sim[m,p,1:K]) # Public sector proportion on logit scale
      Z[m,p,t] = rnorm(1, Z_tmp[m,p,t], 0.1)
      # Proportions
      P_sim[m, p, t, 1] =  exp(Z[m, p,t])/(1+exp(Z[m, p,t]))
      P_sim[m, p, t, 2] = 1 - P_sim[m, p, t, 1]
    }
  }
}


P_sim_df <- plyr::adply(P_sim, c(2,3,4))
index_subnat_table <- tibble(index_country = as.numeric(matchcountry), index_subnat = as.numeric(1:P))
colnames(P_sim_df) <- c('index_subnat','index_year', 'index_sector', 1:M)
P_sim_df <- P_sim_df %>%
  mutate_all(as.numeric) %>%
  pivot_longer(cols = all_of(4:(3+M)), names_to = 'index_method', values_to = 'Observed') %>%
  rowwise() %>%
  left_join(index_subnat_table) %>%
  left_join(index_sector_table) %>%
  select(!index_sector) %>%
  pivot_wider(names_from = Sector, values_from = Observed) 

P_sim_df <- P_sim_df %>%
  mutate(index_country = rep(C, nrow(P_sim_df)))

# Transform exactly 1 and 0 values away from boundary using lemon-squeezer approach ---------
P_sim_df <- P_sim_df %>%
  dplyr::mutate(Private = (Private*(nrow(P_sim_df)-1)+0.5)/nrow(P_sim_df)) %>%   # Y and SE transformation to account for (0,1) limits (total in sector)
  dplyr::mutate(Public = (Public*(nrow(P_sim_df)-1)+0.5)/nrow(P_sim_df)) %>%
  dplyr::select(index_country, index_subnat, index_method, index_year, Public, Private) #, count_NA, remainder)

# samps <- tibble(index_country = rep(1,25),
#                 index_year = seq(1,50, by=2)) %>% 
#   mutate(index_country = as.numeric(index_country), index_year = as.numeric(index_year))

samps <- tibble(index_country = rep(1,5),
                index_year = seq(1,10, by=2)) %>% 
  mutate(index_country = as.numeric(index_country), index_year = as.numeric(index_year))


P_sim_df_sample <- left_join(samps, P_sim_df)

# Get observed data 
P_df<- P_sim_df_sample %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

ggplot() +
  geom_point(data = P_df , aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_df, aes(x=index_year, y=Observed, colour=Sector, lty=Sector)) +
  #geom_hline(data = alpha_temp, aes(yintercept = invlogit.alpha)) +
  facet_wrap(~interaction(Method, index_subnat), ncol = 5)


