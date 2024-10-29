# set seed
set.seed(1209)
library(tidyverse)

M = 5
P = 5
C = 1
H = 9
t = 20
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix
method_index_table <- tibble(Method = n_method, index_method = 1:length(n_method))

# Set up the parameters to construct the proportions
alpha_sim <- matrix(NA, nrow = M, ncol = P)
beta_c_sim <- matrix(NA, nrow = M, ncol = C)
delta_sim <- array(NA, dim=c(M, P, H))

sd_beta <- c(1, 0.5, 0.3 , 1.5, 1.1)

# Simulate the country-level parameters 
for(c in 1:ncol(beta_c_sim)) {
  for(m in 1:M) {
    beta_c_sim[m,c] <-  rnorm(1, mean = 0, sd=sd_beta[m])
  }
}

# simulate the province-level parameters 
#matchcountry <- c(rep(1, 7), rep(2, 4), rep(3, 4), rep(4,5))
matchcountry = rep(1,5)
sigma_alpha_sim <- c(1, 1, 1.2, 1, 0.6)
for(m in 1:M) {
  for(p in 1:P) {
    alpha_sim[m,p] <- rnorm(1, mean = beta_c_sim[m,matchcountry[p]], sd=sigma_alpha_sim[m])
  }
}


# Simulate FOD spline coefficients 
sd_delta <- c(0.6, 0.8, 0.5, 0.3, 0.75)

for(m in 1:M){
  for(p in 1:P) {
    for(h in 1:(H-1)) {
      delta_sim[m,p,h] <- rnorm(1, mean=0, sd=sd_delta[m])
    }
    delta_sim[m,p,H] <- -sum(delta_sim[m,p,c(1:(H-1))])
  }
}

# Z parameters 
P_sim <- Z <- Z_tmp <- array(NA, dim=c(P, t, M))

all_years = 1:20

# # testing splines ---------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B
D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh 
sigma_y = 0.5

for(m in 1:M) {
  for(p in 1:P) {
    for(t in 1:t) { 
      Z_tmp[p,t,m] <- alpha_sim[m,p] + sum(Zih[t,]*delta_sim[m,p,1:H])
      Z[p,t,m] <- rnorm(1,Z_tmp[p,t,m],sigma_y)
      P_sim[p,t,m] <- exp(Z[p,t,m])/(1+exp(Z[p,t,m]))
    }
  }
}

P_sim_df <- plyr::adply(P_sim, c(2,3))
index_subnat_table <- tibble(index_country = as.character(matchcountry), index_subnat = as.character(1:P))

colnames(P_sim_df) <- c('index_year', 'index_method', 1:P)
P_sim_df <- P_sim_df %>%
  pivot_longer(cols = all_of(3:(2+P)), names_to = 'index_subnat', values_to = 'Public') %>%
  rowwise() %>%
  mutate(Private = 1-Public) %>%
  left_join(index_subnat_table)

# Transform exactly 1 and 0 values away from boundary using lemon-squeezer approach ---------
P_sim_df <- P_sim_df %>%
  dplyr::mutate(Private = (Private*(nrow(P_sim_df)-1)+0.5)/nrow(P_sim_df)) %>%   # Y and SE transformation to account for (0,1) limits (total in sector)
  dplyr::mutate(Public = (Public*(nrow(P_sim_df)-1)+0.5)/nrow(P_sim_df)) %>%
  dplyr::select(index_country, index_subnat, index_method, index_year, Public, Private) #, count_NA, remainder)


# samps <- tibble(index_country = c(rep(1,5), 
#                                   rep(2,3), 
#                                   rep(3,2), 
#                                   rep(4,5)), 
#                 index_year = c(4,12, 15, 20, 21, 8, 16, 20, 13, 16, 2, 4, 12, 17, 20)) %>%
#   mutate(index_country = as.character(index_country), 
#          index_year = as.character(index_year))

samps <- tibble(index_country = rep(1,4),
                index_year = c(4, 8, 12, 15)) %>% #20, 21, 8, 16, 20, 13, 16, 2, 4, 12, 17, 20)) %>%
  mutate(index_country = as.character(index_country), index_year = as.character(index_year))


P_sim_df_sample <- left_join(samps, P_sim_df)

# Get observed data 
P_df<- P_sim_df_sample %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Method, pch=Sector)) +
  geom_line(data = P_df, aes(x=index_year, y=Observed, colour=Method, lty=Sector)) +
  facet_wrap(~ interaction(Method, index_subnat))

indexid=2
alpha_temp <- tibble(alpha = alpha_sim[,indexid], index_method = 1:5, Method = n_method) %>%
  mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha)))

ggplot() +
  geom_point(data = P_df %>% filter(index_subnat==indexid), aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_df %>% filter(index_subnat==indexid), aes(x=index_year, y=Observed, colour=Sector, lty=Sector)) +
  #geom_hline(data = alpha_temp, aes(yintercept = invlogit.alpha)) +
  facet_wrap(~Method)


