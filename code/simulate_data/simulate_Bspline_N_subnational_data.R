# set seed
set.seed(1209)
library(tidyverse)
library(splines)

M = 5
P = 4
C = 1
n_years = 20
n_method <- c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills" ) # As per the method correlation matrix
method_index_table <- tibble(Method = n_method, index_method = 1:length(n_method))
index_sector_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)

beta_c_new <- c(0.1, 0.2, 0.15, 0, -0.2)

# simulate the province-level parameters 
matchcountry = rep(1,P)
sigma_alpha_sim <- rep(0.2, 6)
for(m in 1:M) {
  for(p in 1:P) {
    alpha_sim[m,p] <- rnorm(1, mean = beta_c_new[m], sd=sigma_alpha_sim[m])
  }
}

X <- seq(1, 20.5, by=0.5) # generating inputs
B <- t(splines::bs(X, degree=3, knots=c(seq(1, 20, by=4)),  intercept = TRUE)) # creating the B-splines
num_data <- length(X); num_basis <- nrow(B)
P_sim <-  array(NA, dim=c(num_data, 2, 5, P))
Y <- Y_true <- array(NA, dim=c(num_data, 5, P))
a <- array(NA, dim=c(num_basis, 5, P))

for(p in 1:P) {
  for(m in 1:M) {
    a[,m,p] <- rnorm(num_basis, 0, 1) # coefficients of B-splines
    Y_true[,m,p] <- as.vector(beta_c_new[m]*X + a[,m,p]%*%B) # generating the output
    Y[,m,p] <- Y_true[,m,p] + rnorm(length(X),0, 0.2) # adding noise
    
    # Proportions
    P_sim[,1,m,p] =  exp(Y[,m,p])/(1+exp(Y[,m,p]))
    P_sim[,2,m,p]  = 1 -  P_sim[,1,m,p] 
  }
}


P_sim_df <- plyr::adply(P_sim, c(2,3,4)) #HERE
colnames(P_sim_df) <- c('index_sector', 'index_method', 'index_subnat', 1:length(X))
P_sim_df <- P_sim_df %>%
  mutate_all(as.numeric) %>%
  pivot_longer(cols = all_of(4:(3+length(X))), names_to = 'index_year', values_to = 'Observed') %>%
  rowwise() %>%
  left_join(index_sector_table) %>%
  select(!index_sector) %>%
  pivot_wider(names_from = Sector, values_from = Observed) 


# Transform exactly 1 and 0 values away from boundary using lemon-squeezer approach ---------
P_sim_df <- P_sim_df %>%
  dplyr::mutate(Private = (Private*(nrow(P_sim_df)-1)+0.5)/nrow(P_sim_df)) %>%   # Y and SE transformation to account for (0,1) limits (total in sector)
  dplyr::mutate(Public = (Public*(nrow(P_sim_df)-1)+0.5)/nrow(P_sim_df)) %>%
  dplyr::select(index_subnat, index_method, index_year, Public, Private) %>% #, count_NA, remainder)
  mutate_all(as.numeric)

samps <- tibble(index_country = rep(1,10),
                index_year = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)) %>% 
  mutate(index_country = as.numeric(index_country), index_year = as.numeric(index_year))

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

alpha_temp <- as_tibble(alpha_sim) 
colnames(alpha_temp) <- 1:P

alpha_temp <- alpha_temp %>%
  mutate(index_method = 1:M) %>%
  pivot_longer(cols=`1`:`6`, names_to = 'index_subnat', values_to = 'alpha') %>%
  mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha)))

ggplot() +
  geom_point(data = P_df , aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_df, aes(x=index_year, y=Observed, colour=Sector, lty=Sector)) +
  facet_wrap(~interaction(Method, index_subnat), ncol = 5) 


