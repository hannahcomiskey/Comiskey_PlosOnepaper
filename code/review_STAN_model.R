library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/read_in_subnational_SEdata.R")
source("code/set_up_globalrunjags.R")

# # testing splines ---------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B
D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh 
H <- dim(Zih)[2]

# Fit model -----------------------------------------
fit <- readRDS(file='results/global_subnat/STAN_model_zih_deltak_factorcov_obsdata.RDS')

traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]", "alpha_pms[3,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("L_Sigma[1,1]", "L_Sigma[2,2]", "L_Sigma[3,3]"), inc_warmup = FALSE, nrow = 3)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("Q"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)


code <- get_stancode(fit)
cat(code)

check_divergences(fit)

# Rhats 
rhats <-  bayesplot::rhat(fit)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)


# N-eff
ratios <- bayesplot::neff_ratio(fit)
print(ratios)
mcmc_neff(ratios)

plot(fit, pars=c('mu_lt'))

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)


# Get parameter estimates ------------------
method_index_table <- tibble(Method = n_method, index_method = 1:5)
sector_index_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

# get alpha params
alpha_samps <- rstan::extract(fit, pars = "alpha_pms")
alpha_samps <- alpha_samps$alpha_pms
alpha_mean <- t(apply(alpha_samps, c(2,3), mean))
colnames(alpha_mean) <- c(1:5)
alpha_mean <- as_tibble(alpha_mean)
alpha_df <- alpha_mean %>% 
  mutate(index_subnat = 1:nrow(alpha_mean)) %>%
  pivot_longer(cols = c(`1`:`5`), names_to = 'index_method', values_to = 'alpha_mean')

delta_samps <- rstan::extract(fit, pars = "delta_k")
delta_samps <- delta_samps$delta_k
delta_k.mean <- apply(delta_samps, c(2,3,4), mean)
knots.k <- c(1994.188,1998.375,2002.562,2006.750,2010.938,2015.125,2019.312)

## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,Zih[,1], type= "n", xaxt="n",
     xlab = "Year",
     ylab ="Basis Function",
     xlim = range(all_years))
axis(1, at = min(all_years):max(all_years))
abline(v=knots.k, col = seq(1, K), lwd = 1)
for (k in 1:H){
  lines(all_years,Zih[,k], type= "l", col = k, lwd = 1)
}

random_spline_comp <- Zih %*% delta_k.mean[5,109,]
## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,random_spline_comp[,1], type= "n", xaxt="n",
     xlab = "Year", ylab ="Basis Function",
     xlim = range(all_years))
axis(1, at = min(all_years):max(all_years))
abline(v=knots.k, col = seq(1, H), lwd = 1)
lines(all_years, random_spline_comp, type= "l", col = k, lwd = 1)

delta_df <- plyr::adply(delta_k.mean, c(1,2))
colnames(delta_df) <- c('index_method', 'index_subnat', 1:9)

# Check delta.k sums to 0 
delta_df %>%
  filter(index_method==1 & index_subnat==1) %>% 
  select(`1`:`9`) %>% 
  unlist() %>% as.vector() %>% sum()

# Assemble estimated proportions from parameters 
M_count=5
P_count = nrow(country_subnat_tbl)
S_count=2
z <- array(NA, dim=c(M_count, P_count, n_all_years))
P <- array(NA, dim=c(M_count, P_count, S_count, n_all_years))
intercept = rep(1, n_all_years)
inv_logit <- function(x) {
  inv.logit <- exp(x)/(1+exp(x))
  return(inv.logit)
}
for(m in 1:M_count){ 
  for(p in 1:P_count){
    alpha_pms <- as.vector(unlist(alpha_df %>% filter(index_method==m & index_subnat==p) %>% select(alpha_mean)))
    delta_k <- as.vector(unlist(delta_df  %>% filter(index_method==m & index_subnat==p) %>% select(`1`:`9`)))
    for(t in 1:n_all_years) {
      z[m,p,t] <- alpha_pms + sum(Zih[t,]*delta_k)
    }
    P[, ,1,] = inv_logit(z)
    P[, ,2,] =1 - P[, ,1,]
  } 
} 

P_mean <- plyr::adply(P, c(1,2,3,4))
colnames(P_mean) <- c('index_method', 'index_subnat', 'index_sector', 'index_year', 'Mean')
P_mean <- P_mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(index_country_subnat_tbl) %>% 
  left_join(method_index_table) %>% 
  left_join(sector_index_table) %>%
  left_join(year_index_table)

FP_long <- FP_source_data_wide %>%
  rowwise() %>%
  mutate(Private = Commercial_medical + Other) %>%
  #select(!c(Public_n, Private_n, check_total)) %>%
  select(!c(count_SE.NA, DEFT, Commercial_medical.SE, Other.SE, Public.SE, Commercial_medical, Other, n_Public, n_Commercial_medical, n_Other, check_total)) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')


# # Plot means vs observed values
ggplot() +
  geom_point(data = FP_long %>% filter(index_subnat==200), aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_mean %>% filter(index_subnat==200), aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  facet_wrap(~Method)



