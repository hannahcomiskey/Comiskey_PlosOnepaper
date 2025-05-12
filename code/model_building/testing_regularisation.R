library(rstan)
library(shinystan)
library(tidyverse)
library(tidybayes)
library(bayesplot)
source('code/stan_utility.R')
source('code/load_functions.R')

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

load("data/simulated_data/simulated_data_all_N_kstar_Kenya_50years.RData")

model <- "
data {  
  int<lower=1> n_years; // Number of years
  int<lower=1> n_obs; // Number of observations
  int<lower=1> H; // Number of knots
  int<lower=1> K; // Number of spline cofficients (H+1)
  int<lower=1> P_count; // Number of provinces
  int<lower=1> C_count; // Number of countries
  int<lower=1> M_count; // Number of methods
  int<lower=1> S_count; // Number of sectors
  array[P_count] int<lower=1, upper=K> kstar; // Spline index K star for estimation
  vector[K] Bik[P_count, n_years]; // Basis functions
  real zero;
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  real<lower=1> nu_global ; // degrees of freedom for the half-t prior for tau
  real<lower=1> nu_local ; // degrees of freedom for the half-t priors for lambdas
  real<lower=0> scale_global ; // scale for the half-t prior for tau
  real<lower=0> slab_scale ; # slab scale for the regularized horseshoe
  real<lower=0> slab_df ; # slab degrees of freedom for the regularized horseshoe
  matrix[C_count, M_count] beta_c; // expected mean trend
  }

parameters {   // The parameters accepted by the model.
  vector[H] delta_k[P_count, M_count]; // variation associated with time
  vector<lower=0>[M_count] sigma_delta; // variance of mean trend
  vector<lower=0>[M_count] sigma_alpha; // variance of mean trend
  vector<lower=0>[M_count] sigma_y;
  vector[M_count] alpha_raw[P_count]; // non-centered parameter for hierarchy
  real<lower=0> caux;
  real logsigma;
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[H] aux1_local[P_count, M_count];
  vector<lower=0>[H] aux2_local[P_count, M_count];
}

transformed parameters { 
  vector[M_count] alpha_pms[P_count]; // expected mean trend
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] P[M_count, P_count]; // logit observation
  vector<lower=0>[H] lambda_tilde[P_count, M_count]; // 'truncated' local shrinkage parameter
  vector<lower=0>[H] lambda[P_count, M_count]; // local shrinkage parameter
  real<lower=0> sigma_tau;  // noise std
  real<lower=0> cstar; # slab scale
  real<lower=0> tau_delta; // global shrinkage parameter

  sigma_tau = exp(logsigma);
  tau_delta = aux1_global*sqrt(aux2_global)*scale_global*sigma_tau;
  cstar = slab_scale*sqrt(caux);
  
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      alpha_pms[p,m] = beta_c[matchcountry[p], m] + sigma_alpha[m]*alpha_raw[p,m];
      // Spline coefficients
      beta_k[m,p,kstar[p]] = zero; // set spline coefficient to 0
      for(j in (kstar[p]+1):K) {
        beta_k[m,p,j] = beta_k[m,p,j-1] + delta_k[p,m, j-1];
      } // after kstar
      for(j in 1:(kstar[p]-1)) { // Estimating spline coefficient here
        int t = kstar[p] - j;
        beta_k[m,p,t] = beta_k[m,p,t+1] - delta_k[p,m, t];
      } // before kstar
      lambda[p,m] = aux1_local[p,m].*sqrt(aux2_local[p,m]);
      lambda_tilde[p,m]  = sqrt((pow(cstar,2)*pow(lambda[p,m],2)) ./ (pow(cstar,2) + pow(sigma_tau,2)*pow(lambda[p,m], 2)));

      // Latent variable
      for(t in 1:n_years) {
        z[m,p,t] = alpha_pms[p,m] +  dot_product(Bik[p, t, 1:K],beta_k[m,p]); // Public sector proprtion on logit scale
      }
      
      // Proportions
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigma_delta ~ normal(0, 2);
  sigma_alpha ~ normal(0, 2);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
  logsigma ~ normal(0, 2);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  sigma_y ~ normal(0, 2);
  
  // Hierarchical estimation of intercept
  for(p in 1:P_count){
    alpha_raw[p] ~ normal(0, 1); // sharing info across methods within a province so each province public/private sector has an intercept.
  } // end P loop
  
  for(m in 1:M_count) {
    for(p in 1:P_count){  
      // Sum-to-0 constraint on spline coefficients 
      sum(beta_k[m,p]) ~ normal(alpha_pms[p,m], inv_sqrt(1 - inv(K)));
      for(h in 1:H) {
        delta_k[p,m,h] ~ normal(0, tau_delta*lambda_tilde[p,m,h]);
      }
      aux1_local[p,m] ~ normal(0, 1);
      aux2_local[p,m] ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
    }
  }


  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

"

# Get logit of parameters and variance -----------------------------------------
small_samp <- P_sim_df_sample %>% filter(index_subnat==1)

mydata <- small_samp[,c("Public")] %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)))

# # testing splines ------------------------------------------------------------
all_years <- -5:55
B <- bs_bbase_precise(all_years)
Bik <- B$B.ik
K <-dim(Bik)[2]
H = K-1
kstar = B$Kstar

year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

small_samp <- small_samp %>% 
  rename(Year = index_year) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(year_index_table)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(small_samp$index_subnat))
simmatchmethod <- as.vector(as.numeric(small_samp$index_method))
simmatchyears <- as.vector(as.numeric(small_samp$index_year))
simmatchcountry <- rep(1, P)
n_all_years <- length(all_years)
M_count = 5
P = unique(small_samp$index_subnat)

Bik_array <- array(NA, dim=c(P,n_all_years,K))
for(i in 1:P) {
  Bik_array[i,,] <- Bik
  
}

# let p0=3, D=H and n=n_obs 
# See 3.12 Vehtari paper.
scale_global = 3/((H-3)*sqrt(nrow(logit.data)))
  
## The required data ------------------------------
inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  Bik = Bik_array,
                  n_years = n_all_years,
                  n_obs = nrow(logit.data),
                  K = K,
                  H = H,
                  kstar= array(kstar, dim=c(P)),
                  zero = 0,
                  C_count = C,
                  P_count = P,
                  M_count = M,
                  S_count = 2,
                  matchsubnat = simmatchsubnat,
                  matchcountry = simmatchcountry,
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears,
                  scale_global=scale_global,
                  beta_c = matrix(beta_c_new, nrow=1, ncol=5)
                  )

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "delta_k",
          "beta_k",
          "sigma_delta",
          "sigma_beta",
          "sigma_alpha",
          "lambda_tilde",
          "tau_delta",
          "cstar",
          "sigma_tau",
          "P")

# Run stan model ------------------

fit <- stan(
  data = inputdata,    # named list of data
  model_code = model,
  pars = pars,
  iter = 10000,         # total number of iterations per chain
  warmup = 2000,
  thin=4,
  chains=3,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99)
)
saveRDS(fit, 'results/regularisation_model_testing_smallsamp.RData')

# fit <- readRDS( 'results/regularisation_model_testing.RData')
# code <- get_stancode(fit)
# cat(code)
# 
# check_n_eff(fit)
# 
# check_rhat(fit)
# 
# check_divergences(fit)
# 
# check_treedepth(fit)
# 
# check_energy(fit)
# 
# check_div(fit)
# 
# shinystan::launch_shinystan(fit)
# 
# model_samps <- rstan::extract(fit)
# 
# traceplot(fit, pars = c("alpha_pms"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)
# traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
# 
# 
# # Check divergences
# c_dark <- c("#8F272780")
# green <- c("#00FF0080")
# 
# partition <- partition_div(fit)
# div_params <- partition[[1]]
# nondiv_params <- partition[[2]]
# 
# par(mar = c(4, 4, 0.5, 0.5))
# plot(nondiv_params$`alpha_pms[2,5]`, nondiv_params$`sigma_alpha[5]`,
#      col=c_dark, pch=16, cex=0.8)
# points(div_params$`alpha_pms[2,5]`, div_params$`sigma_alpha[5]`,
#        col=green, pch=16, cex=0.8)
# 
# 
# # Plot beta spline coefficients
# betak_samps <- model_samps$beta_k
# dim(betak_samps)
# betak_samps.mean <- apply(betak_samps, c(2,3,4), mean)
# dim(betak_samps.mean)
# 
# betak_test <- tibble(Beta_k = betak_samps.mean[1,3,], H=1:13)
# 
# ggplot() +
#   geom_line(data = betak_test, aes(x=H, y=Beta_k))
# 
# random_spline_comp <- Bik %*% betak_samps.mean[2,3,]
# ## Plot basis
# par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
# plot(all_years,random_spline_comp[,1], type= "n", xaxt="n",
#      xlab = "Year", ylab ="Basis Function",
#      xlim = range(all_years))
# axis(1, at = min(all_years):max(all_years))
# abline(v=knots.k, col = seq(1, H), lwd = 1)
# lines(all_years, random_spline_comp, type= "l", col = k, lwd = 1)
# 
# subnatid = 5
# # # Get Z estimates
# #
# # # From simulation
# # z_tmp <- array(dim=c(M, P, n_years))
# # for(m in 1:M) {
# #   for(p in 1:P) {
# #     for(t in 1:n_years) {
# #       z_tmp[m,p,t] = alpha_sim[m,p] + sum(betak_samps.mean[m,p,]*Bik[t,])
# #     }
# #   }
# # }
# # z_tmp <- plyr::adply(z_tmp, .margins=c(1,2,3))
# # colnames(z_tmp) <- c('index_method', 'index_subnat', 'index_year', 'Z')
# # z_tmp <- z_tmp %>%
# #   mutate(across(everything(), as.numeric)) %>%
# #   left_join(method_index_table)
# # # Plot means vs observed values
# # ggplot() +
# #   geom_line(data = z_tmp %>% filter(index_subnat==5), aes(x=index_year, y=Z)) +
# #   facet_wrap(~Method)
# #
# # # From model
# # Z_samps <- model_samps$z
# # dim(Z_samps)
# # Z_samps.mean <- apply(Z_samps, c(2,3,4), mean)
# # Z_samps.mean <- plyr::adply(Z_samps.mean, .margins=c(1,2,3))
# # colnames(Z_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Z')
# # Z_samps.mean <- Z_samps.mean %>%
# #   mutate(across(everything(), as.numeric)) %>%
# #   left_join(method_index_table)
# # # Plot means vs observed values
# # ggplot() +
# #   geom_point(data = Z_samps.mean %>% filter(index_subnat==subnatid), aes(x=index_year, y=Z)) +
# #   geom_line(data = Z_samps.mean %>% filter(index_subnat==subnatid), aes(x=index_year, y=Z)) +
# #   facet_wrap(~Method)
# 
# # Get P estimates
# P_samps <- model_samps$P
# dim(P_samps)
# P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
# P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
# colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private')
# P_samps.mean <- P_samps.mean %>%
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')
# 
# sector_index_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)
# P_samps <- model_samps$P
# dim(P_samps)
# P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
# P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
# colnames(P_samps.quantile) <- c('index_method', 'index_subnat','index_sector', 'index_year', 'lower_95', 'upper_95')
# P_samps.quantile <- P_samps.quantile %>%
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   left_join(sector_index_table)
# 
# P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
# 
# # Get observed data
# P_df<- P_sim_df_sample %>%
#   mutate(across(everything(), as.numeric)) %>%
#   left_join(method_index_table) %>%
#   pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')
# 
# subnatid = 3
# # Plot means vs observed values
# ggplot() +
#   geom_point(data = P_df %>% filter(index_subnat==subnatid), aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
#   geom_line(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
#   geom_ribbon(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
#   facet_wrap(~Method)
# 
# 
