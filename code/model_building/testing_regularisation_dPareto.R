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
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  vector[M_count] beta_c; // expected mean trend
  }

parameters {   // The parameters accepted by the model.
  vector<lower=0>[M_count] sigma_delta; // variance of mean trend
  vector<lower=0>[M_count] sigma_alpha; // variance of mean trend
  vector<lower=0>[M_count] sigma_y;
  vector[M_count] alpha_raw[P_count]; // non-centered parameter for hierarchy
  vector<lower=0>[H] lambda; // local shrinkage parameter
  vector<lower=0>[H] nu; // local shrinkage parameter
  vector[H] delta_k[P_count, M_count]; // variation associated with time
}

transformed parameters { 
  vector[M_count] alpha_pms[P_count]; // expected mean trend
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] P[M_count, P_count]; // logit observation

  for(m in 1:M_count){ 
    for(p in 1:P_count){
      alpha_pms[p,m] = beta_c[m] + sigma_alpha[m]*alpha_raw[p,m];
      // Spline coefficients
      beta_k[m,p,kstar[p]] = zero; // set spline coefficient to 0
      for(j in (kstar[p]+1):K) {
        beta_k[m,p,j] = beta_k[m,p,j-1] + delta_k[p,m, j-1];
      } // after kstar
      for(j in 1:(kstar[p]-1)) { // Estimating spline coefficient here
        int t = kstar[p] - j;
        beta_k[m,p,t] = beta_k[m,p,t+1] - delta_k[p,m, t];
      } // before kstar
      
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
  sigma_y ~ normal(0, 2);
  
  // Hierarchical estimation of intercept
  for(p in 1:P_count){
    alpha_raw[p] ~ normal(0, 1); // sharing info across methods within a province so each province public/private sector has an intercept.
  } // end P loop
  
  for(h in 1:H) {
    lambda[h] ~ gamma(1, 1);
    nu[h] ~ exponential(square(lambda[h])/2);
  }
  for(m in 1:M_count) {
    for(p in 1:P_count){  
      sum(beta_k[m,p]) ~ normal(alpha_pms[p,m], inv_sqrt(1 - inv(K)));  // Sum-to-0 constraint on spline coefficients 
      for(h in 1:H) {
        delta_k[p,m,h] ~ normal(0, sigma_delta[m]*nu[h]);
      }
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
all_years <- seq(min(small_samp$index_year)-5, max(small_samp$index_year)+5, by=1)
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
P = unique(small_samp$index_subnat)
simmatchsubnat <- as.vector(as.numeric(small_samp$index_subnat))
simmatchmethod <- as.vector(as.numeric(small_samp$index_method))
simmatchyears <- as.vector(as.numeric(small_samp$index_year))
simmatchcountry <- as.vector(1)
n_all_years <- length(all_years)
M_count = 5

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
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears,
                  scale_global=scale_global,
                  beta_c = beta_c_new
                  )

## Parameters to look at ------------------------------
pars <- c("alpha_pms", # required for P
          "delta_k",
          "beta_k",
          "sigma_delta",
          "sigma_alpha",
          "lambda",
          "nu",
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

fit <- readRDS( 'results/regularisation_model_testing_smallsamp.RData')

code <- get_stancode(fit)
cat(code)

check_n_eff(fit)

check_rhat(fit)

check_divergences(fit)

check_treedepth(fit)

check_energy(fit)

check_div(fit)

shinystan::launch_shinystan(fit)

model_samps <- rstan::extract(fit)

traceplot(fit, pars = c("alpha_pms"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("nu[1,1,1]", "nu[1,2,2]", "nu[1,3,4]", "nu[1,4,5]", "nu[1,5,11]" ), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("lambda[1,1,1]", "lambda[1,2,2]", "lambda[1,3,4]", "lambda[1,4,5]", "lambda[1,5,11]" ), inc_warmup = FALSE, nrow = 5)


# Check divergences
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$`alpha_pms[1,1]`, log(nondiv_params$`sigma_alpha[1]`),
     col=c_dark, pch=16, cex=0.8)
points(div_params$`alpha_pms[1,1]`, log(div_params$`sigma_alpha[1]`),
       col=green, pch=16, cex=0.8)


# Plot beta spline coefficients
betak_samps <- model_samps$beta_k
dim(betak_samps)
betak_samps.mean <- apply(betak_samps, c(2,3,4), mean)
dim(betak_samps.mean)

betak_test <- tibble(Beta_k = betak_samps.mean[5,1,], H=1:13)
knots.k <- B$knots.k
ggplot() +
  geom_line(data = betak_test, aes(x=H, y=Beta_k))

random_spline_comp <- Bik %*% betak_samps.mean[5,1,]

## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,random_spline_comp[,1], type= "n", xaxt="n",
     xlab = "Year", ylab ="Basis Function",
     xlim = range(all_years))
axis(1, at = min(all_years):max(all_years))
abline(v=knots.k, col = seq(1, H), lwd = 1)
lines(all_years, random_spline_comp, type= "l", lwd = 1)

alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  mutate(index_subnat = 1:6) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  rowwise() %>%
  mutate(Parameter = paste0("alpha_pms[",index_subnat,",",index_method,"]")) %>%
  filter(index_subnat==1)

fit %>%
  spread_draws(`alpha_pms[1,1]`, #`alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
               `alpha_pms[1,2]`, #`alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
               `alpha_pms[1,3]`, #`alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
               `alpha_pms[1,4]`, #`alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
               `alpha_pms[1,5]` # `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`
  ) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, #`alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
                        `alpha_pms[1,2]`, #`alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
                        `alpha_pms[1,3]`, #`alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
                        `alpha_pms[1,4]`, #`alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
                        `alpha_pms[1,5]` # `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`
                        ), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = 'visualisations/simulated_data/kstar/all_N/ncp_reg_GDP_alpha_pm_density_KenyaSim.pdf')

# Get P estimates
P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private')
P_samps.mean <- P_samps.mean %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

sector_index_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)
P_samps <- model_samps$P
dim(P_samps)
P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
colnames(P_samps.quantile) <- c('index_method', 'index_subnat','index_sector', 'index_year', 'lower_95', 'upper_95')
P_samps.quantile <- P_samps.quantile %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table)

P_samps_df <- left_join(P_samps.mean, P_samps.quantile) %>%
  left_join(year_index_table)

# Get observed data
P_df<- small_samp %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

subnatid = 1
# Plot means vs observed values
ggplot() +
  geom_point(data = P_df %>% filter(index_subnat==subnatid), aes(x=Year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=Year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=Year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~Method)
ggsave(filename = 'visualisations/simulated_data/kstar/all_N/ncp_reg_GDP_plotting_p.pdf')


