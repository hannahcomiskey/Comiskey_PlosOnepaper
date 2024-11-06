library(rstan)
library(shinystan)
library(tidyverse)
library(tidybayes)
library(bayesplot)

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

load("data/simulated_data_N_alpha.RData")

# Get logit of parameters and variance --------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  mutate(Public.SE = 0.1) %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# # testing splines ---------------------------------
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B
D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh 

# Set up model inputs -----------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)
D=3
M_count = 5
OD_count  = D*(M_count-D)+ D*(D-1)/2


# Fit model -----------------------------------------
fit <- readRDS('results/STAN_model_test_zih_simdata_N_alpha_N_delta.RDS')

# shiny::launch_shinystan(fit)

traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigmabeta_tau"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigmabeta_Omega"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("beta_c"), inc_warmup = FALSE, nrow = 5)

pairs(fit, pars = "sigma_y")

code <- get_stancode(fit)
cat(code)

check_divergences(fit)

# Rhats 
rhats <- rhat(fit)
print(rhats[grep('delta', names(rhats))])
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)


# N-eff
ratios <- bayesplot::neff_ratio(fit)
print(ratios)
mcmc_neff(ratios)

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)

fit_summary <- summary(fit)
View(print(fit_summary$summary))

# Get parameter estimates 

model_samps <- rstan::extract(fit)

# Check delta.k sums to 0 
delta.k_samps <- model_samps$delta_k
dim(delta.k_samps)
sum(delta.k_samps[1,1,1,1:9])

# Review spline coefficients
delta_k.mean <- apply(delta.k_samps, c(2,3,4), mean)

# testing splines ------------------------------------------------------------
all_years <- -10:30
B <- splines::bs(all_years, df=10, degree=3, intercept = FALSE)
K <-dim(B)[2]
B.ik <- B
D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Zih <- B.ik%*%Q.kh 
knots.k <- c(-5, 0, 5, 10, 15, 20, 25)
year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

P_sim_df_sample <- P_sim_df_sample %>% 
  rename(Year = index_year) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(year_index_table)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- matchcountry
n_all_years <- length(all_years)

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

random_spline_comp <- Zih %*% delta_k.mean[2,3,]
## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,random_spline_comp[,1], type= "n", xaxt="n",
     xlab = "Year", ylab ="Basis Function",
     xlim = range(all_years))
axis(1, at = min(all_years):max(all_years))
abline(v=knots.k, col = seq(1, H), lwd = 1)
lines(all_years, random_spline_comp, type= "l", col = k, lwd = 1)

# Check out variance-covariance matrices 
View(get_posterior_mean(fit))

# Compare posterior to inputs 

get_variables(fit)

fit %>%
  spread_draws(`sigma_y[1]`, `sigma_y[2]`, `sigma_y[3]`, `sigma_y[4]`, `sigma_y[5]`) %>%
  pivot_longer(cols = c(`sigma_y[1]`, `sigma_y[2]`, `sigma_y[3]`, `sigma_y[4]`, `sigma_y[5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  stat_halfeye(aes(y = Parameter, x = sample)) 
ggsave(filename = 'visualisations/simulated_data/independent_deltas/sigma_y_density.pdf')

beta_c_sim <- t(beta_c_sim)
colnames(beta_c_sim) <- 1:5
beta_c_df <- as_tibble(beta_c_sim) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_country = rep(1, each=1),
         Parameter = paste0("beta_c[",rep(1, each=5),",",rep(1:5, 1),"]"))

fit %>%
  spread_draws(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`) %>%
               # `beta_c[2,1]`, `beta_c[2,2]`, `beta_c[2,3]`, `beta_c[2,4]`, `beta_c[2,5]`,
               # `beta_c[3,1]`, `beta_c[3,2]`, `beta_c[3,3]`, `beta_c[3,4]`, `beta_c[3,5]`,
               # `beta_c[4,1]`, `beta_c[4,2]`, `beta_c[4,3]`, `beta_c[4,4]`, `beta_c[4,5]`
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  #stat_halfeye(aes(y = Parameter, x = sample)) +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=beta_c_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = 'visualisations/simulated_data/independent_deltas/beta_c_density_fixedQ.pdf')

#HERE
alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  mutate(index_subnat = rep(1:P, each=5),
         Parameter = paste0("alpha_pms[",rep(1:P, each=5),",",rep(1:5, 5),"]"))

fit %>%
  spread_draws(`alpha_pms[1,1]`, `alpha_pms[1,2]`, `alpha_pms[1,3]`, `alpha_pms[1,4]`, `alpha_pms[1,5]`,
               `alpha_pms[2,1]`, `alpha_pms[2,2]`, `alpha_pms[2,3]`, `alpha_pms[2,4]`, `alpha_pms[2,5]`,
               `alpha_pms[3,1]`, `alpha_pms[3,2]`, `alpha_pms[3,3]`, `alpha_pms[3,4]`, `alpha_pms[3,5]`,
               `alpha_pms[4,1]`, `alpha_pms[4,2]`, `alpha_pms[4,3]`, `alpha_pms[4,4]`, `alpha_pms[4,5]`,
               `alpha_pms[5,1]`, `alpha_pms[5,2]`, `alpha_pms[5,3]`, `alpha_pms[5,4]`, `alpha_pms[5,5]`) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[1,2]`, `alpha_pms[1,3]`, `alpha_pms[1,4]`, `alpha_pms[1,5]`,
                        `alpha_pms[2,1]`, `alpha_pms[2,2]`, `alpha_pms[2,3]`, `alpha_pms[2,4]`, `alpha_pms[2,5]`,
                        `alpha_pms[3,1]`, `alpha_pms[3,2]`, `alpha_pms[3,3]`, `alpha_pms[3,4]`, `alpha_pms[3,5]`,
                        `alpha_pms[4,1]`, `alpha_pms[4,2]`, `alpha_pms[4,3]`, `alpha_pms[4,4]`, `alpha_pms[4,5]`,
                        `alpha_pms[5,1]`, `alpha_pms[5,2]`, `alpha_pms[5,3]`, `alpha_pms[5,4]`, `alpha_pms[5,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = 'visualisations/simulated_data/independent_deltas/alpha_pm_density_fixedQ.pdf')


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

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)

# Get alpha intercepts 
alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_samps.mean <- as_tibble(apply(alpha_samps, c(2,3), mean)) 
colnames(alpha_samps.mean) <- c(1:20)
alpha_samps.mean <- alpha_samps.mean %>%
  pivot_longer(cols=everything(), names_to = 'index_subnat', values_to = 'alpha_mean')
alpha_samps.mean <- alpha_samps.mean %>%
  mutate(index_method = rep(1:5, each=P))
alpha_samps.mean <- alpha_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  mutate(invlogit.alpha = exp(alpha_mean)/(1+exp(alpha_mean)))

apply(model_samps$sigma_alpha, 2, mean)

# Get observed data 
P_df<- P_sim_df_sample %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

subnatid = 4
# Plot means vs observed values
ggplot() +
  geom_point(data = P_df %>% filter(index_subnat==subnatid), aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  geom_hline(data = alpha_samps.mean %>% filter(index_subnat==subnatid), aes(yintercept = invlogit.alpha)) +
  facet_wrap(~Method)



