library(rstan)
library(shinystan)
library(tidyverse)
library(tidybayes)
library(bayesplot)

source('code/stan_utility.R')

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

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())
load("data/simulated_data/simulated_data_MVN_kstar_Kenya_new.RData")

# source('code/load_functions.R')
# source('code/2sector_code/read_in_subnational_2sector_data.R')
# source('code/2sector_code/set_up_2sector_bivar_globalrunjags.R')

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  mutate(Public.SE = 0.1) %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- as.vector(as.numeric(P_sim_df_sample$index_country))
n_all_years <- length(all_years)
M_count = 5

Bik_array <- array(NA, dim=c(P,n_all_years,K))
for(i in 1:P) {
  Bik_array[i,,] <- Bik
  
}

# Fit model -----------------------------------------
fit <- readRDS('results/STAN_model_test_kstar_simdata_MVN_alpha_N_delta_KenyaSim_new.RDS')

code <- get_stancode(fit)
cat(code)

check_n_eff(fit)

check_rhat(fit)

check_divergences(fit)

check_treedepth(fit)

check_energy(fit)

check_div(fit)

shinystan::launch_shinystan(fit)

# Traceplots
traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigmabeta_tau"), inc_warmup = FALSE, nrow = 6)
# traceplot(fit, pars = c("sigmabeta_Omega"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("beta_c"), inc_warmup = FALSE, nrow = 5)

# Check divergences 
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$`alpha_pms[1,1]`, nondiv_params$`beta_c[1,1]`,
     col=c_dark, pch=16, cex=0.8)
points(div_params$`alpha_pms[1,1]`, div_params$`beta_c[1,1]`,
       col=green, pch=16, cex=0.8)


# Rhats 
rhats <- rhat(fit)
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
beta.k_samps <- model_samps$beta_k
dim(beta.k_samps)
sum(beta.k_samps[1,1,1,1:13])

# Review spline coefficients
beta_k.mean <- apply(beta.k_samps, c(2,3,4), mean)

# testing splines ------------------------------------------------------------
all_years <- -10:30
B <- bs_bbase_precise(all_years)
Bik <- B$B.ik
K <-dim(Bik)[2]
H = K-1
kstar = B$Kstar

year_index_table <- tibble(Year = all_years, index_year = 1:length(all_years))

P_sim_df_sample <- P_sim_df_sample %>% 
  rename(Year = index_year) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(year_index_table)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
simmatchcountry <- as.vector(as.numeric(P_sim_df_sample$index_country))
n_all_years <- length(all_years)
M_count = 5

# Check out variance-covariance matrices 
View(get_posterior_mean(fit))

# Compare posterior to inputs 
get_variables(fit)

# fit %>%
#   spread_draws(`sigma_alpha[1]`, `sigma_alpha[2]`, `sigma_alpha[3]`, `sigma_alpha[4]`, `sigma_alpha[5]` ) %>%
#   pivot_longer(cols = c(`sigma_alpha[1]`, `sigma_alpha[2]`, `sigma_alpha[3]`, `sigma_alpha[4]`, `sigma_alpha[5]`), names_to='Parameter', values_to = 'sample') %>%
#   ggplot() +
#   stat_halfeye(aes(y = Parameter, x = sample)) +
#   geom_vline(xintercept = 1, col='red')
# ggsave(filename = 'visualisations/simulated_data/kstar/all_N/sigma_alpha_density.pdf')

beta_c_new <- t(as.matrix(beta_c_new))
colnames(beta_c_new) <- 1:5
beta_c_df <- as_tibble(beta_c_new) %>%
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
ggsave(filename = 'visualisations/simulated_data/kstar/MVN/beta_c_density.pdf')

#HERE
alpha_sim <- t(alpha_sim)
colnames(alpha_sim) <- 1:5
alpha_df <- as_tibble(alpha_sim) %>%
  mutate(index_subnat = 1:P) %>%
  pivot_longer(cols=c(`1`:`5`), names_to = 'index_method', values_to ='value') %>%
  rowwise() %>%
  mutate(Parameter = paste0("alpha_pms[",index_subnat,",",index_method,"]"))

fit %>%
  spread_draws(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
               `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
               `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
               `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
               `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`
  ) %>%
  pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[2,1]`, `alpha_pms[3,1]`, `alpha_pms[4,1]`, `alpha_pms[5,1]`, `alpha_pms[6,1]`,
                        `alpha_pms[1,2]`, `alpha_pms[2,2]`, `alpha_pms[3,2]`, `alpha_pms[4,2]`, `alpha_pms[5,2]`, `alpha_pms[6,2]`,
                        `alpha_pms[1,3]`, `alpha_pms[2,3]`, `alpha_pms[3,3]`, `alpha_pms[4,3]`, `alpha_pms[5,3]`, `alpha_pms[6,3]`,
                        `alpha_pms[1,4]`, `alpha_pms[2,4]`, `alpha_pms[3,4]`, `alpha_pms[4,4]`, `alpha_pms[5,4]`, `alpha_pms[6,4]`,
                        `alpha_pms[1,5]`, `alpha_pms[2,5]`, `alpha_pms[3,5]`, `alpha_pms[4,5]`, `alpha_pms[5,5]`, `alpha_pms[6,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
  facet_wrap(~Parameter)
ggsave(filename = 'visualisations/simulated_data/kstar/MVN/alpha_pm_density_KenyaSim.pdf')

# fit %>%
#   spread_draws(`alpha_pms[1,1]`, `alpha_pms[1,2]`, `alpha_pms[1,3]`, `alpha_pms[1,4]`, `alpha_pms[1,5]`, `alpha_pms[1,6]`,
#                `alpha_pms[2,1]`, `alpha_pms[2,2]`, `alpha_pms[2,3]`, `alpha_pms[2,4]`, `alpha_pms[2,5]`, `alpha_pms[2,6]`,
#                `alpha_pms[3,1]`, `alpha_pms[3,2]`, `alpha_pms[3,3]`, `alpha_pms[3,4]`, `alpha_pms[3,5]`, `alpha_pms[3,6]`,
#                `alpha_pms[4,1]`, `alpha_pms[4,2]`, `alpha_pms[4,3]`, `alpha_pms[4,4]`, `alpha_pms[4,5]`, `alpha_pms[4,6]`,
#                `alpha_pms[5,1]`, `alpha_pms[5,2]`, `alpha_pms[5,3]`, `alpha_pms[5,4]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`
#                ) %>%
#   pivot_longer(cols = c(`alpha_pms[1,1]`, `alpha_pms[1,2]`, `alpha_pms[1,3]`, `alpha_pms[1,4]`, `alpha_pms[1,5]`, `alpha_pms[1,6]`,
#                         `alpha_pms[2,1]`, `alpha_pms[2,2]`, `alpha_pms[2,3]`, `alpha_pms[2,4]`, `alpha_pms[2,5]`, `alpha_pms[2,6]`,
#                         `alpha_pms[3,1]`, `alpha_pms[3,2]`, `alpha_pms[3,3]`, `alpha_pms[3,4]`, `alpha_pms[3,5]`, `alpha_pms[3,6]`,
#                         `alpha_pms[4,1]`, `alpha_pms[4,2]`, `alpha_pms[4,3]`, `alpha_pms[4,4]`, `alpha_pms[4,5]`, `alpha_pms[4,6]`,
#                         `alpha_pms[5,1]`, `alpha_pms[5,2]`, `alpha_pms[5,3]`, `alpha_pms[5,4]`, `alpha_pms[5,5]`, `alpha_pms[5,6]`), names_to='Parameter', values_to = 'sample') %>%
#   ggplot() +
#   geom_density(aes(x = sample), fill='grey', alpha=0.8) +
#   geom_vline(data=alpha_df, aes(colour = Parameter, xintercept = value), show.legend = FALSE) +
#   facet_wrap(~Parameter)
# ggsave(filename = 'visualisations/simulated_data/kstar/all_N/alpha_pm_density_KenyaSim.pdf')


# Get P estimates 
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
colnames(P_samps.mean) <- c('index_method', 'index_subnat', 'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps <- model_samps$P
dim(P_samps)
P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
colnames(P_samps.quantile) <- c('index_method', 'index_subnat','index_sector', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) %>%
  # left_join(index_country_subnat_tbl) %>%
  left_join(year_index_table)

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)

# Get alpha intercepts 
alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_samps.mean <- t(as_tibble(apply(alpha_samps, c(2,3), mean)))
colnames(alpha_samps.mean) <- c(1:P)   # c(1:M)

alpha_samps.mean <- as_tibble(alpha_samps.mean) %>% 
  mutate(index_method = 1:5) %>%
  pivot_longer(cols=`1`:`6`, names_to = 'index_subnat', values_to = 'alpha_mean') %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  # left_join(index_country_subnat_tbl) %>%
  mutate(invlogit.alpha = exp(alpha_mean)/(1+exp(alpha_mean)))


# Get observed data 
P_df<- P_sim_df_sample %>%
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

# for(i in 1:P) {
  ggplot() +
    geom_point(data = P_df , aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
    geom_line(data = P_df, aes(x=index_year, y=Observed, colour=Sector, lty=Sector)) +
    geom_line(data = P_samps_df, aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
    geom_ribbon(data = P_samps_df, aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    facet_wrap(~interaction(Method, index_subnat), ncol = 5)
  # ggsave(filename = paste0('visualisations/simulated_data/kstar/MVN/prop_plot_',i,'_KenyaSim.pdf'))

# }

# P_df<- FP_source_data_wide %>% 
#   select(Country, Region, Method, average_year, Public, Private) %>%
#   pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')
# 
# P_df_SE <- FP_source_data_wide %>% 
#   select(Country, Region, Method, average_year, Public.SE, Private.SE) %>%
#   pivot_longer(cols = c(Public.SE, Private.SE), names_to = 'Sector', values_to = 'SE') %>%
#   mutate(Sector = str_replace(Sector, '.SE', ''))
# 
# P_df <- P_df %>% 
#   left_join(P_df_SE) %>%
#   mutate(lower_95 = Observed - 2*SE,
#          upper_95 = Observed + 2*SE)
# 
# for(i in n_subnat) {
#   c = P_df %>% filter(Region==i) %>% select(Country) %>% distinct() %>% unlist() %>% as.vector()
#   # Plot means vs observed values
#   ggplot() +
#     geom_point(data = P_df %>% filter(Region==i), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
#     geom_errorbar(data = P_df %>% filter(Region==i), aes(x=average_year, ymin = lower_95, ymax = upper_95, colour=Sector, pch=Sector)) +
#     geom_line(data = P_samps_df %>% filter(Region==i), aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
#     geom_ribbon(data = P_samps_df %>% filter(Region==i), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
#     facet_wrap(~Method)
#   ggsave(filename = paste0('visualisations/country_data/',i,'_',c, '.pdf'))
# 
# }

