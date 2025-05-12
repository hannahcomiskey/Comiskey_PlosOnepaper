library(tidyverse)
library(tidybayes)
library(mcmcplots)
library(bayesplot)
library(coda)
library(shinystan)
library(rstan)

# Source code --------------------------------------
source('code/stan_utility.R')
source('code/load_functions.R')
source('code/DBDA2E-utilities.R')
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")

# Get logit of parameters and variance --------------------
mydata <- FP_source_data_wide[,c("Public", "Public.SE")]
logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)),
         logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
         logit.Public.SE = sqrt(logit.Public.Var))
M =  length(n_method)

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
year_index_table <- tibble(average_year = all_years, index_year = 1:length(all_years))

# fit model -----------------------------------------
vispath = 'visualisations/country_data/bspline/STAN/'

# fit model -----------------------------------------
fit <- readRDS('results/STAN_model_basis_expand_4country.RDS')

code <- get_stancode(fit)
cat(code)

check_all_diagnostics(fit)

check_energy(fit)
check_divergences(fit)
check_n_eff(fit)
check_rhat(fit)
check_treedepth(fit)

stan_ac(fit, pars = 'tau')
stan_ac(fit, pars = 'sigma_y')

stan_plot(fit,  pars = 'tau')

stan_hist(fit,  pars = 'tau')

stan_dens(fit,  pars = 'tau')

# Get parameter estimates 
model_samps <- rstan::extract(fit)

fit_mcmc <- coda::as.mcmc(fit)

denplot(fit_mcmc, parms = c("tau", "sigma_y"))

# Traceplots
traceplot(fit, pars = c("a[1,5,1]", "a[2,5,3]", "a[3,5,9]", "a[4,5,6]", "a[5,5,4]", "a[6,5,4]"), inc_warmup = FALSE, nrow=6)
traceplot(fit, pars = c("a0[1,1]", "a0[5,2]", "a0[15,3]", "a0[20,4]", "a0[25,5]"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("tau"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_y"), inc_warmup = FALSE, nrow = 5)

# Check divergences 
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$`a0[1,1]`, #nondiv_params$`sigma_a0[1]`,
     col=c_dark, pch=16, cex=0.8)
points(div_params$`a0[1,1]`, #div_params$`sigma_a0[1]`,
       col='blue', pch=16, cex=0.8)

# Checking for funnel in beta parameters
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("a0[1,", k, "]", sep="")
  # name2 <- paste("sigma_a0[", k, "]", sep="")
  plot(nondiv_params[name][[1]], #log(nondiv_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name) #, ylab="log(tau)")
  points(div_params[name][[1]], #log(div_params[name2][[1]]),
         col="blue", pch=16, cex=0.8)
}
par(mfrow=c(1, 1))


# Rhats 
rhats <- rhat(fit)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats) + yaxis_text(hjust = 1)


# N-eff
ratios <- bayesplot::neff_ratio(fit)
mcmc_neff(ratios)

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)

fit_summary <- summary(fit)

# Review spline coefficients
beta.k_samps <- model_samps$a
beta_k.mean <- apply(beta.k_samps, c(2,3,4), mean)
dim(beta_k.mean)
sum(beta_k.mean[1,1,])

# Posterior predictive check
y_sim <- fit$BUGSoutput$sims.list$Y_sim
y_sim_mean <- apply(y_sim, 2, mean)
y_obs <- as.vector(unlist(logit.data[,c("logit.Public")]))

pdf(paste0(vispath, "posterior_predictive_observed_4country_bspline.pdf")) 
plot(y_sim_mean, y_obs)
abline(a=0, b=1, col='red')
dev.off()

pdf(paste0(vispath,"posterior_predictive_check_4country_bspline.pdf"))
# plot the observed data hist
hist(y_obs, freq = FALSE, breaks = 50) # observed data 
n_samples = nrow(y_sim) 
for (r in 1:n_samples){ 
  lines(density(y_sim[r,]), col = rgb(135, 206, 235, max= 255, alpha=25)) 
}
dev.off()

# Get P estimates --------------------------------------------------------------

P_samps <- model_samps$P
dim(P_samps)
P_samps.mean <- apply(P_samps, c(2,3,4,5), mean)
P_samps.mean <- plyr::adply(P_samps.mean, .margins=c(1,2,4))
colnames(P_samps.mean) <- c('index_subnat', 'index_method', 'index_year', 'Public', 'Private') 
P_samps.mean <- P_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Mean')

P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
dim(P_samps.quantile)
colnames(P_samps.quantile) <- c('index_subnat', 'index_method', 'index_sector', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) %>%
  left_join(index_country_subnat_tbl) %>%
  left_join(year_index_table)

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)
saveRDS(P_samps_df, file='results/STAN/P_samps/STAN_bspline_4country_Psamps_df.RDS')

# Get observed data 
P_df<- FP_source_data_wide %>%
  select(Country, Region, Method, index_year, average_year, Public, Private) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

P_SEdf<- FP_source_data_wide %>%
  select(Country, Region, Method, index_year, average_year, Public.SE, Private.SE) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public.SE, Private.SE), names_to = 'Sector', values_to = 'SE') %>%
  mutate(Sector = str_replace(Sector, '.SE', ''))

P_df <- left_join(P_df, P_SEdf) %>%
  rowwise() %>%
  mutate(lower_95 = Observed - 2*SE,
         upper_95 = Observed + 2*SE)

# Plot means vs observed values
for(c in n_country) {
  ggplot() +
    geom_point(data = P_df %>% filter(Country==c), aes(x=average_year, y=Observed, colour=Sector, pch=Sector)) +
    geom_errorbar(data = P_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax=upper_95, colour=Sector)) +
    geom_line(data = P_samps_df %>% filter(Country==c), aes(x=average_year, y=Mean, colour=Sector, lty=Sector)) +
    geom_ribbon(data = P_samps_df %>% filter(Country==c), aes(x=average_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
    theme(title = element_text(size=20), axis.text.x = element_text(angle = 90), strip.text.x = element_text(size=10), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) +
    theme(legend.position = "bottom", legend.title = element_text(size = 20), legend.text = element_text(size = 20))+
    facet_wrap(~interaction(Region, Method), ncol=5)
  ggsave(filename = paste0("p_plot_",c,".pdf"), path = vispath, height=12, width=15) 
  
}

