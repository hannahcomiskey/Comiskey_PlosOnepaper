library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")
source('code/stan_utility.R')

fit <- readRDS('results/STAN_model_KenyaCameroon_allN_NCP_kstar_sum0.RDS')

code <- get_stancode(fit)
cat(code)

check_all_diagnostics(fit)

model_samps <- rstan::extract(fit)

# Traceplots
traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("beta_c"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("beta_k[1,1,1]",  "beta_k[2,3,5]",  "beta_k[3,2,2]", "beta_k[4,5,7]", "beta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_beta"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)

check_divergences(fit)

# Check divergences 
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params$`beta_c[1,2]`, nondiv_params$`sigma_beta[1]`,
     col=c_dark, pch=16, cex=0.8)
points(div_params$`beta_c[1,1]`, div_params$`sigma_beta[1,1]`,
       col=green, pch=16, cex=0.8)

# Checking for funnel in beta parameters
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("beta_c[1,", k, "]", sep="")
  name2 <- paste("sigma_beta[", k, "]", sep="")
  plot(nondiv_params[name][[1]], log(nondiv_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name, ylab="log(tau)")
  points(div_params[name][[1]], log(div_params[name2][[1]]),
         col="blue", pch=16, cex=0.8)
}
par(mfrow=c(1, 1))


# Plot beta spline coefficients
betak_samps <- model_samps$beta_k
dim(betak_samps)
betak_samps.mean <- apply(betak_samps, c(2,3,4), mean)
dim(betak_samps.mean)

betak_test <- tibble(Beta_k = betak_samps.mean[1,3,], H=1:13)

ggplot() +
  geom_line(data = betak_test, aes(x=H, y=Beta_k))

# Set up indexing --------------------------------------------------------------
method_index_table <- tibble(index_method = 1:length(n_method), Method = n_method)
sector_index_table <- tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
subnat_index_table <- cbind(country_subnat_tbl, index_subnat=1:10)

# Check priors vs posteriors samples -------------------------------------------
fit %>%
  spread_draws(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`,
               `beta_c[2,1]`, `beta_c[2,2]`, `beta_c[2,3]`, `beta_c[2,4]`, `beta_c[2,5]`) %>%
  pivot_longer(cols = c(`beta_c[1,1]`, `beta_c[1,2]`, `beta_c[1,3]`, `beta_c[1,4]`, `beta_c[1,5]`,
                        `beta_c[2,1]`, `beta_c[2,2]`, `beta_c[2,3]`, `beta_c[2,4]`, `beta_c[2,5]`), names_to='Parameter', values_to = 'sample') %>%
  ggplot() +
  geom_density(aes(x = sample), fill='grey', alpha=0.8) +
  facet_wrap(~Parameter, ncol=5)
ggsave(filename = "betac_plot_KenyaCameroon.pdf", path = 'visualisations/country_data/allN/', height=12, width=15) 


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

P_samps <- model_samps$P
dim(P_samps)
P_samps.quantile <- apply(P_samps, c(2,3,4,5), quantile, probs=c(0.025, 0.975), na.rm=TRUE)
P_samps.quantile <- plyr::adply(P_samps.quantile, .margins=c(2,3,4,5))
colnames(P_samps.quantile) <- c('index_method', 'index_subnat','index_sector', 'index_year', 'lower_95', 'upper_95') 
P_samps.quantile <- P_samps.quantile %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  left_join(sector_index_table) %>%
  left_join(subnat_index_table)

P_samps_df <- left_join(P_samps.mean, P_samps.quantile)


# Get observed data 
P_df<- FP_source_data_wide %>%
  select(!c(Public.SE, Public_n, Private.SE, Private_n)) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

#%>% filter(index_year<30)
# Plot means vs observed values
ggplot() +
  geom_point(data = P_df, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df , aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, Region, Country), ncol=5)
ggsave(filename = "p_plot_KenyaCameroon.pdf", path = 'visualisations/country_data/allN/', height=12, width=15) 


