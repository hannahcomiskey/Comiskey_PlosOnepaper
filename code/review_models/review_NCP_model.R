library(rstan)
library(tidyverse)
library(tidybayes)

# Source code --------------------------------------
source("code/load_functions.R")
source("code/2sector_code/read_in_subnational_2sector_data.R")
source("code/2sector_code/set_up_2sector_bivar_globalrunjags.R")
source('code/stan_utility.R')

fit <- readRDS('results/STAN_model_kstar_reg_NCP_allN_Kenya.RDS')

code <- get_stancode(fit)
cat(code)

check_all_diagnostics(fit)

model_samps <- rstan::extract(fit)

traceplot(fit, pars = c("alpha_pms"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("alpha_raw"), inc_warmup = FALSE, nrow = 6)

traceplot(fit, pars = c("sigma_alpha"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("sigma_delta"), inc_warmup = FALSE, nrow = 5)

check_divergences(fit)

# Check divergences
c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(div_params$`delta_k[1,1,1]`, log(div_params$`lstar_sq[1,1,1]`),
       col='blue', pch=16, cex=0.8)
points(nondiv_params$`delta_k[1,1,1]`, log(nondiv_params$`lstar_sq[1,1,1]`),
       col=c_dark, pch=16, cex=0.8)

# Checking for funnel in beta parameters
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("beta_c[1,", k, "]", sep="")
  name2 <- paste("sigma_beta[", k, "]", sep="")
  plot(div_params[name][[1]], log(div_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name, ylab="log(tau)")
  points(nondiv_params[name][[1]], log(nondiv_params[name2][[1]]),
         col="blue", pch=16, cex=0.8)
}
par(mfrow=c(1, 1))

# Checking for funnel in beta parameters
par(mfrow=c(5, 1))
for (k in 1:length(n_method)) {
  name <- paste("alpha_pms[1,", k, "]", sep="")
  name2 <- paste("sigma_alpha[", k, "]", sep="")
  plot(div_params[name][[1]], log(div_params[name2][[1]]),
       col="#8F272780", pch=16, cex=0.8,
       xlab=name, ylab="log(tau)")
  points(nondiv_params[name][[1]], log(nondiv_params[name2][[1]]),
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

# Get P estimates
subnat_index_table <- FP_source_data_wide %>% ungroup() %>% select(Region, index_subnat) %>% distinct() 
sector_index_table <- tibble(Sector = c('Public', 'Private'), index_sector = 1:2)
method_index_table <- tibble(Method= n_method, index_method = 1:length(n_method))
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
plot_data <- FP_source_data_wide %>%
  select(Country, Region, Method, average_year, index_year, Public, Private) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

plot_SE <- FP_source_data_wide %>%
  select(Country, Region, Method, average_year, index_year, Public.SE, Private.SE) %>%
  pivot_longer(cols = c(Public.SE, Private.SE), names_to = 'Sector', values_to = 'SE') %>%
  mutate(Sector = str_replace(Sector, '.SE', ''))

plot_data <- plot_data %>% left_join(plot_SE)

# Plot means vs observed values
ggplot() +
  geom_point(data = plot_data, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_errorbar(data = plot_data, aes(x=index_year, ymin = Observed-2*SE, ymax=Observed+2*SE,  colour=Sector))+
  geom_line(data = P_samps_df, aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df, aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  facet_wrap(~interaction(Method, Region), ncol=5)


alpha_samps <- model_samps$alpha_pms
dim(alpha_samps)
alpha_mean <- apply(alpha_samps, c(2,3), mean)
colnames(alpha_mean) <- subnat_index_table$Region

alpha_mean <- as_tibble(alpha_mean) %>%
  mutate(index_method = 1:length(n_method)) %>%
  pivot_longer(cols=subnat_index_table$Region, names_to = 'Region', values_to = 'alpha') %>%
  mutate(invlogit.alpha = exp(alpha)/(1+exp(alpha))) %>%
  left_join(method_index_table)

ggplot() +
  geom_point(data = plot_data, aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_errorbar(data = plot_data, aes(x=index_year, ymin = Observed-2*SE, ymax=Observed+2*SE,  colour=Sector))+
  geom_hline(data=alpha_mean, aes(yintercept = invlogit.alpha)) +
  facet_wrap(~interaction(Method, Region), ncol = 5)

