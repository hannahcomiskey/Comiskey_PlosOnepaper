library(rstan)
library(shinystan)
library(tidyverse)
library(tidybayes)
library(bayesplot)
source('code/stan_utility.R')

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

load("data/simulated_data/simulated_data_bspline.RData")

# Get logit of parameters and variance -----------------------------------------
mydata <- P_sim_df_sample[,c("Public")] %>%
  rowwise() %>%
  mutate(Public = ifelse(Public < 0.0001 , 0, Public))

logit.data <- mydata %>%
  rowwise() %>%
  mutate(logit.Public = log(Public/(1-Public)))

# # testing splines ------------------------------------------------------------
X <- seq(1, 20.5, by=0.5) # generating inputs
B <- t(splines::bs(X, degree=3, knots=c(seq(1, 20, by=4)),  intercept = TRUE)) # creating the B-splines
num_data <- length(X)
num_basis <- nrow(B)

# Set up model inputs ----------------------------------------------------------
simmatchsubnat <- as.vector(as.numeric(P_sim_df_sample$index_subnat))
simmatchmethod <- as.vector(as.numeric(P_sim_df_sample$index_method))
simmatchyears <- as.vector(as.numeric(P_sim_df_sample$index_year))
n_all_years <- length(X)
M_count = 5

## The required data ------------------------------
inputdata <- list(Y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                  X = X,
                  num_years = length(X),
                  num_data = nrow(logit.data),
                  num_knots = nrow(B),
                  P_count = P,
                  M_count = M,
                  S_count = 2,
                  B = B,
                  matchsubnat = simmatchsubnat,
                  matchmethod = simmatchmethod,
                  matchyears = simmatchyears
                  )

## Parameters to look at ------------------------------
pars <- c("Y_hat",
          "a0",
          "a",
          "tau",
          "P")

# Run stan model ------------------

fit <- stan(
  data = inputdata,    # named list of data
  file = 'model/STAN/Bspline_model.stan',
  pars = pars,
  iter = 10000,         # total number of iterations per chain
  warmup = 2000,
  thin=4,
  chains=1,
  save_warmup = FALSE,
  control=list(adapt_delta=0.99, max_treedepth=12)
)

saveRDS(fit, 'results/model_testing_Bspline.RData')

# fit <- readRDS('results/model_testing_Bspline.RData')
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
# plot(nondiv_params$`alpha_pms[3,3]`, log(nondiv_params$`sigma_alpha[3]`),
#      col=c_dark, pch=16, cex=0.8)
# points(div_params$`alpha_pms[3,3]`, log(div_params$`sigma_alpha[3]`),
#        col=green, pch=16, cex=0.8)
# 
# 
# par(mar = c(4, 4, 0.5, 0.5))
# plot(nondiv_params$`delta_k[3,3,6]`, log(nondiv_params$`sigma_delta[3]`),
#      col=c_dark, pch=16, cex=0.8)
# points(div_params$`delta_k[3,3,6]`, log(div_params$`sigma_delta[3]`),
#        col=green, pch=16, cex=0.8)
# 
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
