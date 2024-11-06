library(rstan)
library(tidyverse)
library(tidybayes)
library(bayesplot)

# Source simulated data --------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1, auto_write = TRUE)

load("data/simulated_data_new.RData")

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
fit <- readRDS('results/STAN_model_test_spline_simdata.RDS')

traceplot(fit, pars = c("alpha_pms[3,1]", "alpha_pms[3,2]", "alpha_pms[3,3]", "alpha_pms[3,4]", "alpha_pms[3,5]", "alpha_pms[3,6]"), inc_warmup = FALSE, nrow = 6)
traceplot(fit, pars = c("sigmabeta_Omega[1,1]",
                        "sigmabeta_Omega[2,2]",
                        "sigmabeta_Omega[3,3]",
                        "sigmabeta_Omega[4,4]",
                        "sigmabeta_Omega[5,5]"), inc_warmup = FALSE, nrow = 5)
traceplot(fit, pars = c("delta_k[1,1,1]",  "delta_k[2,3,5]",  "delta_k[3,2,2]", "delta_k[4,5,7]", "delta_k[5,4,6]"), inc_warmup = FALSE, nrow = 6)


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

plot(fit, pars=c('L_Sigma_beta'))
plot(fit, pars=c('L_Sigma_delta'))

pairs(fit,pars=c("beta_c"))

mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))

np_cp <- nuts_params(fit)
mcmc_nuts_energy(np_cp)

fit_summary <- summary(fit)
View(print(fit_summary$summary))

# Get parameter estimates 

model_samps <- rstan::extract(fit)

# Check delta.k sums to 0 
delta.k_samps <- model_samps$a_k
dim(delta.k_samps)
sum(delta.k_samps[1,1,1,1:10])

# Review spline coefficients
delta_k.mean <- apply(delta.k_samps, c(2,3,4), mean)
knots.k <- c(4.625,8.250,11.875,15.500,19.125,22.750,26.375)

## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,B.ik[,1], type= "n", xaxt="n",
     xlab = "Year",
     ylab ="Basis Function",
     xlim = range(all_years),
     ylim=c(0,1))
axis(1, at = min(all_years):max(all_years))
abline(v=knots.k, col = seq(1, K), lwd = 1)
for (k in 1:H){
  lines(all_years,B.ik[,k], type= "l", col = k, lwd = 1)
}

random_spline_comp <- B.ik %*% delta_k.mean[2,3,]
## Plot basis
par(lwd = 3, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3, mfrow = c(1,1))
plot(all_years,random_spline_comp[,1], type= "n", xaxt="n",
     xlab = "Year", ylab ="Basis Function",
     xlim = range(all_years))
axis(1, at = min(all_years):max(all_years))
abline(v=knots.k, col = seq(1, H), lwd = 1)
lines(all_years, random_spline_comp, type= "l", col = k, lwd = 1)


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
  mutate(index_method = rep(1:5, each=20))
alpha_samps.mean <- alpha_samps.mean %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  mutate(invlogit.alpha = exp(alpha_mean)/(1+exp(alpha_mean)))

# Get observed data 
P_df<- P_sim_df_sample %>% 
  mutate(across(everything(), as.numeric)) %>%
  left_join(method_index_table) %>%
  pivot_longer(cols = c(Public, Private), names_to = 'Sector', values_to = 'Observed')

subnatid = 1
# Plot means vs observed values
ggplot() +
  geom_point(data = P_df %>% filter(index_subnat==subnatid), aes(x=index_year, y=Observed, colour=Sector, pch=Sector)) +
  geom_line(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=index_year, y=Mean, colour=Sector, lty=Sector)) +
  geom_ribbon(data = P_samps_df %>% filter(index_subnat==subnatid), aes(x=index_year, ymin=lower_95, ymax = upper_95, fill=Sector), alpha=0.2) +
  geom_hline(data = alpha_samps.mean %>% filter(index_subnat==subnatid), aes(yintercept = invlogit.alpha)) +
  facet_wrap(~Method)



