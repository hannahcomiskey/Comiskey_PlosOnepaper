library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
source("code/load_functions.R")
source("code/read_in_IPUMS_subnational_SEdata.R")

# simulate data ----------------------------------------------
## all years including projection period
all_years <- seq(from = 1990, to = 2030.5, by=0.5)
n_all_years <- length(all_years)

## obs years
n_years = 5
year = sort(sample(all_years[1:(n_all_years - 15)], size = n_years))
x <- match(year,all_years) # index for obs years

subset_C = "Nepal" 
subset_R = "Central"

# Get data ---------
dat <- FP_source_data_wide %>% 
  filter(Country==subset_C & Region==subset_R) %>% 
  select(Country, Method, Region, average_year, Public, Public.SE) %>%
  arrange(Region)

n_method <- tibble(Region = unique(dat$Method), index_method = 1:length(unique(dat$Method)))
n_subnat <- tibble(Region = unique(dat$Region), index_subnat = 1:length(unique(dat$Region)))
year_index <- tibble(average_year = all_years, index_year = 1:length(all_years))
dat <- dat %>%
  left_join(n_subnat) %>% # Region indexing
  left_join(n_method)

dat <- dat %>% 
  mutate(logit.pub = log(Public/(1-Public))) %>% 
  mutate(logit.pubSE = log(Public.SE/(1-Public.SE))) %>%
  mutate(logit.ymin = logit.pub - 2*logit.pubSE) %>%
  mutate(logit.ymax = logit.pub + 2*logit.pubSE) %>%
  mutate(ymin = Public - 2*Public.SE) %>%
  mutate(ymax = Public + 2*Public.SE) %>%
  mutate(index_year = match(average_year, all_years))

# Plot data 
ggplot() +
  geom_line(data = dat, aes(x=average_year, y=Public, group=Region)) +
  facet_wrap(~Method)

## JAGS model for spline (forward/backward approach)
fit_spline1 <- '
model{

 # Variance structure
 inv.sigma_delta[1:M_count,1:M_count] ~ dwish(natdf*natRmat[1:M_count,1:M_count],natdf)
 
 for(m in 1:M_count) {
  for(t in 1:n_all_years) {
    logit(mu[m,t]) <- alpha[m] + inprod(B.ik[t,],beta.k[m,])
  }
  ## Beta = 0 at Tstar.... => delta(Tstar+1) = B(Tstar+1), and -delta(Tstar-1)=B(Tstar-1)
    for(j in 1:(K-1)){
      delta.k[m,j] ~ dmnorm(mu_delta,inv.sigma_delta) # change between spline coeff for country c and method m
    } # end H loop
  
    beta.k[m,Kstar] <- 0

    for(i in 1:(Kstar-1)) { # before Kstar
      beta.k[m,(Kstar-i)] <- beta.k[m,(Kstar-i)+1] - delta.k[m,(Kstar-i)]
    }
    for(i in (Kstar+1):K) { # After Kstar+1
      beta.k[m,i] <- beta.k[m,i-1] + delta.k[m,i-1]
    }
  } # end m loop 
  
  # Likelihood
  for (k in 1:n_obs) {
    y[k] ~ dnorm(mu[matchmethod[k],matchyear[k]], tau_y[k])T(0,1)
    tau_y[k] <- 1/(se_prop[k]^2) 
    }
  
  ## priors
  for(m in 1:M_count) {
    alpha[m] ~ dnorm(alpha.mean[m],0.1)
  }
} 
'

######################################################################
# Read in national alpha estimates ---------------------------------
######################################################################
alpha_med <- readRDS(file="data/alpha_estimates/mean_alpha_intercepts.RDS")
alpha_sd <- readRDS(file="data/alpha_estimates/sd_alpha_intercepts.RDS") 

######################################################################
# Read in country indexing ---------------------------------
######################################################################
national_country_index_table <- readRDS("data/alpha_estimates/country_index_table.RDS")

######################################################################
# Get group size for each sector, method combination (UIP) -----------
######################################################################
n_sector <- c("Public", "Commercial_medical", "Other") # Names of categories
alpha_groupsize <- readRDS("data/alpha_estimates/Bspline_correlations_groupsize.RDS") %>%
  pivot_wider(names_from = "Method", values_from = "n") %>%
  arrange(factor(sector_category, levels = n_sector)) # order sectors properly

######################################################################
# Match intercepts to country names  ---------------------------------
######################################################################
dimnames(alpha_med)[[3]] <- dimnames(alpha_sd)[[3]] <- as.list(unlist(national_country_index_table$Country))
myalpha_med <- alpha_med[,,subset_C]
myalpha_sd <- alpha_sd[,,subset_C]

# Get T_star and match_Tstar-----------------------------------------------
T_star <- dat %>%
  group_by(Region) %>%
  filter(index_year==max(index_year)) %>%
  select(Region, index_subnat, average_year, index_year) %>%
  arrange(index_subnat) %>%
  ungroup() %>%
  select(index_subnat, index_year, average_year) %>%
  distinct()

# splines -----------------------------------------------------------------
nseg=12
Kstar <- vector()
B.ik <- array(dim = c(length(n_subnat$Region), length(all_years),nseg+3))
knots.all <- matrix(nrow = length(n_subnat$Region), ncol=nseg+3)

for(r in 1:nrow(T_star)) {
  index_mc <- T_star$average_year[r]
  res <- bs_bbase_precise(all_years, lastobs=index_mc, nseg = nseg) # number of splines based on segments, here choosing 10
  B.ik[r,,] <- res$B.ik
  Kstar[r] <- res$Kstar
  knots.all[r,] <- res$knots.k
}

K <- dim(res$B.ik)[2]
H <- K-1

jags_data1 <- list(y = dat$Public,
                   sigma_y = dat$Public.SE,
                   alpha.mean = myalpha_med[1,],
                   Kstar = Kstar,
                   B.ik = B.ik,
                   K = K,
                   n_all_years = n_all_years,
                   n_obs = nrow(dat),
                   mu_delta = rep(0, length(n_subnat$Method)),
                   M_count = length(n_subnat$Method),
                   R_count = length(n_subnat$Region),
                   matchmethod = dat$index_method,
                   matchyear = dat$index_year)

# parameters to save
jags_pars <- c("mu",
               "alpha",
               "beta.k",
               "sigma_delta",
               "delta.k")

# run model
mod1 <- jags(data = jags_data1,
            parameters.to.save=jags_pars,
            model.file = textConnection(fit_spline1),
            n.iter = 20000,
            n.burnin = 10000,
            n.thin = 5)

plot(mod1)

##create an object containing the posterior samples
m1<- mod1$BUGSoutput$sims.matrix
post_samples <- as_tibble(m1)

## plot alpha vs the truth
ggplot(post_samples, aes(x = `alpha[5]`)) +
  stat_halfeye() +
  geom_vline(mapping = aes(xintercept = alpha.inj,  colour = "True value"))

## format data for plotting results
sample_draws <- tidy_draws(m1)
n_iter <- nrow(sample_draws)
pred_summary1 <- sample_draws %>%
  select(`mu[1,1]`:`mu[8,81]`) %>%
  tidyr::pivot_longer(`mu[1,1]`:`mu[8,81]`,
                      names_to = "n",
                      values_to = "mu_pred") %>%
  dplyr::mutate(index_year = rep(rep(1:n_all_years, each=8), n_iter)) %>%
  dplyr::mutate(index_subnat = rep(rep(1:8, n_all_years), n_iter)) %>%
  dplyr::group_by(index_year, index_subnat) %>%
  tidybayes::median_qi(mu_pred) %>%
  left_join(year_index) %>%
  left_join(n_subnat) %>%
  dplyr::mutate(
    pred_mu = mu_pred,
    lwr_95 = .lower,
    upr_95 = .upper) %>%
  dplyr::select(average_year, Region, pred_mu, lwr_95, upr_95)

pred_summary1 <- pred_summary1 %>% 
  mutate(invlogit.pred_mu = exp(pred_mu)/(exp(pred_mu)+1)) %>%
  mutate(invlogit.lwr_95 = exp(lwr_95)/(exp(lwr_95)+1)) %>%  
  mutate(invlogit.upr_95 = exp(upr_95)/(exp(upr_95)+1)) 

# ALpha data and estimates
alpha_est <- mod1$BUGSoutput$median$alpha
invlogit.alpha <- exp(alpha_est)/(1+exp(alpha_est))

# Splines 
basis_data <- as_tibble(B.ik[1,,]) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "spline", values_to = "response") %>% 
  #mutate(response = response-1) %>%
  mutate(spline = gsub("V","", spline)) %>%
  mutate(spline = as.numeric(spline)) %>%
  arrange(spline) %>%
  mutate(index_year = rep(1:length(all_years), ncol(B.ik[1,,]))) %>%
  mutate(year = rep(all_years, ncol(B.ik[1,,]))) %>%
  mutate(prop.response = exp(response)/(1+exp(response)))

## plot obs + model estimates and CIs
ggplot() +
  geom_line(data = pred_summary1, aes(x = average_year, y = invlogit.pred_mu), colour="Red") +
  geom_ribbon(data = pred_summary1, aes(x = average_year, ymin =  invlogit.lwr_95, ymax =  invlogit.upr_95), fill="Red", alpha = 0.2) +
  geom_point(data = dat, aes(x = average_year, y = Public)) + #here
  geom_errorbar(data= dat, aes(x=average_year, ymin = ymin, ymax=ymax), width = 1.5) +
  #geom_hline(yintercept =  invlogit.alpha, lty=3) +
  #scale_x_continuous(breaks = 1990.5:2030.5) +
  scale_y_continuous(breaks = seq(from= 0, to =1, by=0.1)) +
  facet_wrap(~Region)
  
## look at deltas - transform back to spline coefficients
del.k<-mod1$BUGSoutput$median$delta.k
del.k

midpoints <- res$knots.k[-length(res$knots.k)] + diff(res$knots.k) / 2

deltas_data <- cbind(year = res$knots.k, delta =c(0,del.k))

beta.k<-mod1$BUGSoutput$median$beta.k
betas_data <- cbind(year = res$knots.k, beta=beta.k)

B.ik[which(all_years==2019.5),]

## plot deltas
ggplot() +
  geom_line(data = as_tibble(betas_data), aes(x = year, y = beta), colour="Magenta") +
  geom_line(data = as_tibble(deltas_data), aes(x = year, y = delta), colour="Red") +
  geom_line(data = pred_summary1, aes(x = year, y = invlogit.pred_mu), colour="Blue") +
  geom_vline(xintercept=res$knots.k) +
  geom_point(data = dat, aes(x = average_year, y = Public)) + #here
  geom_errorbar(data= dat, aes(x=average_year, ymin = ymin, ymax=ymax), width = 1.5) +
  geom_hline(yintercept =  invlogit.alpha, lty=3) +
  geom_line(data=basis_data, aes(x=year, y=response-1, group=spline), show.legend = FALSE)

# Explanation of bump in 2020 --------------------
inv.logit(alpha_est) # 2019.5
inv.logit(alpha_est + sum(B.ik[(which(all_years==2020.5)),]*beta.k)) 
inv.logit(alpha_est + sum(B.ik[(which(all_years==2021.5)),]*beta.k)) 
inv.logit(alpha_est + sum(B.ik[(which(all_years==2022.5)),]*beta.k)) # Three years before return to tstar level

  
  
  