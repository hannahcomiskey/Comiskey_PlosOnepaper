# Multicountry subnational dataset test
testthat::test_that("Get subnational model params example", {
  data <- clean_fp_source_data(mcmsector::subnat_bivar_data, deft_lookup = mcmsector::DEFT_DHS_database)
  
  inputs <- run_preprocessing_pipeline(raw_df = data, 
                                       area_classification = mcmsector::Country_and_area_classification_inclFP2020,
                                       methods = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills"),
                                       all_years = seq(1990, 2030.5, by = 0.5),
                                       nseg = 10)
  mydata <- inputs$data[,c("Public", "Public.SE")]
  logit.data <- mydata %>%
    dplyr::rowwise() %>%
    dplyr::mutate(logit.Public = log(Public/(1-Public)),
           logit.Public.Var = ((1/(Public*(1-Public)))^2)*Public.SE^2,
           logit.Public.SE = sqrt(logit.Public.Var))
  n_method = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills")
  M =  length(n_method)
  inputdata <- list(y = as.vector(unlist(logit.data[,c("logit.Public")])), # using total proportions as collapsing over sectors
                    se_prop = as.vector(unlist(logit.data[,c("logit.Public.SE")])),
                    Bik = inputs$splines$B_ik,
                    n_years = inputs$meta$n_years,
                    n_obs = nrow(logit.data),
                    K = inputs$splines$K,
                    H = inputs$splines$H,
                    kstar = inputs$splines$Kstar,
                    C_count = length(inputs$meta$n_country),
                    P_count = length(inputs$meta$n_subnat),
                    M_count = length(c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills")),
                    Omega = diag(M) * 0.01 + diag(M),
                    matchcountry= inputs$meta$index_country_subnat_tbl$index_country,
                    matchsubnat = inputs$data$index_subnat,
                    matchmethod = inputs$data$index_method,
                    matchyears = inputs$data$index_year)
  pars <- c( "alpha_pms", # required for P
             "inv.Sigma.alpha_cms",
             "inv.Sigma.alpha_pms",
             "sigma_delta",
             "beta.k")
  
  model_path <- system.file("model", "ncp_MVN_SEmodel_JAGS.txt", package = "mcmsector")
  
  mod <- R2jags::jags.parallel(data=inputdata,
                       parameters.to.save=pars,
                       model.file = model_path,
                       n.iter = 10,         # total number of iterations per chain
                       n.burnin = 1,
                       n.thin=3)
  
  method_index_table <- tibble::tibble(index_method = 1:length(n_method), Method = n_method)
  sector_index_table <- tibble::tibble(index_sector = 1:2, Sector = c('Public', 'Private'))
  
  all_years = seq(1990, 2030.5, by = 0.5)
  
  year_index_table <- tibble::tibble(average_year = all_years, index_year = 1:length(all_years))
  subnat_index_table <- inputs$data %>% dplyr::select(Country, Region, index_country, index_subnat) %>% dplyr::distinct()
  
  params <- extract_mcmsector_params(mod = mod, 
                                     n_method = length(n_method),
                                     n_subnat = length(inputs$meta$n_subnat))
  testthat::expect_true(is.array(params$alpha_pms))
  testthat::expect_true(is.array(params$beta_k))
})  
