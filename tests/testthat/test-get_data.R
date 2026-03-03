# Multicountry subnational dataset test
testthat::test_that("Get subnational data example", {
  data <- clean_fp_source_data(mcmsector::subnat_bivar_data,
                               deft_lookup = mcmsector::DEFT_DHS_database)
  
  inputs <- run_preprocessing_pipeline(raw_df = data, 
                                       area_classification = mcmsector::Country_classification,
                                       methods = c("Female Sterilization", "Implants", "Injectables", "IUD", "OC Pills"),
                                       all_years = seq(1990, 2030.5, by = 0.5),
                                       nseg = 10)
  testthat::expect_true(is.data.frame(data))
  testthat::expect_true(is.data.frame(inputs$data))
  testthat::expect_true(is.data.frame(inputs$T_star))
  testthat::expect_true(is.list(inputs$meta))
  testthat::expect_true(is.list(inputs$splines))
  
})
