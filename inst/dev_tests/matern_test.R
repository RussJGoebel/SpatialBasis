library(downscaling)

test_regions <- build_region_list(target_grid)

n_basis <- c(6000,6000,7000,8000,9000,10000)

results <- matrix(nrow = length(test_regions),ncol = length(n_basis))

for(jj in seq_along(n_basis)){

  basis <- make_matern_fourier_basis(max_frequency = 1000, K = n_basis[jj])
  prior <- make_matern_fourier_prior(basis)
  spat_field <- spatial_field(basis,prior,integration_rule_qmc_average(generate_qmc_unit_square(256)))

  prep <- prepare_spatial_data(spat_field,soundings,"SIF_757nm", parallel = TRUE, n_workers = 32)

  fit <- fit_spatial_field(prep)

  test_region_design <- integrate_region_list_parallel(test_regions,rule = integration_rule_qmc_average(generate_qmc_unit_square(256)), basis = basis,n_workers = 32)

  results[,jj] <- test_region_design %*% fit$posterior_mean

  message("Done with ",jj, "basis function test")

}
