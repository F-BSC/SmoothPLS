test_that("funcPLS execution and consistency", {
  nind <- 30
  df <- generate_X_df(nind = nind, curve_type = 'cat', seed = 123)
  Y <- generate_Y_df(df, curve_type = 'cat',
                     beta_real_func_or_list = beta_1_real_func,
                     seed = 123)$Y_noised
  basis <- create_bspline_basis(0, 100, nbasis = 10)

  # Run FPLS
  res_fpls <- funcPLS(df_list = list(df), Y = Y, basis_obj = basis,
                      regul_time_obj = seq(0, 100, 1), curve_type_obj = 'cat',
                      print_nbComp = FALSE, plot_rmsep = FALSE,
                      print_steps = FALSE, plot_reg_curves = FALSE)

  expect_s3_class(res_fpls$plsr_model, "mvr")
  expect_equal(length(res_fpls$reg_obj), 2) # Intercept + 1 curve
  expect_s3_class(res_fpls$reg_obj[[2]], "fd")
})
