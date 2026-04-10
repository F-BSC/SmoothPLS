#### Smooth_PLS ####

# tests/testthat/test-smooth_PLS.R

test_that("smoothPLS univariate equivalence with discrete PLS (Theorem)", {
  # 1. Setup data
  nind <- 50
  df <- generate_X_df(nind = nind, curve_type = 'cat', seed = 123)
  Y_df <- generate_Y_df(df, curve_type = 'cat',
                        beta_real_func_or_list = beta_1_real_func, seed = 123)
  Y <- Y_df$Y_noised

  basis <- create_bspline_basis(0, 100, nbasis = 10)

  # 2. Run SmoothPLS
  spls_res <- smoothPLS(df_list = df, Y = Y, basis_obj = basis,
                        curve_type_obj = 'cat', orth_obj = TRUE,
                        print_steps = FALSE, plot_rmsep = FALSE,
                        print_nbComp = FALSE)

  # 3. Manual Discrete PLS on Lambda for comparison
  # We need the orthonormal basis to match the Equivalence Theorem assumptions
  ortho_basis <- gram_schmidt_orthonormalize(basis)
  Lambda <- evaluate_lambda(df, ortho_basis, curve_type = 'cat')

  # Discrete PLS model
  pls_discrete <- pls::plsr(Y ~ as.matrix(Lambda),
                            jackknife = TRUE, validation = 'LOO',
                            intercept = TRUE)

  # 4. Check Identity of scores (t_h)
  expect_equal(as.numeric(spls_res$plsr_model$scores),
               as.numeric(pls_discrete$scores),
               tolerance = 1e-5)
})

test_that("smoothPLS score orthogonality (Proposition 3 & 4)", {
  df <- generate_X_df(nind = 30, curve_type = 'cat', seed = 123)
  Y_df <- generate_Y_df(df, curve_type = 'cat',
                        beta_real_func_or_list = beta_2_real_func)

  basis <- create_bspline_basis(0, 100, nbasis = 10)
  spls_res <- smoothPLS(df_list = df, Y = Y_df$Y_noised, basis_obj = basis,
                        curve_type_obj = 'cat',
                        orth_obj = TRUE, plot_rmsep = FALSE,
                        print_nbComp = FALSE)

  scores <- spls_res$plsr_model$scores
  n_comp <- spls_res$nbCP_opti

  # If we have at least 2 components, they must be orthogonal
  if(n_comp >= 2) {
    cor_matrix <- cor(scores)
    diag(cor_matrix) <- 0 # Ignore self-correlation
    expect_true(all(abs(cor_matrix) < 1e-10))
  }
})

test_that("smoothPLS_predict consistency", {
  nind <- 20
  df <- generate_X_df(nind = nind, curve_type = 'cat', seed = 42)
  Y_df <- generate_Y_df(df, curve_type = 'cat',
                        beta_real_func_or_list = beta_1_real_func, seed = 42)

  basis <- create_bspline_basis(0, 100, nbasis = 10)
  spls_res <- smoothPLS(df_list = df, Y = Y_df$Y_noised, basis_obj = basis,
                        curve_type_obj = 'cat',
                        plot_rmsep = FALSE, print_nbComp = FALSE)

  # Prediction on the same training data
  Y_hat <- smoothPLS_predict(df_predict_list = df,
                             delta_list = spls_res$reg_obj,
                             curve_type_obj = 'cat')

  # Predictions should be identical to the fitted values of the underlying model
  fitted_pls <- as.numeric(spls_res$plsr_model$fitted.values[,,
                                                            spls_res$nbCP_opti])

  expect_equal(as.numeric(Y_hat), fitted_pls, tolerance = 1e-4)
})

test_that("smoothPLS handles multivariate mixed inputs", {
  nind <- 10
  df_cat <- generate_X_df(nind = nind, curve_type = 'cat')
  df_num <- generate_X_df(nind = nind, curve_type = 'num')
  Y <- rnorm(nind)

  basis <- create_bspline_basis(0, 100, nbasis = 5)

  # Testing the full multivariate pipeline
  expect_output(
    smoothPLS(df_list = list(df_cat, df_num), Y = Y,
              basis_obj = basis,
              curve_type_obj = list('cat', 'num'), plot_rmsep = FALSE,
              print_steps = FALSE, print_nbComp = TRUE,
              plot_reg_curves = FALSE)
  )
})
