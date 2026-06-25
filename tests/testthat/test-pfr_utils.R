
# Tests for lambda_augmented_uni

test_that("lambda_augmented_uni correctly adds the intercept", {

  # 1. Standard test with a matrix
  lam_matrix <- matrix(c(2, 4, 6, 8), nrow = 2)
  res_matrix <- lambda_augmented_uni(lam_matrix)

  expect_equal(ncol(res_matrix), ncol(lam_matrix) + 1) # Check column addition
  #expect_equal(colnames(res_matrix)[1], "intercept")   # Check column name
  expect_true(all(res_matrix[, 1] == 1))     # Check if filled with 1s
  expect_equal(res_matrix[, 2], lam_matrix[, 1])       # Check if original data is intact

  # 2. Test with a simple vector (common edge case)
  lam_vec <- c(10, 20, 30)
  res_vec <- lambda_augmented_uni(lam_vec)

  expect_equal(dim(res_vec), c(3, 2)) # A vector of size 3 should become a 3x2 matrix
  #expect_equal(colnames(res_vec)[1], "intercept")
})


# Tests for calc_penalty_matrix


test_that("calc_penalty_matrix computes the penalty matrix correctly", {

  skip_if_not_installed("fda")

  nbasis <- 5
  test_basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)

  # 1. Test with default value (LDO = 2) on a basis.fd object
  pen_mat <- calc_penalty_matrix(test_basis)

  expect_true(is.matrix(pen_mat))
  expect_equal(dim(pen_mat), c(nbasis, nbasis))
  expect_true(isSymmetric(unname(pen_mat)))

  # 2. Explicit test with LDO = 2
  pen_mat_explicit <- calc_penalty_matrix(test_basis, LDO = 2)
  expect_equal(pen_mat, pen_mat_explicit)

  # 3. Test error messages
  expect_error(
    calc_penalty_matrix(test_basis, LDO = "D2"),
    "LDO must be a numeric integer representing the derivative order."
  )

  # 4. Test with a list of 'fd' objects
  coef1 <- matrix(c(1, 0, 0, 0, 0), ncol = 1)
  coef2 <- matrix(c(0, 1, 0, 0, 0), ncol = 1)
  fd1 <- fda::fd(coef1, test_basis)
  fd2 <- fda::fd(coef2, test_basis)

  fd_list <- list(fd1, fd2)
  pen_mat_list <- calc_penalty_matrix(fd_list, LDO = 2)

  expect_true(is.matrix(pen_mat_list))
  expect_equal(dim(pen_mat_list), c(2, 2))
  expect_true(isSymmetric(pen_mat_list))

  # 5. Test invalid input
  expect_error(
    calc_penalty_matrix(list(1, 2, 3), LDO = 2),
    "basis_obj must be a 'basisfd' object, an 'fd' object, or a list"
  )

})



# Tests for augment_penalty_matrix

test_that("augment_penalty_matrix adds zero padding correctly", {

  # 1. Create a dummy 2x2 penalty matrix
  orig_mat <- matrix(c(4, -1, -1, 4), nrow = 2, ncol = 2)

  # 2. Augment the matrix
  aug_mat <- augment_penalty_matrix(orig_mat)

  # 3. Check new dimensions (should be 3x3)
  expect_equal(nrow(aug_mat), nrow(orig_mat) + 1)
  expect_equal(ncol(aug_mat), ncol(orig_mat) + 1)

  # 4. Check that the first row and first column are purely zeros
  expect_true(all(aug_mat[1, ] == 0))
  expect_true(all(aug_mat[, 1] == 0))

  # 5. Check that the original matrix is perfectly preserved in the bottom-right
  expect_equal(aug_mat[-1, -1], orig_mat, ignore_attr = TRUE)

  # 6. Check the naming
  #expect_equal(colnames(aug_mat)[1], "intercept")
  #expect_equal(rownames(aug_mat)[1], "intercept")
})


# Tests for fit_pfr_uni and cv_pfr_uni

test_that("fit_pfr_uni and cv_pfr_uni work properly on simulated data", {

  # 1. Simulate dummy data
  set.seed(42)
  n_obs <- 30
  n_basis <- 5

  # X is 30 x 6 (1 intercept + 5 basis evaluations)
  X_dummy <- cbind(1, matrix(rnorm(n_obs * n_basis), nrow = n_obs))

  # Y is just a random vector
  Y_dummy <- rnorm(n_obs)

  # R is the penalty matrix (6x6, augmented with 0 for intercept)
  R_dummy <- diag(c(0, rep(1, n_basis)))

  # --- Test fit_pfr_uni ---
  fit_res <- fit_pfr_uni(X_dummy, Y_dummy, R_dummy, lambda = 10)

  expect_type(fit_res, "list")
  expect_length(fit_res$beta_hat, n_basis + 1)
  expect_length(fit_res$Y_hat, n_obs)
  expect_equal(length(fit_res$residuals), n_obs)

  # --- Test cv_pfr_uni ---
  lam_grid <- 10^seq(-2, 2, length.out = 10)
  cv_res <- cv_pfr_uni(X_dummy, Y_dummy, R_dummy, lambda_grid = lam_grid)

  expect_s3_class(cv_res, "cv_pfr_uni")
  expect_true(cv_res$optimal_lambda %in% lam_grid)
  expect_true(is.numeric(cv_res$min_rmsep))
  expect_equal(nrow(cv_res$cv_results), length(lam_grid))
})
