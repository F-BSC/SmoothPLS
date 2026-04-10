# tests/testthat/test-synthetic_data.R

test_that("generate_X_df correctly generates CFD data", {
  nind <- 10
  start <- 0
  end <- 100

  # Generation for categorical data
  df <- generate_X_df(nind = nind, start = start, end = end, curve_type = 'cat')

  expect_s3_class(df, "data.frame")
  expect_equal(length(unique(df$id)), nind)
  # Check time bounds
  expect_equal(min(df$time), start)
  expect_equal(max(df$time), end)
  # Check binary states
  expect_true(all(df$state %in% c(0, 1)))
})

test_that("generate_X_df correctly generates SFD data", {
  nind <- 5
  # Generation for scalar functional data
  df <- generate_X_df(nind = nind, curve_type = 'num', noise_sd = 0.1)

  expect_s3_class(df, "data.frame")
  expect_equal(length(unique(df$id)), nind)
  # For SFD, column should be 'value'
  expect_true("value" %in% names(df))
})

test_that("generate_Y_df calculates Y with correct dimensions and noise", {
  # Setup minimal CFD data
  df <- generate_X_df(nind = 20, curve_type = 'cat', seed = 123)

  # Use beta_1_real_func
  beta_f <- beta_1_real_func
  beta_0 <- 5.4321

  # 1. Test without noise (NotS_ratio = 0)
  Y_no_noise <- generate_Y_df(df, curve_type = 'cat',
                              beta_real_func_or_list = beta_f,
                              beta_0_real = beta_0, NotS_ratio = 0)

  expect_equal(nrow(Y_no_noise), 20)
  expect_equal(Y_no_noise$Y_real, Y_no_noise$Y_noised)

  # 2. Test with noise ratio
  Y_noisy <- generate_Y_df(df, curve_type = 'cat',
                           beta_real_func_or_list = beta_f,
                           beta_0_real = beta_0, NotS_ratio = 0.2)

  # Y_noised should differ from Y_real
  expect_false(identical(Y_noisy$Y_real, Y_noisy$Y_noised))
})

test_that("generate_X_df_multistates respects N_states", {
  N_states <- 3 # [cite: 880]
  lambdas <- lambda_determination(N_states)
  transition_df <- transfer_probabilities(N_states)

  df_multi <- generate_X_df_multistates(nind = 5, N_states = N_states,
                                        lambdas = lambdas,
                                        transition_df = transition_df)

  # Check that all states are within the expected range
  expect_true(all(df_multi$state >= 1 & df_multi$state <= N_states))
})

test_that("Synthetic data reproducibility via seed", {
  # Two calls with same seed should yield identical dataframes
  df1 <- generate_X_df(nind = 5, curve_type = 'cat', seed = 42)
  df2 <- generate_X_df(nind = 5, curve_type = 'cat', seed = 42)

  expect_identical(df1, df2)
})

test_that("beta functions return expected types and shapes", {
  t <- seq(0, 100, 10)
  # All provided real beta functions should handle vector inputs
  expect_type(beta_1_real_func(t), "double")
  expect_length(beta_1_real_func(t), length(t))
  expect_type(beta_2_real_func(t), "double")
  expect_length(beta_2_real_func(t), length(t))
  expect_type(beta_3_real_func(t), "double")
  expect_length(beta_3_real_func(t), length(t))
})

#### State_indicator_functions ####
test_that("state_indicator correctly transforms multi-state data", {
  # 1. Create a minimal multi-state individual
  # t=0: State A, t=10: State B, t=20: State A, t=30: End
  df_multi <- data.frame(
    id = 1,
    time = c(0, 10, 20, 30),
    state = c("A", "B", "A", "A")
  )

  res <- state_indicator(df_multi, id_col = 'id', time_col = 'time')

  # Should have columns: id, time, state_A, state_B
  expect_true(all(c("state_A", "state_B") %in% names(res)))

  # Check values for state_A
  # At t=0: state is A (1), t=10: state is B (0), t=20: state is A (1)
  expect_equal(res$state_A, c(1, 0, 1, 1))
  expect_equal(res$state_B, c(0, 1, 0, 0))
})

test_that("cat_data_to_indicator preserves timing and data structure", {
  # Create data with 3 states
  df_raw <- data.frame(
    id = rep(1, 5),
    time = c(0, 5, 15, 25, 30),
    state = c(1, 2, 3, 1, 1)
  )

  # Full processing pipeline
  processed_list <- cat_data_to_indicator(df_raw, id_col = 'id',
                                            time_col = 'time')

  # 1. Check list structure
  expect_type(processed_list, "list")
  expect_length(processed_list, 3) # 3 states
  expect_named(processed_list, c("state_1", "state_2", "state_3"))

  # 2. Check for duplicate removal (remove_duplicate_states)
  # For state_2: active only at t=5.
  # In the list, it should have the entry where it starts and where it ends.
  s2_df <- processed_list$state_2
  expect_equal(nrow(s2_df), 4) # t=0 (0), t=5 (1), t=15 (0), t=30 (0)

  # 3. Check individual dataframe integrity
  expect_true(all(s2_df$state %in% c(0, 1)))
})

test_that("cat_data_to_indicator handles character states correctly", {
  df_char <- data.frame(
    id = 1,
    time = c(0, 10, 20),
    state = c("rest", "run", "run")
  )

  res <- cat_data_to_indicator(df_char)
  expect_named(res, c("state_rest", "state_run"))
})

test_that("cat_data_to_indicator stops on invalid column number", {
  # Function expects exactly 3 columns: id, time, state
  df_bad <- data.frame(id = 1, time = 0, state = 1, extra = 99)
  expect_error(cat_data_to_indicator(df_bad),
               "The dataframe should have 3 columns")
})
