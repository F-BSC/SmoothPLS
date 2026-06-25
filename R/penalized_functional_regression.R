#### Utils ####

#' Augment a lambda matrix with an intercept column
#'
#' @description
#' Prepends a column of ones to a given matrix, effectively adding an
#' intercept term for regression models.
#'
#' @param lambda A numeric matrix or data frame containing the lambda values.
#'
#' @return A matrix with an additional first column named "intercept"
#' filled with ones.
#' @export
#'
#' @author Francois Bassac
#'
#' @examples
#' lam <- matrix(1:4, nrow = 2)
#' lambda_augmented_uni(lam)
lambda_augmented_uni <- function(lambda) {

  lambda_matrix <- as.matrix(lambda)

  one_vector_lambda <- cbind(1, lambda_matrix)
  #colnames(one_vector_lambda)[1] <- "intercept"

  return(one_vector_lambda)
}

#' Compute the penalty matrix for a functional basis
#'
#' @description
#' Calculates the penalty matrix associated with a specific Linear Differential
#' Operator (LDO). By default, it computes the inner product of the second
#' derivatives of the basis functions, which is typically used to enforce
#' smoothness in penalized functional regression.
#'
#' @param basis_obj A functional data object (`fd`), a basis object (`basis.fd`),
#' or a list of `fd` objects representing individual basis functions.
#' @param LDO A numeric integer defining the Linear Differential Operator (LDO).
#' Defaults to 2 for the second derivative penalty.
#'
#' @return A square symmetric numeric matrix representing the penalty matrix.
#' @export
#'
#' @author Francois Bassac
#' @importFrom fda eval.penalty inprod int2Lfd
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_basis' is a basisfd object
#' pen_mat <- calc_penalty_matrix(my_basis, LDO = 2)
#' }
calc_penalty_matrix <- function(basis_obj, LDO = 2) {

  # 1. Check that LDO is numeric
  if (!is.numeric(LDO)) {
    stop("LDO must be a numeric integer representing the derivative order.")
  }
  LDO <- as.integer(LDO)

  # 2. Case: basis_obj is a list of 'fd' objects
  if (is.list(basis_obj) && all(vapply(basis_obj, inherits, logical(1), "fd"))) {
    n_basis <- length(basis_obj)
    penalty_matrix <- matrix(0, nrow = n_basis, ncol = n_basis)


    # Convert the integer LDO into an actual Lfd object for inprod
    Lfd_obj <- fda::int2Lfd(LDO)

    # Compute the inner product of the derivatives for each pair of functions
    for (i in seq_len(n_basis)) {
      for (j in i:n_basis) {
        # inprod computes the integral of the product of the LDO derivatives
        val <- fda::inprod(basis_obj[[i]], basis_obj[[j]],
                           Lfdobj1 = Lfd_obj, Lfdobj2 = Lfd_obj)
        penalty_matrix[i, j] <- val
        penalty_matrix[j, i] <- val # The penalty matrix is strictly symmetric
      }
    }
    return(penalty_matrix)
  }

  # 3. Case: basis_obj is a single 'fd' object
  if (inherits(basis_obj, "fd")) {
    basis_obj <- basis_obj$basis
  }

  # 4. Case: basis_obj is a 'basis.fd' object
  if (inherits(basis_obj, "basisfd")) {
    penalty_matrix <- fda::eval.penalty(basisobj = basis_obj, Lfdobj = LDO)
    return(penalty_matrix)
  }

  # Fallback if the user passes an unsupported object type
  stop("basis_obj must be a 'basisfd' object, an 'fd' object, or a list of 'fd' objects.")
}


#' Augment a penalty matrix with zero padding for the intercept
#'
#' @description
#' Expands a penalty matrix by adding a first row and a first column of zeros.
#' This is required in penalized regression models to ensure that the intercept
#' term is not penalized, while keeping the matrix dimensions aligned with the
#' augmented design matrix.
#'
#' @param pen_matrix A square numeric matrix representing the penalty matrix.
#'
#' @return A square matrix with dimensions `nrow(pen_matrix) + 1` by
#' `ncol(pen_matrix) + 1`, where the first row and first column are exactly zero.
#' @export
#'
#' @author Francois Bassac
#'
#' @examples
#' pen_mat <- matrix(c(2, -1, -1, 2), nrow = 2)
#' augment_penalty_matrix(pen_mat)
augment_penalty_matrix <- function(pen_matrix) {

  # Ensure the input is treated as a matrix
  pen_matrix <- as.matrix(pen_matrix)

  # Add a column of 0s on the left
  mat_with_col <- cbind(0, pen_matrix)

  # Add a row of 0s on the top
  aug_matrix <- rbind(0, mat_with_col)

  # Name the intercept row and column for consistency and readability
  #rownames(aug_matrix)[1] <- "intercept"
  #colnames(aug_matrix)[1] <- "intercept"

  return(aug_matrix)
}

#' Fit Univariate Penalized Functional Regression
#'
#' @description
#' Solves the penalized functional regression for a single, specific penalty
#' parameter (lambda). It computes the Ridge-like estimator using base R
#' matrix operations for maximum speed.
#'
#' @param X The design matrix (usually basis evaluations, augmented with an intercept).
#' @param Y The numeric response vector.
#' @param R The penalty matrix (augmented with zero-padding for the intercept).
#' @param lambda A single numeric value for the penalty parameter.
#'
#' @return A list containing the estimated coefficients (`beta_hat`),
#' the fitted values (`Y_hat`), and the residuals (`residuals`).
#' @export
#'
#' @author Francois Bassac
fit_pfr_uni <- function(X, Y, R, lambda) {

  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)

  # A = (X'X + lambda * R)
  A <- XtX + lambda * R

  # Compute coefficients: beta = A^(-1) X'Y
  A_inv <- solve(A)
  beta_hat <- A_inv %*% XtY

  # Predictions and residuals
  Y_hat <- X %*% beta_hat
  residuals <- Y - Y_hat

  return(list(
    beta_hat = as.vector(beta_hat),
    Y_hat = as.vector(Y_hat),
    residuals = as.vector(residuals)
  ))
}

#' Leave-One-Out Cross-Validation for Univariate PFR
#'
#' @description
#' Finds the optimal penalty parameter (lambda) over a specified grid using
#' the exact Leave-One-Out Cross-Validation (LOOCV). It uses the hat matrix
#' diagonal trick to compute the LOOCV error without refitting the model N times.
#'
#' @param X The design matrix (augmented with an intercept).
#' @param Y The numeric response vector.
#' @param R The penalty matrix (augmented).
#' @param lambda_grid A numeric vector of penalty values to test.
#'
#' @return An object of class `cv_pfr_uni` containing the optimal lambda,
#' the minimum RMSEP, and a data frame of all grid results.
#' @export
#'
#' @author Francois Bassac
cv_pfr_uni <- function(X, Y, R, lambda_grid) {

  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  n <- nrow(X)

  rmsep_vec <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    A <- XtX + lam * R

    # solve() can fail if matrix is singular (e.g. lambda = 0 with many splines).
    # tryCatch prevents the whole loop from crashing.
    A_inv <- tryCatch(solve(A), error = function(e) NULL)

    if (is.null(A_inv)) {
      rmsep_vec[i] <- Inf
      next
    }

    beta_hat <- A_inv %*% XtY
    Y_hat <- X %*% beta_hat
    res <- Y - Y_hat

    # Exact LOOCV trick using the diagonal of the Hat matrix: diag(X * A^-1 * X')
    # rowSums(...) is a highly optimized way to extract the diagonal
    h_diag <- rowSums((X %*% A_inv) * X)

    # LOOCV residuals: e_i / (1 - h_ii)
    res_loo <- res / (1 - h_diag)

    # Root Mean Squared Error of Prediction (RMSEP)
    rmsep_vec[i] <- sqrt(mean(res_loo^2))
  }

  best_idx <- which.min(rmsep_vec)

  res_obj <- list(
    optimal_lambda = lambda_grid[best_idx],
    min_rmsep = rmsep_vec[best_idx],
    cv_results = data.frame(lambda = lambda_grid, rmsep = rmsep_vec)
  )

  class(res_obj) <- "cv_pfr_uni"
  return(res_obj)
}

#' Plot method for cv_pfr_uni objects
#'
#' @description
#' Visualizes the LOOCV RMSEP curve as a function of the penalty parameter lambda.
#'
#' @param x An object of class `cv_pfr_uni`.
#' @param ... Additional graphical parameters.
#'
#' @export
#' @author Francois Bassac
plot.cv_pfr_uni <- function(x, ...) {

  # Using a logarithmic scale for the X axis is standard for lambda grids
  plot(x$cv_results$lambda, x$cv_results$rmsep,
       type = "b", log = "x",
       col = "steelblue", pch = 16, lwd = 2,
       xlab = expression(lambda),
       ylab = "LOOCV RMSEP",
       main = "LOOCV Error vs Penalty Parameter")

  # Add a vertical red dashed line at the optimal lambda
  abline(v = x$optimal_lambda, col = "firebrick", lty = 2, lwd = 2)

  # Add a legend
  legend("topright", legend = paste("Optimal lambda:",
                                    signif(x$optimal_lambda, 4)),
         text.col = "firebrick", bty = "n")
}


#' Fit a Penalized Functional Regression Model
#'
#' @description
#' Fits a univariate Penalized Functional Regression (PFR) model. It automatically
#' handles data formatting, penalty matrix computation, and hyperparameter
#' tuning via exact Leave-One-Out Cross-Validation (LOOCV).
#'
#' @param Y A numeric vector representing the scalar response variable.
#' @param X_func The functional predictor. Can be a raw data matrix or an `fd` object.
#' @param basis_obj A `basis.fd` object defining the functional basis.
#' @param type A character string specifying the type of functional data:
#' `"cfd"` (Continuous Functional Data) or `"sfd"` (Step/Sparse Functional Data).
#' @param lambda_grid A numeric vector of penalty parameters to evaluate during CV.
#' @param LDO An integer defining the Linear Differential Operator for the penalty. Defaults to 2.
#'
#' @return An object of class `pfr` containing the fitted model, optimal lambda,
#' cross-validation results, and estimated functional coefficients.
#' @export
#'
#' @author Francois Bassac
#' @importFrom fda inprod Data2fd
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' model <- pfr(Y = my_target, X_func = my_curves, basis_obj = my_basis, type = "cfd")
#' }
pfr <- function(Y, X_func, basis_obj, type = "cfd",
                lambda_grid = 10^seq(-5, 5, length.out = 30), LDO = 2) {


  # 1. INPUT VALIDATION

  if (!is.numeric(Y)) stop("Y must be a numeric vector.")
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0)
    stop("lambda_grid must be a numeric vector.")
  if (!type %in% c("cfd", "sfd")) stop("type must be either 'cfd' or 'sfd'.")

  # Ensure Y is a standard vector and drop empty dimensions
  Y <- as.vector(Y)
  n_obs <- length(Y)

  # 2. DATA PROCESSING & DISPATCH ('cfd' vs 'sfd')
  # Objective: Extract the coefficient matrix (C) of the functional data

  if (inherits(X_func, "fd")) {
    # If the user already provides an fd object, extract coefficients directly
    coef_matrix <- t(X_func$coefs)

  } else {
    # If the user provides a raw matrix, we need to convert it based on 'type'
    X_func <- as.matrix(X_func)
    if (nrow(X_func) != n_obs) {
      stop("The number of rows in X_func must match the length of Y.")
    }

    if (type == "cfd") {
      # Standard continuous functional data smoothing
      # Assuming points are evenly spaced in the basis range
      argvals <- seq(basis_obj$rangeval[1],
                     basis_obj$rangeval[2],
                     length.out = ncol(X_func))
      fd_obj <- fda::Data2fd(argvals = argvals,
                             y = t(X_func),
                             basisobj = basis_obj)
      coef_matrix <- t(fd_obj$coefs)

    } else if (type == "sfd") {
      # Specific logic for SFD (Step Functional Data)
      # [!] Replace this block with your exact SmoothPLS logic for SFD
      # For now, we assume it's handled similarly or requires a specific projection
      warning("SFD processing is currently using the default continuous projection.")
      argvals <- seq(basis_obj$rangeval[1],
                     basis_obj$rangeval[2],
                     length.out = ncol(X_func))
      fd_obj <- fda::Data2fd(argvals = argvals,
                             y = t(X_func), basisobj = basis_obj)
      coef_matrix <- t(fd_obj$coefs)
    }
  }


  # 3. MATHEMATICAL ENGINE (Design Matrix & Penalty)

  # Calculate the inner product matrix of the basis functions (J)
  # J_ij = integral( phi_i(t) * phi_j(t) dt )
  J_matrix <- fda::inprod(basis_obj, basis_obj)

  # The true design matrix X for functional regression is C %*% J
  X_design <- coef_matrix %*% J_matrix

  # Augment X with an intercept
  X_aug <- lambda_augmented_uni(X_design)

  # Calculate and augment the penalty matrix (R)
  R_matrix <- calc_penalty_matrix(basis_obj, LDO = LDO)
  R_aug <- augment_penalty_matrix(R_matrix)


  # 4. MODEL FITTING & CROSS-VALIDATION

  # Run the highly optimized LOOCV to find the best lambda
  cv_results <- cv_pfr_uni(X = X_aug,
                           Y = Y,
                           R = R_aug,
                           lambda_grid = lambda_grid)
  best_lambda <- cv_results$optimal_lambda

  # Fit the final model using the optimal lambda
  final_fit <- fit_pfr_uni(X = X_aug, Y = Y, R = R_aug, lambda = best_lambda)

  # Separate intercept from functional coefficients
  intercept <- final_fit$beta_hat[1]
  beta_coefs <- final_fit$beta_hat[-1]

  # Reconstruct the functional beta(t) as an 'fd' object
  beta_fd <- fda::fd(as.matrix(beta_coefs), basis_obj)

  # 5. RETURN S3 OBJECT

  res <- list(
    call = match.call(),
    type = type,
    intercept = intercept,
    beta_fd = beta_fd,           # The functional coefficient beta(t)
    beta_coefs = beta_coefs,     # Raw coefficients
    fitted.values = final_fit$Y_hat,
    residuals = final_fit$residuals,
    optimal_lambda = best_lambda,
    min_rmsep = cv_results$min_rmsep,
    cv_object = cv_results       # Keep this to allow plot(model$cv_object)
  )

  class(res) <- "pfr"
  return(res)
}


# Do the univariate pfr on CFD and SFD
# OK for UNI 1-CFD case
mpfr <- function(df_list, Y,
                 basis_obj, regul_time_obj = NULL,
                 curve_type_obj, lambda_grid = 10^seq(-5, 10, length.out = 15),
                 id_col_obj = 'id', time_col_obj = 'time', int_mode = 1,
                 print_steps = FALSE,
                 plot_rmsep = TRUE,
                 plot_reg_curves = FALSE,
                 parallel = TRUE){

  # Step 1 : assertion
  if(print_steps){
    cat("=> Input format assertions.\n")
  }
  assert_obj = assert_multivariate_smoothPLS_inputs(
    df_list = df_list,
    Y = Y,
    basis_obj = basis_obj,
    regul_time_obj = regul_time_obj,
    curve_type_obj = curve_type_obj,
    orth_obj = FALSE,
    id_col_obj = id_col_obj,
    time_col_obj = time_col_obj)

  N_curves = assert_obj$N_curves
  basis_list = assert_obj$basis_list
  regul_time_list = assert_obj$regul_time_list
  curve_type_list = assert_obj$curve_type_list
  id_col_list = assert_obj$id_col_list
  time_col_list = assert_obj$time_col_list
  orth_list = assert_obj$orth_list

  if(print_steps){
    cat("=> Input format assertions OK.\n")
  }

  # Step 2 create orth_basis_list
  if(print_steps){
    cat("=> Create list of basis functions. \n")
  }

  basis_list_obj = orthonormalize_basis_list(basis_list = basis_list,
                                              orth_list = FALSE)


  # Step 3 build df_processed_list and curves_names_list
  if(print_steps){
    cat("=> Data objects formatting.\n")
  }

  if(N_curves == 1 && mode(df_list[[1]]) != 'list' && ncol(df_list) == 3){
    df_list = list(df_list)
  }

  new_list_obj = build_new_data_list(df_list = df_list,
                                     N_curves = N_curves,
                                     orth_basis_list = basis_list_obj,
                                     basis_list = basis_list,
                                     curve_type_list = curve_type_list,
                                     id_col_list = id_col_list,
                                     time_col_list = time_col_list,
                                     regul_time_list = regul_time_list)

  df_processed_list = new_list_obj$df_processed_list
  curves_names_list = new_list_obj$curves_names_list
  new_curves_type_list = new_list_obj$new_curves_type_list
  new_basis_list = new_list_obj$new_basis_list
  #new_orth_basis_list = new_list_obj$new_orth_basis_list
  new_id_col_list = new_list_obj$new_id_col_list
  new_time_col_list = new_list_obj$new_time_col_list
  new_regul_time_list = new_list_obj$new_regul_time_list

  # Step 4 Build Lambda matrix
  if(print_steps){
    cat("=> Evaluate Lambda matrix.\n")
  }

  for(i in 1:length(df_processed_list)){

    if(print_steps){
      cat(paste0("==> Lambda for : ", curves_names_list[i], ".\n"))
    }


    lambda = evaluate_lambda(df = df_processed_list[[i]],
                             basis = new_basis_list[[i]],
                             curve_type = new_curves_type_list[[i]],
                             int_mode = int_mode,
                             id_col = new_id_col_list[[i]],
                             time_col = new_time_col_list[[i]],
                             regul_time = new_regul_time_list[[i]],
                             parallel = parallel)
    #dim(lambda)
    if(i != 1){
      Lambda = cbind(Lambda, lambda)
    }else{
      Lambda = lambda
    }
  }


  # Step 5 Build augmented objects
  X_aug = lambda_augmented_uni(Lambda)
  R = calc_penalty_matrix(new_basis_list[[1]])
  R_aug = augment_penalty_matrix(R)


  # Step 6 LOOCV
  cv_res <- cv_pfr_uni(X_aug, Y, R_aug, lambda_grid = lambda_grid)

  if(plot_rmsep){plot(cv_res)}


  # Step 7 model fit with best parameters
  modele_final <- fit_pfr_uni(X_aug, Y, R_aug,
                              lambda = cv_res$optimal_lambda)

  # Improve for multivariate!
  if(plot_reg_curves){
    plot(fda::fd(coef = modele_final$beta_hat[-1], basisobj = basis_obj))
  }

  # Step 8 build functional coefficients
  delta_list = list(
    modele_final$beta_hat[1],
    fda::fd(coef = modele_final$beta_hat[-1], basisobj = basis_obj)
  )
  names(delta_list) = c("Intercept", curves_names_list)


  mpfr_obj = list(cv_res, modele_final, delta_list)
  names(mpfr_obj) = c("cv_res", "modele_final", "reg_obj")

  return(mpfr_obj)

}

# WIP
build_reg_curve_mpfr_old <- function(mpfr_model, curves_names_list,
                                 print_steps = TRUE){


  N_curves_processed = length(curves_names_list)

  delta_list = list()
  for(i in 1:N_curves_processed){
    if(print_steps){
      cat(paste0("==> Build regression curve for : ",
                 curves_names_list[[i]], "\n"))
    }
    d_i = evaluate_reg_curve_PFR_uni(mpfr_model = mpfr_model,
                                      curve_name = curves_names_list[[i]])
    delta_list = append(delta_list, list(d_i))
  }

  delta_0 = mpfr_model$beta_hat[1]


  delta_spls = list(delta_0)
  for(i in 1:length(delta_list)){
    delta_spls = append(delta_spls, list(delta_list[[i]]))
  }
  names(delta_spls) = c("Intercept", curves_names_list)
  return(delta_spls)
}

# WIP
evaluate_reg_curve_PFR_uni_old <- function(mpfr_model, curve_name){
  # nb_comp=length(v_i_list)

  delta = fda::fd(coef=rep(0,v_i_list[[1]]$basis$nbasis),
                  basisobj = v_i_list[[1]]$basis)

  if(is.null(nb_comp)){
    # take ALL the components >> all the v_i(t)
    nb_stop = length(v_i_list)
  }else if(nb_comp > length(v_i_list)){
    stop("evaluate_reg_curve_PFR_uni() :
           nb_comp superior than the number of PLS steps!")
  }else{
    nb_stop = nb_comp
  }
  for(i in 1:nb_stop){
    delta = delta + plsr_model$Yloadings[i] * v_i_list[[i]]
  }

  return(delta)
}



#' Evaluate the functional regression coefficient for univariate PFR
#'
#' @description
#' Reconstructs the continuous functional coefficient curve beta(t) using
#' the estimated basis coefficients from the penalized functional regression model.
#'
#' @param fit_obj The fitted model object from `fit_pfr_uni` (must contain `beta_hat`).
#' @param basis_obj The `basisfd` object used to expand the functional predictor.
#'
#' @return An `fd` object representing the smooth functional coefficient beta(t).
#' @export
#'
#' @author Francois Bassac
evaluate_reg_curve_PFR_uni <- function(fit_obj, basis_obj) {

  # 1. Vérification de sécurité
  if (!inherits(basis_obj, "basisfd")) {
    stop("basis_obj must be a valid fda 'basisfd' object.")
  }

  # 2. Extraction des coefficients
  # Le premier élément est l'intercept, on prend tout le reste
  b_hat <- fit_obj$beta_hat[-1]

  # 3. Reconstruction de la courbe avec fda::fd
  # fda::fd s'attend à recevoir une matrice colonne pour les coefficients
  beta_fd <- fda::fd(coef = as.matrix(b_hat), basisobj = basis_obj)

  return(beta_fd)
}

build_reg_curve_mpfr <- function(modele_final, curves_names_list,
                                 basis_obj, print_steps = TRUE) {

  if (print_steps) {
    cat(paste0("==> Build regression curve for : ", curves_names_list[[1]], "\n"))
  }

  # 1. Évaluation de la courbe fonctionnelle
  delta_curve <- evaluate_reg_curve_PFR_uni(fit_obj = modele_final,
                                            basis_obj = basis_obj)

  # 2. Récupération de l'intercept
  delta_0 <- modele_final$beta_hat[1]

  # 3. Formatage de la liste de sortie
  delta_list <- list(delta_0, delta_curve)
  names(delta_list) <- c("Intercept", curves_names_list[[1]])

  return(delta_list)
}

# Multivariate

#### ETAPE 1 : LES USINES A BLOCS (Block Factories) ####

#' Build block for Non-Functional Data (NFD / Scalars)
#'
#' @param X_nfd A numeric matrix or data.frame of scalar predictors (n x Q)
#' @return A list with the design block Z and the penalty base R (Identity)
build_block_nfd <- function(X_nfd) {
  Z_block <- as.matrix(X_nfd)

  # Pour les scalaires, on ne peut pas dériver.
  # La pénalité de base est donc une matrice Identité (Ridge classique).
  # On pourra la multiplier par 0 plus tard si on ne veut pas pénaliser.
  R_block <- diag(ncol(Z_block))

  return(list(
    type = "nfd",
    Z = Z_block,
    R = R_block,
    n_coefs = ncol(Z_block) # Q
  ))
}

#' Build block for Scalar Functional Data (SFD)
#'
#' @param coef_matrix Matrix of coefficients (n x d) of the curves on their basis
#' @param basis_obj The fda basis object
#' @param LDO Linear Differential Operator for the penalty (default 2)
#' @return A list with the design block Z and the penalty base R
build_block_sfd <- function(coef_matrix, basis_obj, LDO = 2) {
  # Interaction temporelle : integral(phi * psi)
  J_matrix <- fda::inprod(basis_obj, basis_obj)

  Z_block <- coef_matrix %*% J_matrix
  R_block <- calc_penalty_matrix(basis_obj, LDO = LDO)

  return(list(
    type = "sfd",
    Z = Z_block,
    R = R_block,
    n_coefs = ncol(Z_block) # q
  ))
}

#' Build block for Categorical Functional Data (CFD)
#'
#' @param df_cfd The dataframe containing the CFD for all individuals
#' @param basis_obj The fda basis object
#' @param reference_state Character. The state to drop (K-1). If NULL, keeps all K states.
#' @param ... Additional arguments to pass to your evaluate_lambda function
#' @return A list with the concatenated design block Z and the base penalty R
build_block_cfd <- function(df_cfd, basis_obj, reference_state = NULL, LDO = 2, ...) {

  # Identifier tous les états disponibles
  # (On suppose ici que tu as une colonne d'état, adapte selon le nom de ta variable)
  all_states <- unique(df_cfd$state)

  # Gestion de l'option K ou K-1
  states_to_keep <- all_states
  if (!is.null(reference_state)) {
    states_to_keep <- setdiff(all_states, reference_state)
  }

  K_prime <- length(states_to_keep)

  # Boucle pour calculer l'aire active (Lambda) de chaque état conservé
  Z_list <- list()
  for (k in 1:K_prime) {
    current_state <- states_to_keep[k]

    # On filtre le dataframe pour ne garder que l'état en cours
    df_state <- df_cfd[df_cfd$state == current_state, ]

    # Ton incroyable fonction de calcul d'aire active
    Lambda_k <- evaluate_lambda(df = df_state, basis = basis_obj, ...)

    Z_list[[k]] <- Lambda_k
  }

  # Concaténation horizontale de tous les blocs Lambda (Z = [L1, L2, ..., LK'])
  Z_block <- do.call(cbind, Z_list)

  # La matrice de pénalité de base pour UNE courbe de cet état
  R_base <- calc_penalty_matrix(basis_obj, LDO = LDO)

  return(list(
    type = "cfd",
    Z = Z_block,
    R = R_base,       # On renvoie la base, l'Assembleur la dupliquera K_prime fois !
    n_coefs = ncol(Z_block), # K_prime * q
    K_prime = K_prime,
    states = states_to_keep
  ))
}


#### ETAPE 2 : L'ASSEMBLEUR GLOBAL (The Matrix Assembler) ####


#' Build the Global Design Matrix (D or Z)
#'
#' @description
#' Concatenates the intercept and all individual predictor blocks horizontally.
#'
#' @param block_list A list of blocks generated by the build_block_* functions.
#' @param n_obs The number of observations (rows).
#'
#' @return A matrix of dimension n x (1 + P_total)
build_global_design <- function(block_list, n_obs) {

  # 1. Initialiser avec la colonne de l'intercept
  D_global <- matrix(1, nrow = n_obs, ncol = 1)
  colnames(D_global) <- "Intercept"

  # 2. Concaténer horizontalement tous les blocs Z
  for (i in seq_along(block_list)) {
    D_global <- cbind(D_global, block_list[[i]]$Z)
  }

  return(D_global)
}


#' Helper function: Create a block-diagonal matrix from a list of matrices
#' (Pure base R for speed and zero dependencies)
bdiag_base <- function(mat_list) {
  # Calcul des dimensions totales
  dims <- sapply(mat_list, dim)
  total_rows <- sum(dims[1, ])
  total_cols <- sum(dims[2, ])

  # Initialisation d'une matrice vide
  res <- matrix(0, nrow = total_rows, ncol = total_cols)

  # Remplissage de la diagonale
  current_row <- 1
  current_col <- 1
  for (mat in mat_list) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    res[current_row:(current_row + nr - 1), current_col:(current_col + nc - 1)] <- mat
    current_row <- current_row + nr
    current_col <- current_col + nc
  }

  return(res)
}


#' Build the Augmented Global Penalty Matrix (R_0_lambda)
#'
#' @description
#' Constructs the block-diagonal penalty matrix dynamically. Supports both
#' global lambda per predictor and specific lambdas per CFD state.
#'
#' @param block_list A list of blocks generated by the build_block_* functions.
#' @param lambda_list A list of numeric vectors containing the lambdas for each block.
#'
#' @return A symmetric square matrix representing R_0_lambda
build_global_penalty <- function(block_list, lambda_list) {

  if (length(block_list) != length(lambda_list)) {
    stop("lambda_list must have exactly one element per predictor block.")
  }

  penalties_to_diag <- list()

  for (i in seq_along(block_list)) {
    block <- block_list[[i]]
    lams <- lambda_list[[i]] # Ceci est un vecteur de taille 1 ou K'

    if (block$type == "cfd") {
      # Cas 1 : Un seul lambda global pour tout le CFD
      if (length(lams) == 1) {
        R_penalized <- lams * block$R
        cfd_blocks <- replicate(block$K_prime, R_penalized, simplify = FALSE)
        penalties_to_diag[[i]] <- bdiag_base(cfd_blocks)

        # Cas 2 : Un lambda spécifique par état
      } else if (length(lams) == block$K_prime) {
        cfd_blocks <- list()
        for (k in 1:block$K_prime) {
          cfd_blocks[[k]] <- lams[k] * block$R
        }
        penalties_to_diag[[i]] <- bdiag_base(cfd_blocks)

        # Erreur de saisie utilisateur
      } else {
        stop(paste("Block", i, "is a CFD with", block$K_prime,
                   "states. You must provide either 1 or", block$K_prime, "lambdas."))
      }

    } else {
      # Pour NFD et SFD, il ne faut strictement qu'un seul lambda
      if (length(lams) != 1) {
        stop(paste("Block", i, "is of type", block$type, "and requires exactly 1 lambda."))
      }
      penalties_to_diag[[i]] <- lams * block$R
    }
  }

  # 2. Assembler la méga-matrice sur la diagonale
  R_lambda_global <- bdiag_base(penalties_to_diag)

  # 3. Augmenter la matrice (R_0_lambda) pour protéger l'intercept
  R_0_lambda <- augment_penalty_matrix(R_lambda_global)

  return(R_0_lambda)
}
