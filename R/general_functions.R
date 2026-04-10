#### Basis ####

#' create_bspline_basis
#'
#' @param start start time
#' @param end end time
#' @param nbasis number of basis functions, default 10
#' @param norder order of the basis function, default cubic splines 4
#'
#' @returns a basis fd object
#' @export
#'
#' @importFrom fda create.bspline.basis
#'
#' @examples
#' b0 = create_bspline_basis(0, 10, 10, 4)
#' plot(b0)
#'
#' b1 = create_bspline_basis(0, 10, 10, 2)
#' plot(b1)
#'
#' b2 = create_bspline_basis(0, 10, 10, 1)
#' plot(b1)
#'
#'
#' @author Francois Bassac
create_bspline_basis <- function(start, end, nbasis=10, norder=4){
  #This function creates a bspline basis.
  b = fda::create.bspline.basis(c(start, end),
                                nbasis = nbasis,
                                norder = norder)
  return(b)
}


#' from_basis_to_fdlist
#'
#' This function transform if necessary the input basis into a list of fd
#' functions.
#' If basis is a basis object from fda, the output fd_list is the list of the
#' different basis functions as fd functions.
#' If basis is already a list of fd functions, nothing changes.
#'
#' @param basis a basis fd object or a list of fd functions.
#'
#' @returns a list of fd functions
#' @export
#'
#' @importFrom fda fd
#'
#' @examples
#' basis = create_bspline_basis(start = 0, end = 10, nbasis = 5, norder = 4)
#' plot(basis)
#'
#' basis_list = from_basis_to_fdlist(basis)
#'
#' plot(basis_list[[1]], col = 1)
#' for(i in 2:length(basis_list)){
#'   plot(basis_list[[i]], add=TRUE, col = i)
#' }
#'
#' basis_list_2 = from_basis_to_fdlist(basis_list)
#'
#' @author Francois Bassac
from_basis_to_fdlist <- function(basis){

  # If necessary, convert into a fd list
  if (inherits(basis, "basisfd")) {
    nbasis <- basis$nbasis
    rangeval <- basis$rangeval
    coefs <- diag(nbasis)
    fd_list <- list()
    for (i in 1:nbasis) {
      fd_list[[i]] <- fda::fd(coefs[, i], basis)
    }
  } else if (is.list(basis) && all(sapply(basis, function(x) inherits(x, "fd")))) {
    fd_list <- basis
  } else {
    stop("from_basis_to_fdlist() : basis have to be a basisobj from fda or a list of fd objects")
  }
  return(fd_list)
}

#' obj_list_creation
#'
#' This functions creates a list of basis containing the
#' value of the number of states or curves sharing the same basis.
#'
#' @param N_rep a value of the number of states sharing the same basis
#' @param obj a object
#'
#' @returns a list of N_rep appended obj
#' @export
#'
#' @examples
#' N_rep = 4
#' start = 1
#' end = 51
#' nbasis = 13
#' norder = 3
#'
#' basis = create_bspline_basis(start, end, nbasis, norder)
#' basis_list = obj_list_creation(N_rep, basis)
#' obj_list_creation(N_rep, 0:100)
#'
#' @author Francois Bassac
obj_list_creation <- function(N_rep, obj){
  # This functions creates a list of basis.
  for(i in 1:N_rep){
    if(i==1){
      obj_list = list(obj)
    }else{
      obj_list <- append(obj_list, list(obj))
    }
  }
  return(obj_list)
}

#' evaluate_metric
#'
#' This function evaluates the metric of a certain basis.
#' The metric is the inprod of the basis functions.
#'
#' @param basis basis to evaluate the metric
#'
#' @returns a matrix of dimension nbasis X nbasis
#' @export
#'
#' @importFrom fda norder.bspline inprod
#'
#' @examples
#' basis = create_bspline_basis(start=0, end=10, nbasis=10, norder=4)
#' metric = evaluate_metric(basis)
#'
#' basis1 = create_bspline_basis(start=0, end=20, nbasis=10, norder=1)
#' metric1 = evaluate_metric(basis1)
#'
#' @author Francois Bassac
evaluate_metric <- function(basis){
  # This function evaluate the metric of a basis.

  if(basis$type == 'bspline' && fda::norder.bspline(basis) == 1){
    end = basis$rangeval[2]
    metric = diag(rep((end/basis$nbasis), basis$nbasis))
  }else{
    metric = fda::inprod(basis, basis)
  }
  return(metric)
}

#' assemble_basis_metric
#'
#' This function assemble the metrics of all the basis of the basis list.
#' This function only assemble the needed basis, especially if length(curve_to_keep) != N_states
#'
#' @param basis_list a list of basis fd object
#' @param curves_to_keep a list of the states curves to keep
#'
#' @returns a matrix of the metric to consider
#' @export
#'
#' @importFrom fda norder.bspline inprod
#'
#' @examples
#' basis1 = fda::create.bspline.basis(c(0,100), nbasis=10, norder=4)
#' basis2 = fda::create.bspline.basis(c(0,100), nbasis=15, norder=1)
#' basis3 = fda::create.fourier.basis(c(0,100), nbasis=7)
#' assemble_basis_metric(list(basis1, basis2, basis3), list(1,2,4))
#' assemble_basis_metric(list(basis1, basis2, basis3), list(1,2))
#'
#' @author Francois Bassac
assemble_basis_metric <- function(basis_list, curves_to_keep=NULL){
  # This function assemble the metrics of all the basis of the basis list.
  # This function only assemble the needed basis, especially if length(curve_to_keep) != N_states

  #basis1 = fda::create.bspline.basis(c(0,100), nbasis=10, norder=4)
  #basis2 = fda::create.bspline.basis(c(0,100), nbasis=15, norder=1)
  #basis3 = fda::create.fourier.basis(c(0,100), nbasis=7)
  # assemble_basis_metric(list(basis1, basis2, basis3), list(1,2,4))
  # assemble_basis_metric(list(basis1, basis2, basis3), list(1,2))

  if(is.null(curves_to_keep)){
    curves_to_keep = c(1:length(basis_list))
  }

  for(i in c(1:length(curves_to_keep))){
    if(i == 1){
      # i = 1
      metric = evaluate_metric(basis = basis_list[[i]])
    }else{
      metric0 = evaluate_metric(basis = basis_list[[i]])
      metric <- block_diag(metric, metric0)
    }
  }
  return(metric)
}

#### Orthogonalization ####

#' is_orthonormal
#'
#' Check if a basis function is orthonormal
#'
#' @param basis A basis object from fda package or a list of fd functions.
#' @param tol a float, tolerance parameter, default 1e-10)
#'
#' @return A boolean, TRUE if orthonormal, FALSE if not
#'
#' @importFrom fda inprod fd
#' @export
#'
#' @author Francois Bassac
is_orthonormal <- function(basis, tol = 1e-10) {

  if (inherits(basis, "basisfd")) {
    gram_matrix <- fda::inprod(basis, basis)
    nbasis <- basis$nbasis
  }
  # if list(fd.objects)
  else if (is.list(basis) && all(sapply(basis, inherits, "fd"))) {
    nbasis <- length(basis)

    coefs <- sapply(basis, function(x) x$coefs)
    base_obj <- basis[[1]]$basis
    # Gram = t(C) %*% Inprod(Phi, Phi) %*% C
    gram_matrix <- t(coefs) %*% fda::inprod(base_obj, base_obj) %*% coefs
  } else {
    stop("is_orthonormal() : needs 'basisfd' or 'fd' object list.")
  }

  # Vérification
  identity_matrix <- diag(nbasis)
  return(all(abs(gram_matrix - identity_matrix) < tol))
}

#' is_orthogonal
#'
#' Check if a basis function is orthogonal
#'
#' @param basis A basis object from fda package or a list of fd functions.
#' @param tol a float, tolerance parameter, default 1e-10)
#'
#' @return A boolean, TRUE if orthogonal, FALSE if not
#'
#' @importFrom fda inprod fd
#' @export
#'
#' @author Francois Bassac
is_orthogonal <- function(basis, tol = 1e-10) {

  if (inherits(basis, "basisfd")) {
    gram_matrix <- fda::inprod(basis, basis)
  } else if (is.list(basis) && all(sapply(basis, inherits, "fd"))) {
    coefs <- sapply(basis, function(x) x$coefs)
    base_obj <- basis[[1]]$basis
    gram_matrix <- t(coefs) %*% fda::inprod(base_obj, base_obj) %*% coefs
  } else {
    stop("is_orthogonal() : needs 'basisfd' or 'fd' object list")
  }

  diag(gram_matrix) <- 0
  return(all(abs(gram_matrix) < tol))
}

#' gram_schmidt_orthonormalize
#'
#' Orthonormalizer a basis functions with Gram-Schmidt algorithm.
#'
#' @param basis A basis object from fda package or a list of fd functions.
#' @param output_type A character to choose the output format. "fdlist"
#' (default) or "funlist" for R functions.
#' @param tol a float, tolerance parameter, default 1e-12
#'
#' @return A list of orthonormalized functions fd or func(t)
#' @export
#'
#' @importFrom fda inprod fd eval.fd
#'
#' @examples
#'
#' start = 0
#' end = 10
#' basis = create_bspline_basis(start, end, nbasis = 10, norder = 4)
#'
#' basis_orth = gram_schmidt_orthonormalize(basis, "fdlist")
#'
#' @author Francois Bassac
gram_schmidt_orthonormalize <- function(basis, output_type = "fdlist",
                                         tol = 1e-12) {

  # 1. Checks
  if (inherits(basis, "basisfd")) {
    nbasis <- basis$nbasis
  } else {
    stop("orthonormalize_basis() : 'basis' must be fda::basisfd object.")
  }

  # 2. Gram matrix
  G <- fda::inprod(basis, basis)

  # 3. Spectral decomposition
  eig <- eigen(G)

  # Eigenvalues cleaning (numerical stability)
  eig$values[eig$values < tol] <- 0
  if (any(eig$values == 0)) {
    # If zero, then real linear dependance
    nb_dep = sum(eig$values == 0)
    warning(paste(nb_dep, "base function(s) seam(s) : linear dependance."))
  }

  # 4. S Passage matrix C
  nonzero <- eig$values > 0
  S <- matrix(0, nbasis, nbasis)
  S[, nonzero] <- eig$vectors[, nonzero] %*% diag(1 / sqrt(eig$values[nonzero]))

  # 5. List of FD object creation
  orthonormal_list <- lapply(1:nbasis, function(i) {
    fda::fd(S[, i], basis)
  })

  # 6. Output format
  if (output_type == "fdlist") {
    return(orthonormal_list)
  } else if (output_type == "funlist") {
    return(lapply(orthonormal_list, function(fd_obj) {
      function(t) { as.numeric(fda::eval.fd(t, fd_obj)) }
    }))
  } else {
    stop("output_type have to be 'fdlist' or 'funlist'.")
  }
}

#' transition_matrix
#'
#' Build transition matrix between to basis.
#'
#' @param basis1 First basis, basis obj or list of fd functions
#' @param basis2 Second basis, basis obj or list of fd functions
#'
#' @return Transition matrix P such as basis1 = P * basis2
#' @export
#'
#' @importFrom fda inprod fd
#'
#' @author Francois Bassac
transition_matrix <- function(basis1, basis2) {

  # Auxiliar function to convert a basis object into a fd list.
  convert_to_fd_list <- function(basis) {
    if (inherits(basis, "basisfd")) {
      nbasis <- basis$nbasis
      coefs <- diag(nbasis)
      fd_list <- list()
      for (i in 1:nbasis) {
        fd_list[[i]] <- fda::fd(coefs[, i], basis)
      }
      return(fd_list)
    } else if (is.list(basis) && all(sapply(basis,
                                            function(x) inherits(x, "fd")))) {
      return(basis)
    } else {
      stop("transition_matrix() : basis have to be a basisobj from fda
         or a list of fd objects")
    }
  }

  # Convert the 2 basis
  fd_list1 <- convert_to_fd_list(basis1)
  fd_list2 <- convert_to_fd_list(basis2)

  n1 <- length(fd_list1)
  n2 <- length(fd_list2)

  # Build transition matrix
  # P[i,j] = <basis1[i], basis2[j]>
  P <- matrix(0, n1, n2)

  for (i in 1:n1) {
    for (j in 1:n2) {
      P[i, j] <- fda::inprod(fd_list1[[i]], fd_list2[[j]])
    }
  }

  return(P)
}

#' test_basis_properties
#'
#' General function to test basis functions
#'
#' @param basis Basis to test, basis object or list of fd functions
#' @param name Character, name of the basis
#' @param tol Float, precision, default 1e-10
#'
#' @return a boolean
#' @export
#' @examples
#' start = 0
#' end = 10
#' basis = create_bspline_basis(start, end, nbasis=10, norder=4)
#' test_basis_properties(basis, "cubic splines", tol = 1e-3)
#'
#' basis_f = fda::create.fourier.basis(rangeval=c(start, end), nbasis=5)
#' test_basis_properties(basis_f, "Fourier", tol = 1e-3)
#'
#' @author Francois Bassac
test_basis_properties <- function(basis, name = "Base", tol = 1e-10) {
  cat(paste("=== Tests for", name, "===\n"))
  cat(paste("Orthogonal:", is_orthogonal(basis, tol = tol), "\n"))
  cat(paste("Orthonormal:", is_orthonormal(basis, tol = tol), "\n"))
  cat("\n")
}



#### Integral evaluation ####

#' evaluate_id_func_integral
#'
#' @description
#' Evaluates the integral \eqn{\int_{\tau_i} f(t) dt} where \eqn{\tau_i} are the active
#' intervals of a categorical functional data (states 0 or 1).
#'
#' @param id_df Dataframe for a single individual with at least columns (id, time, state).
#' @param func The R function to integrate.
#' @param id_col Character, name of the id column, default 'id'.
#' @param time_col Character, name of the time column, default 'time'.
#' @param rel_tol Relative tolerance for stats::integrate, default 1e-8.
#' @param subdivisions Max number of subdivisions for integrate, default 100.
#' @param ... Additional arguments (ignored to prevent passing unused params to func).
#'
#' @return A dataframe with the id and the calculated integral value.
#' @export
evaluate_id_func_integral <- function(id_df, func,
                                      id_col = 'id', time_col = 'time',
                                      rel_tol = .Machine$double.eps^0.5,
                                      subdivisions = 1000L, ...) { # Augmenté à 1000

  state_col <- setdiff(names(id_df), c(id_col, time_col))
  id_df <- id_df[order(id_df[[time_col]]), ]

  integral_sum <- 0

  for(i in 1:(nrow(id_df) - 1)) {
    t1 <- id_df[[time_col]][i]
    t2 <- id_df[[time_col]][i + 1]
    current_state <- id_df[[state_col]][i]

    if(as.numeric(current_state) == 1 && t2 > t1) {
      # tryCatch
      segment_integral <- stats::integrate(f = func,
                                           lower = t1,
                                           upper = t2,
                                           rel.tol = rel_tol,
                                           subdivisions = subdivisions,
                                           stop.on.error = FALSE)

      # If fail we keep the best estimation
      integral_sum <- integral_sum + segment_integral$value

      if(segment_integral$message != "OK") {
        warning(paste("Integration issue at id", id_df[[id_col]][1], ":",
                      segment_integral$message))
      }
    }
  }
  return(data.frame(id = id_df[[id_col]][1], integral = integral_sum))
}


#' evaluate_id_func_integral
#'
#' This function evaluate the integral for a state (0, 1) functional data :
#' int( X(t) func(t) )dt.
#' This function works ONLY for a one state CFD!
#'
#' @param id_df a single id dataframe of at least named columns (id, time)
#' @param func the function to integrate
#' @param mode select the integration mode 1 for R function integrate,
#' 2 for pracma::trapz. default value : 1
#' @param id_col col_name of df for the id
#' @param time_col col_name of df for the time
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100

#'
#' @returns a dataframe with the id and the integral value.
#' @export
#'
#' @importFrom stats integrate
#' @importFrom pracma trapz
#'
#' @examples
#' id_df = data.frame(id=rep(1,5), time=seq(0, 40, 10), state=c(0, 1, 1, 0, 1))
#' evaluate_id_func_integral(id_df, function(t){t})
#'
#' @author Francois Bassac
evaluate_id_func_integral_deprecated <- function(id_df, func, mode = 1,
                                      id_col = 'id', time_col = 'time',
                                      nb_pt = 10, subdivisions = 100){
  # This function evaluate the integral of id_df$state * func on time interval.
  # This function works for only one id!

  # IMPROVEMENT add id_col, time_col, state_col?
  # OR find state_col by setdiff(names(id_df), c(id_col, time_col))
  state_col = setdiff(names(id_df), c(id_col, time_col))

  # Order
  id_df <- id_df[order(id_df[[time_col]]), ]

  # Init
  integral_sum <- 0

  # Loop on the id_df
  for(i in 1:(nrow(id_df) - 1)){
    t1 <- id_df[[time_col]][i]
    t2 <- id_df[[time_col]][i + 1]
    current_state <- id_df[[state_col]][i]

    # If current_state is 0, No integral evaluation
    if(as.numeric(current_state) == 0){
      next
    } else {
      # If current_state is 1, we evaluate the integral based on the trapezoidal
      # numerical method.

      if(mode==2){
        # V1
        # We will use at least 10 points
        n_points <- max(nb_pt, ceiling((t2 - t1) * nb_pt))
        n_points <- nb_pt

        points <- seq(t1, t2, length.out = n_points)
        values <- sapply(points, func)
        segment_integral <- pracma::trapz(points, values)
        # Add to the total
        integral_sum <- integral_sum + segment_integral

      }else if(mode==1){
        # V2
        segment_integral = integrate(f=func, lower=t1, upper=t2,
                                     subdivisions = subdivisions,
                                     rel.tol = .Machine$double.eps^0.5)
        #print((paste0("segment_integral : ", segment_integral$value)))
        # Add to the total
        integral_sum <- integral_sum + segment_integral$value

      }

    }

  }

  return(data.frame(id = id_df[[id_col]][1], integral = integral_sum))
}




#' from_fd_to_func
#'
#' This funciton transform a fd object into a function.
#' It require either the fd object OR the coefficient and the basis object.
#'
#' @param fd_obj a fd object to transform into a function
#' @param coef the coefficient of an fd object to transform into a function
#' @param basisobj the basis object of an fd object to transform into a function
#'
#' @returns a function
#' @export
#'
#' @importFrom fda eval.fd fd
#'
#' @examples
#' basis = create_bspline_basis(0, 100, 10, 4)
#' coef = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' func_from_fd = from_fd_to_func(coef = coef, basis = basis)
#'
#' @author Francois Bassac
from_fd_to_func <- function(fd_obj = NULL, coef = NULL, basisobj = NULL) {
  if (!is.null(fd_obj)) {
    # If a fd object is given, we extract the components
    coef <- fd_obj$coefs
    basisobj <- fd_obj$basis
  }

  # Safety : check if coef and basis are given
  if (is.null(coef) || is.null(basisobj)) {
    stop("You have to give either one 'fd' object,
           either two arguments 'coef' and 'basisobj'")
  }

  force(coef)
  force(basisobj)

  fd_func <- function(t) {
    fda::eval.fd(evalarg = t, fdobj = fda::fd(coef = coef, basisobj = basisobj))
  }
  return(fd_func)
}


#### Model evaluation ####
#' r_squared_values
#'
#' This function evaluates the R_2 using values y and y_hat
#'
#' @param y a vector of real values
#' @param y_hat a vector of predicted values
#'
#' @returns a value
#' @export
#'
#' @examples
#' y = c(1, 2, 4, 6, 8, 10)
#' y_hat = c(2, 3, 5, 7, 9, 11)
#' r_squared_values(y, y_hat)
#'
#' @author Francois Bassac
r_squared_values <- function(y, y_hat){
  # This function evaluate the r2 determination coefficient from data.
  SSR <- sum((y - y_hat)^2)
  SST <- sum((y - mean(y))^2)
  return(1 - (SSR/SST))
}

#' press_values
#'
#' This function evaluates the press error using values y and y_hat
#'
#' @param y a vector of real values
#' @param y_hat a vector of predicted values
#'
#' @returns a value
#' @export
#'
#' @examples
#' y = c(1, 2, 4, 6, 8, 10)
#' y_hat = c(2, 3, 5, 7, 9, 11)
#' press_values(y, y_hat)
#'
#' @author Francois Bassac
press_values <- function(y, y_hat) {
  return(sum((y - y_hat)^2))
}

#' mae_values
#'
#' This function evaluates the MAE error based on the values.
#'
#' @param y a vector of real values
#' @param y_hat a vector of predicted values
#'
#' @returns a value
#' @export
#'
#' @examples
#' y = c(1, 2, 4, 6, 8, 10)
#' y_hat = c(2, 3, 5, 7, 9, 11)
#' mae_values(y, y_hat)
#'
#' @author Francois Bassac
mae_values <- function(y, y_hat) {
  return(mean(abs(y - y_hat)))
}

#' evaluate_variance_explained
#'
#' This function return the % of variance explained by Y_hat comparing to Y.
#'
#' @param Y a reference value
#' @param Y_hat a modelized value Y_hat = model(X)
#'
#' @returns a value in %
#' @export
#'
#' @importFrom stats var
#'
#' @examples
#' evaluate_variance_explained(c(1,2,3,4,5,6,7,8,9,10),
#' c(0.9, 1.1, 1.9, 2.4, 5.3, 6.01, 7,45, 9.12, 9.04, 11.6))
#'
#' @author Francois Bassac
evaluate_variance_explained <-function(Y, Y_hat){
  # This function return the % of variance explained by Y_hat comparing to Y.
  Var_Exp = 100*(var(Y_hat)/var(Y))
  return(Var_Exp)
}

#' evaluate_results
#' This function evaluates the PRESS, RMSE, MAE, R2 and the % of variance between
#' Y and Y_hat
#'
#' @param Y a vector of real values
#' @param Y_hat a vector of modelized values
#'
#' @returns a dataframe
#' @export
#'
#' @examples
#' evaluate_results(c(1,2,3,4,5), c(0.9, 2.2, 4, 5.5, 5))
#'
#' @author Francois Bassac
evaluate_results <- function(Y, Y_hat){

  Y_real = as.vector(Y)
  Y_hat_2 = as.vector(Y_hat)

  results = data.frame(matrix(ncol=5, nrow=1))
  colnames(results) = c("PRESS", "RMSE", "MAE" , "R2", "var_Y")

  results$PRESS = press_values(Y_real, Y_hat_2)
  results$RMSE = sqrt(press_values(Y_real, Y_hat_2) / length(Y_real))
  results$MAE = mae_values(Y_real, Y_hat_2)
  results$R2 = r_squared_values(Y_real, Y_hat_2)
  results$var_Y = evaluate_variance_explained(Y_real, Y_hat_2)
  return(results)
}

#' evaluate_curves_distances
#'
#' Either for R func or fd function, this function evaluates the distance
#' between the real curve and each curve fun or fd which are in the
#' func_fd_list.
#'
#' @param real_f a fun of fd function, base function to compare
#' @param regul_time a vector of time regularization values
#' @param fun_fd_list a list of fun or fd functions or a fun or a fd function
#'
#' @export
#'
#' @author Francois Bassac
evaluate_curves_distances <- function(real_f, regul_time, fun_fd_list=NULL){
  # either for R func or fd function, this function evaluates the distance
  # between the curves D1 and D2 and the real one.

  if(is.null(fun_fd_list)){
    stop("evaluate_curves_distances_2() : fun_fd_list have to be a list of fun or fd functions.")
  }

  local_from_fun_or_fd_to_fun <- function(f_fd){
    if(mode(f_fd) == "function"){
      f_fd_fun = f_fd
    } else if(inherits(f_fd, "fd")){
      f_fd_fun = from_fd_to_func(f_fd)
    }else{
      stop("evaluate_curves_distances() : Inputs have to be a R function or a fd function.")
    }
    return(f_fd_fun)
  }

  real_fun = local_from_fun_or_fd_to_fun(real_f)

  if(inherits(fun_fd_list, "fd") | mode(fun_fd_list) =="function"){
    # if only one fd or fun, convert to a list
    fun_fd_list = list(fun_fd_list)
  }

  for(i in 1:length(fun_fd_list)){
    delta_i_fun = local_from_fun_or_fd_to_fun(fun_fd_list[[i]])

    diff_i_f <- function(t){
      real_fun(t) - delta_i_fun(t)
    }

    diff_i_ps = pracma::trapz(regul_time, (diff_i_f(regul_time)))
    diff_i_ps_2 = sqrt(pracma::trapz(regul_time, (diff_i_f(regul_time))^2))

    print(paste0("real_f -> curve_",i," / INPROD  : ", diff_i_ps,
                 " / DIST : ",diff_i_ps_2))

  }
}


#' press_model
#' This function evaluates the PRESS error 'LOO' of a lm or glm model.
#'
#' @param model a model from lm or glm
#'
#' @returns a value
#' @export
#'
#' @author Francois Bassac
press_model <- function(model) {
  return(sum((residuals(model)/(1-lm.influence(model)$hat))^2))
}


#### Data visualisation ####

#' eval_max_min_y
#'
#' This function returns the min and max values of a list of functions and fd
#' objects on regul_time.
#'
#' @param f_list a list of functions and fd objects
#' @param regul_time a vector of time evaluation points.
#'
#' @returns a vector
#' @export
#'
#' @author Francois Bassac
eval_max_min_y <- function(f_list, regul_time){
  y_max = 0
  y_min = 0
  for(i in 1:length(f_list)){
    current_f = f_list[[i]]
    if(mode(current_f)=="function"){
      temp_max = max(current_f(regul_time))
      temp_min = min(current_f(regul_time))
    } else if (inherits(current_f, "fd")){
      temp_max = max(eval.fd(evalarg = regul_time, fdobj = current_f))
      temp_min = min(eval.fd(evalarg = regul_time, fdobj = current_f))
    }
    y_max = max(y_max, temp_max)
    y_min = min(y_min, temp_min)
  }
  return(c(y_min, y_max))
}

#' plot_real_and_smoothed_data_ind
#'
#' @param df_wide a wide dataframe, output of convert_to_wide_format()
#' @param df_fd a list of functional data from df_wide
#' @param time_seq a vector to plot on, default 0:100
#' @param id a value of the id to plot
#' @param col_list a list of color to separate the states
#'
#' @returns a plot
#' @export
#'
#' @import graphics
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas,  transition_df)
#' df_processed = cat_data_to_indicatrice(df)
#'
#' df_regul = list()
#' df_wide = list()
#' for(name in names(df_processed)){
#' print(paste0(name, " regularisation"))
#' df_regul[[name]] = regularize_time_series(df_processed[[name]],
#' time_seq =  c(0:100), curve_type = 'cat', id_col='id', time_col='time')
#'
#' df_wide[[name]] = convert_to_wide_format(df_regul[[name]], id_col='id',
#' time_col='time')
#' }
#'
#' basis = create_bspline_basis(0, 100, 10, 4)
#'
#' df_fd = list()
#' for(name in names(df_wide)){
#' print(paste0(name, " fd transformation"))
#' df_fd[[name]] = fda::Data2fd(argvals = c(0:100),
#' y = t(df_wide[[name]][, -c(1)]), basis)
#' }
#'
#'plot_real_and_smoothed_data_ind(df_wide, df_fd, c(0:100), id=1)
#'
#' @author Francois Bassac
plot_real_and_smoothed_data_ind <- function(df_wide, df_fd,
                                            time_seq =0:100, id=1,
                                            col_list=c('blue', 'red',
                                                       'green', 'yellow')){
  graphics::par(mfrow=c(1,length(df_wide)))
  for(i in 1:length(df_wide)){
    plot(time_seq, df_wide[[i]][df_wide[[i]]$id==id, -c(1)], type = 's',
         ylab = paste0("Real and fd n=", i))
    if(i==1){
      graphics::title(paste0("Real and fd values for the individual n=",id))
      print(paste0("Real and fd values for the individual n=",id))
    }
    plot(df_fd[[i]][id], add = TRUE, col = col_list[[i]], lwd = 2)
  }

  graphics::par(mfrow=c(1,1))
}

#' plot_CFD_individuals
#'
#' This function only plot some individuals.
#' It plots the first plot_individuals of df_to_plot.
#' Works both for single state and multistates data (numerical states 1, 2, 3, etc)
#'
#' @param df_to_plot dataframe whose individuals will be plotted
#' @param n_ind_to_plot number of the first individuals to plot, default 5
#' @param id_col col_name of df_to_plot for the id, not character,
#' default id'
#' @param time_col col_name of df_to_plot for the time, not character,
#' default 'time'
#' @param by_cfda a boolean to use cfda package function plotData, default FALSE
#'
#' @returns a ggplot
#' @export
#'
#' @import ggplot2
#' @importFrom cfda plotData
#'
#' @examples
#' df = generate_X_df()
#' plot_CFD_individuals(df, 5)
#' plot_CFD_individuals(df, 5, by_cfda = TRUE)
#'
#' @author Francois Bassac
plot_CFD_individuals <- function(df_to_plot, n_ind_to_plot = 5,
                                 id_col = 'id', time_col='time',
                                 by_cfda=FALSE){
  # This function only plot some individuals.
  # It plots the first plot_individuals of df_to_plot
  df_plot = df_to_plot[df_to_plot$id %in% c(1:n_ind_to_plot), ]
  df_plot$id = as.factor(df_plot$id)

  state_col = setdiff(names(df_to_plot), c(id_col, time_col))

  if(by_cfda == FALSE){
    ggplot2::ggplot(df_plot,
                    ggplot2::aes(x = .data[[time_col]], y = .data[[state_col]],
                                 color = .data[[id_col]])) +
      ggplot2::geom_step() +
      ggplot2::facet_wrap(~ .data[[id_col]], ncol = 1) +
      ggplot2::scale_y_continuous(breaks = c(0, 1)) +
      ggplot2::labs(title = "Synthetic time series (state 0/1)",
                    x = "Time",
                    y = "State") +
      ggplot2::theme_minimal()
  }else if (by_cfda == TRUE){
    cfda::plotData(df_to_plot[df_to_plot$id %in% c(1:n_ind_to_plot),],
                   addBorder = TRUE)
  }
}

#' plot_model_metrics_base
#'
#' This funciton plots some histograms for train_results and test_results
#'
#' @param train_results a dataframe of train results
#' @param test_results a dataframe of test results
#' @param models_to_plot a list of characters of the models to plot, default
#' c("FPLS", "SmoothPLS", "NaivePLS")
#' @param n_digits a integer for the number of significative numbers to print,
#' default 3
#'
#' @returns a plot
#' @export
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @author Francois Bassac
plot_model_metrics_base <- function(train_results, test_results,
                                    models_to_plot = c("FPLS",
                                                       "SmoothPLS", "NaivePLS"),
                                    n_digits = 3) {

  # Étape 1 : ajouter colonne "Set"
  train_r <- train_results
  test_r <- test_results
  train_r$Set <- "Train"
  test_r$Set <- "Test"

  # Étape 2 : ajouter colonne "Model"
  train_r$Model <- rownames(train_r)
  test_r$Model <- rownames(test_r)

  train_r = train_r[train_r$Model %in% models_to_plot , ]
  test_r = test_r[test_r$Model %in% models_to_plot , ]

  # Étape 3 : fusionner
  all_results <- rbind(train_r, test_r)

  # Étape 4 : passer en format long (base R)
  metrics <- c("PRESS", "RMSE", "MAE", "R2")
  long_results <- data.frame()

  for (metric in metrics) {
    temp <- data.frame(
      Model = all_results$Model,
      Set = all_results$Set,
      Metric = metric,
      Value = all_results[[metric]]
    )
    long_results <- rbind(long_results, temp)
  }

  # Étiquette arrondie
  long_results$Label <- format(round(long_results$Value, n_digits), nsmall = n_digits)

  # Étape 5 : boucle pour plot
  for (metric in metrics) {
    metric_data <- long_results[long_results$Metric == metric, ]

    p <- ggplot(metric_data, aes(x = .data$Model,
                                 y = .data$Value,
                                 fill = .data$Set)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = .data$Label), # Ici aussi
                position = position_dodge(width = 0.9),
                vjust = -0.3, size = 3)+
      labs(title = paste("Comparison", metric, "per model"),
           y = metric, x = "Model") +
      scale_fill_manual(values = c("Train" = "#1b9e77", "Test" = "#d95f02")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    print(p)
  }
}

#' plot_fd_list
#'
#' This function plots on the same figure the fd curves from the fd_list by
#' evaluating them on the given regul_time
#'
#' @param fd_list a list of fd objects
#' @param curves_names a list of the curves names
#' @param regul_time a numeric vector
#'
#' @returns a plot
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom fda eval.fd
#'
#' @author Francois Bassac
plot_fd_list <- function(fd_list, curves_names, regul_time){

  fd_df <- data.frame()

  for (i in seq_along(fd_list)) {
    values <- fda::eval.fd(regul_time, fd_list[[i]]) # Keep only the CFD
    temp_df <- data.frame(
      time = regul_time,
      value = as.vector(values),
      name = curves_names[[i]]
    )
    fd_df <- rbind(fd_df, temp_df)
  }
  ggplot(fd_df, aes(x = .data$time, y = .data$value, color = .data$name)) +
    geom_line(linewidth = 1) +
    labs(title = "Regression curves",
         x = "Time", y = "Value", color = "State") +
    theme_minimal() +
    theme(legend.position = "right")
}

#' select_from_fd_list
#'
#' This function return a list of functions of fd_list without the ones
#' specified in the numeric vector toDrop.
#'
#' @param fd_list a list of fd object
#' @param toDrop a numeric vector
#'
#' @returns a list of fd objects
#' @export
#'
#'
#' @author Francois Bassac
select_from_fd_list <- function(fd_list, toDrop = NULL){

  new_fd_list = list()
  toDrop_order = toDrop[order(toDrop)]
  j = 1
  for(i in 1:length(fd_list)){
    if(i != toDrop_order[j]){
      new_fd_list = append(new_fd_list, fd_list[[i]])
    } else {
      j = j+1
    }
  }
  return(new_fd_list)
}

#### Linear Algebra #####

#' block_diag
#'
#' returns a matrix whoses blocks are A and B.
#' returns :
#' ( A, 0 )
#' ( 0, B )
#'
#' @param A a matrix
#' @param B a matrix
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' A = matrix(c(1, 2, 3, 4), 2)
#' B = matrix(c(5, 6, 7, 8,9, 10), 2)
#' C = block_diag(A, B)
#'
#' @author Francois Bassac
block_diag <- function(A, B) {
  # A & B dimensions
  n1 <- nrow(A); m1 <- ncol(A)
  n2 <- nrow(B); m2 <- ncol(B)

  # Init
  result <- matrix(0, nrow = n1 + n2, ncol = m1 + m2)

  # A in top left
  result[1:n1, 1:m1] <- A

  # B in bottom right
  result[(n1 + 1):(n1 + n2), (m1 + 1):(m1 + m2)] <- B

  return(result)
}


#### parallel future ####

#' get_optimal_cores
#' Determine optimal number of cores for parallel processing
#'
#' @param parallel boolean, user request for parallelization, default TRUE
#' @param computational_load numeric, the size of the task
#'
#' @returns an integer representing the number of cores to use
#' @importFrom future availableCores
#' @noRd
get_optimal_cores <- function(parallel = TRUE, computational_load) {
  if (!parallel) {
    return(1)
  }

  base_threshold <- getOption("SmoothPLS.parallel_threshold", default = 2500)
  max_cores <- max(future::availableCores() - 2, 1)

  # 1 coeur pour chaque bloc de 'base_threshold' calculs
  optimal_cores <- floor(computational_load / base_threshold)

  # On ne dépasse pas la limite de la machine
  nb_cores <- min(max_cores, optimal_cores)

  # Si le calcul est trop petit, on force le séquentiel (1 coeur)
  if (nb_cores <= 1) {
    return(1)
  }

  return(nb_cores)
}
