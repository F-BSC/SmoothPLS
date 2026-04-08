

#### Smooth PLS ####

#### Miscalenious ####

#' build_u_ki_list
#'
#' @param N_states a integer, number of different states
#' @param nbComp a integer, max number of components
#' @param ms_pls_models a list of the intermediate pls models to evaluate the
#' real multi-stats pls components.
#'
#' @returns a list of the u_i^k
#'
#' @author Francois Bassac
build_u_ki_list <- function(N_states, nbComp, ms_pls_models){
  u_ki = list()

  for(k in 1:N_states){
    # For each state
    k_state_u_i = list()
    for(i in 1:nbComp){
      # For every components
      k_state_u_i =append(k_state_u_i,
                          list(ms_pls_models[[i]]$loading.weights[k,1]))
    }
    names(k_state_u_i) = paste0("t_",1:nbComp)

    u_ki = append(u_ki, list(k_state_u_i))
  }

  names(u_ki) = paste0("state_", 1:N_states)
  return(u_ki)
}

#### Multivariate ####

#' smoothPLS_multi
#'
#' This function performs the smooth PLS algorithm for both the univariate case
#' for a Scalar Functionnal Data or a Categorical Functional Data and the
#' multivariate case for a mix of functional data of different nature
#' (SFD or CFD).
#' For some input, if the same value is needed for all the different curves,
#' no need to make a list sith the value per curve.
#'
#' @param df_list a list of dataframe (id, time, value_or_state)
#' @param Y a vector of the scalar response
#' @param basis_obj a basis fd object or a list of basis fd object
#' @param regul_time_obj a vector for time regularization values or a list
#' @param curve_type_obj a list or vector of the differents curves types,
#' 'cat' or 'num.
#' @param orth_obj a list or a vector of booleans if the orthonormalization
#' is needed
#' @param id_col_obj a list or a vector of the id column name
#' @param time_col_obj a list or a vector of time column name
#' @param int_mode a integer of the integration method : 1 for integrate,
#' 2 for pracma::trapz
#' @param print_steps a boolean to print the algorithm steps
#' @param plot_rmsep a boolean to plot the pls model RMSEP
#' @param print_nbComp a boolean to print the optimal number of components
#' @param plot_reg_curves a boolean to plot the regressions curves
#' @param jackknife a boolean for the jackknife input of pls() function,
#' default TRUE
#' @param validation a character for the validation input of pls() function,
#' default 'LOO'
#'
#' @returns a list of the plsr_model and the regression curves (and intercept).
#' @export
#'
#' @importFrom pls plsr RMSEP
#' @importFrom stats coef
#' @importFrom graphics plot title
#'
#' @author Francois Bassac
smoothPLS <- function(df_list, Y,
                      basis_obj, regul_time_obj = NULL,
                      curve_type_obj, orth_obj = TRUE,
                      id_col_obj = 'id', time_col_obj = 'time', int_mode = 1,
                      print_steps = FALSE,
                      plot_rmsep = TRUE, print_nbComp = TRUE,
                      plot_reg_curves = FALSE,
                      jackknife = TRUE,
                      validation = 'LOO'){

  if(print_nbComp){cat("### Smooth PLS ### \n")}
  # 1 assertion

  if(print_steps){
    cat("=> Input format assertions.\n")
  }
  assert_obj = assert_multivariate_smoothPLS_inputs(
    df_list = df_list,
    Y = Y,
    basis_obj = basis_obj,
    regul_time_obj = regul_time_obj,
    curve_type_obj = curve_type_obj,
    orth_obj = orth_obj,
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

  # 2 orthogonalize
  if(print_steps){
    cat("=> Orthonormalize basis.\n")
  }

  orth_basis_list = orthonormalize_basis_list(basis_list = basis_list,
                                              orth_list = orth_list)


  # 3 build df_processed_list and curves_names_list
  if(print_steps){
    cat("=> Data objects formating.\n")
  }

  if(N_curves == 1 && mode(df_list[[1]]) != 'list' && ncol(df_list) == 3){
    df_list = list(df_list)
  }

  new_list_obj = build_new_data_list(df_list = df_list,
                                     N_curves = N_curves,
                                     orth_basis_list = orth_basis_list,
                                     basis_list = basis_list,
                                     curve_type_list = curve_type_list,
                                     id_col_list = id_col_list,
                                     time_col_list = time_col_list,
                                     regul_time_list = regul_time_list)

  df_processed_list = new_list_obj$df_processed_list
  curves_names_list = new_list_obj$curves_names_list
  new_curves_type_list = new_list_obj$new_curves_type_list
  new_basis_list = new_list_obj$new_basis_list
  new_orth_basis_list = new_list_obj$new_orth_basis_list
  new_id_col_list = new_list_obj$new_id_col_list
  new_time_col_list = new_list_obj$new_time_col_list
  new_regul_time_list = new_list_obj$new_regul_time_list

  # 5 Build Lambda matrix
  if(print_steps){
    cat("=> Evaluate Lambda matrix.\n")
  }

  for(i in 1:length(df_processed_list)){

    if(print_steps){
      cat(paste0("==> Lambda for : ", curves_names_list[i], ".\n"))
    }


    lambda = evaluate_lambda(df = df_processed_list[[i]],
                             basis = new_orth_basis_list[[i]],
                             curve_type = new_curves_type_list[[i]],
                             int_mode = int_mode,
                             id_col = new_id_col_list[[i]],
                             time_col = new_time_col_list[[i]],
                             regul_time = new_regul_time_list[[i]])
    #dim(lambda)
    if(i != 1){
      Lambda = cbind(Lambda, lambda)
    }else{
      Lambda = lambda
    }
  }

  # 6 PLSR model
  if(print_steps){
    cat("=> PLSR model.\n")
  }

  plsr_model = pls::plsr(Y ~ as.matrix(Lambda), validation = validation,
                    jackknife = jackknife, intercept = TRUE,
                    center = TRUE)

  # 7 optimal number of components
  nb_comp_pls_opt = which.min(pls::RMSEP(plsr_model)$val[1, , ])-1

  if(print_nbComp){
    cat(paste0("=> Optimal number of PLS components : ", nb_comp_pls_opt, "\n"))
  }

  # Assert a non equal to zero optimal number of component
  if(nb_comp_pls_opt == 0){
    cat("No Optimal number of component!\n")
    cat("Search the optimal number of component EXCLUDING the intercept.\n")
    nb_total_cp = ncol(plsr_model$scores)
    nb_comp_pls_opt = which.min(
      pls::RMSEP(plsr_model)$val[1, , ][c(2:nb_total_cp)]
    )
    cat(paste0("New \'Optimal\' number of PLS components :",
               nb_comp_pls_opt, ".\n"))
  }

  if(plot_rmsep){
    plot(pls::RMSEP(plsr_model),
         main="Cross Validation",
         xlab="Components number")
  }

  # Coefficients
  # b_i : Y = sum b_i lambda_i
  b_i = stats::coef(plsr_model, ncomp = nb_comp_pls_opt)
  # c_i : Y = sum c_i t_i
  c_i = plsr_model$Yscores
  # d_i : Lambda = sum d_i t_i
  d_i = plsr_model$loadings
  # u_i : t_i = Lambda_{i-1} u_i
  u_i = plsr_model$loading.weights

  # Functions
  if(print_steps){
    cat("=> Evaluate SmoothPLS functions and <w_i, p_j> coef.\n")
  }
  functions_gamma_list = build_spls_functions(
    curves_names_list = curves_names_list,
    new_basis_list = new_basis_list,
    new_orth_basis_list = new_orth_basis_list,
    d_i = d_i, u_i = u_i
  )

  p_hij_list = functions_gamma_list$p_hij_list
  w_hij_list = functions_gamma_list$w_hij_list
  gamma_i_list = functions_gamma_list$gamma_i_list
  v_i_list = functions_gamma_list$v_i_list

  # Delta
  if(print_steps){
    cat("=> Build regression functions and intercept.\n")
  }

  delta_spls = build_reg_curve_spls(plsr_model = plsr_model,
                                    curves_names_list = curves_names_list,
                                    v_i_list = v_i_list,
                                    nb_comp_pls_opt = nb_comp_pls_opt,
                                    print_steps = print_steps)

  if(plot_reg_curves){
    for(i in 2:length(delta_spls)){
      plot(delta_spls[[i]])
      title(curves_names_list[[i-1]])
    }
  }

  spls_obj = list(plsr_model, nb_comp_pls_opt, curves_names_list,
                  v_i_list, delta_spls)
  names(spls_obj) = c("plsr_model", "nbCP_opti", "curves_names",
                      "v_i_list", "reg_obj")

  return(spls_obj)

}

#' assert_multivariate_smoothPLS_inputs
#'
#' This function checks the integrity of the input for multivariate_fpls.
#' It returns a list of (basis_list, regul_time_list, curve_type_list,
#' id_col_list, time_col_list)
#'
#' @param df_list a list of dataframes (id, time, value_or_state)
#' @param Y a numeric vector of the response
#' @param basis_obj a list of basis object or a basis object
#' @param regul_time_obj a vector of time regularization values or a list of
#' vectors
#' @param curve_type_obj a character "cat" or 'num' or a list of those values
#' @param orth_obj a boolean, a list or a vector of boolean to orthonormalize
#' or not a basis
#' @param id_col_obj a character of the id column for all the curves or a list
#' of id column character
#' @param time_col_obj a character of the time column for all the curves or a
#' list of time column character
#'
#'
#' @returns a list of (basis_list, regul_time_list, curve_type_list, orth_list,
#' id_col_list, time_col_list)*
#' @export
#'
#' @author Francois Bassac
assert_multivariate_smoothPLS_inputs <- function(df_list, Y, basis_obj,
                                                 regul_time_obj = NULL,
                                                 curve_type_obj = NULL,
                                                 orth_obj = list(TRUE),
                                                 id_col_obj = 'id',
                                                 time_col_obj = 'time'){
  # N_curves
  if(length(df_list) != 3){
    N_curves = length(df_list)
  }else if(length(df_list) == 3){
    # if 3, either 1 curve or 3
    if(mode(df_list[[1]]) == "numeric" || mode(df_list[[1]]) == "character" ){
      N_curves = 1
    }else{
      N_curves = 3
    }
  }

  if(mode(df_list) != "list" && N_curves != 1){
    stop("smoothPLS_multi() : df_list have to be a list of dataframe
         (id, time, value).")
  }else if(mode(df_list) == "list" && N_curves != 1 &&
           !all(sapply(df_list, function(x) ncol(x)==3 ))){
    stop("smoothPLS_multi() : all df of df_list should have 3 columns,
         (id, time, value).")
  }else if(N_curves == 1 && mode(df_list[[1]]) != 'list' && ncol(df_list) != 3){
    stop("naivePLS() : df_list have to be a list of dataframe
         (id, time, value).")
  }

  # Y
  if((!(is.vector(Y) && mode(Y)=="numeric") &
      !(is.list(Y) && all(sapply(Y, function(x) inherits(x, "numeric")))))){
    stop("smoothPLS_multi() : Y should be a numeric vector.")
  }

  # basis_list
  if(inherits(basis_obj, "basisfd")){
    basis_list = obj_list_creation(N_rep = N_curves, obj = basis_obj)
  }else if(mode(basis_obj) == "list" &&
           all(sapply(basis_obj, function(x) inherits(x, "basisfd")))){
    basis_list = basis_obj
  } else {
    stop("smoothPLS_multi() : basis_obj have to be a basis fd object or a list
          of basis fd objects.")
  }

  # regul_time_list
  #if(is.null(regul_time_obj)  &&
  #   (curve_type_obj == 'cat' || curve_type_obj[[1]] == 'cat')){
  #  regul_time_obj = c(
  #    basis_list[[1]]$rangeval[1]:basis_list[[1]]$rangeval[2]
  #    )
  #}

  # regul_time_list assertion update
  if (is.null(regul_time_obj) && any(unlist(curve_type_obj) == 'cat')) {
    # If no regularization time is provided but at least one curve is categorical,
    # we default to a sequence based on the range of the first basis.
    regul_time_obj <- seq(basis_list[[1]]$rangeval[1],
                          basis_list[[1]]$rangeval[2],
                          by = 1)
  }
  if(mode(regul_time_obj) == "numeric"){
    # Because the functions basis_list_creation works :
    regul_time_list =  obj_list_creation(N_rep = N_curves, regul_time_obj)
  } else if(mode(regul_time_obj) == "list" &&
            all(sapply(regul_time_obj, function(x) inherits(x, "numeric")))){
    regul_time_list = regul_time_obj
  } else {
    stop("smoothPLS_multi() : regul_time_obj have to be a numeric sequence or
         a list of numeric sequence.")
  }

  if(is.null(curve_type_obj)){
    stop("smoothPLS_multi() : curve_type_obj have to be either :
    'cat' or 'num' to be applied to all the curves OR
     a list of ('cat', 'num', ....) to be used for each curve.")
  } else if(all(curve_type_obj == 'cat')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(all(curve_type_obj == 'num')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(!(mode(curve_type_obj) == 'list' &&
              all(sapply(curve_type_obj, function(x) (x=='cat' | x=='num'))))){
    stop("smoothPLS_multi() : curve_type_obj have to be either \'cat\' or
         \'num\', for all curves, or a list  or a list specify 'cat' or 'num'
         for each curve.")
  } else {
    curve_type_list = curve_type_obj
  }

  # orth_list
  if(length(orth_obj) == 1){
    orth_list = obj_list_creation(N_rep = N_curves, orth_obj)
  } else if(length(orth_obj) != N_curves){
    stop("smoothPLS_multi() : length(orth_obj) != N_curves. Should be TRUE/FALSE
         or a list of rigth dimension.")
  } else {
    orth_list = orth_obj
  }

  # id_col_obj & time_col_obj
  if(mode(id_col_obj) == "character"){
    # Because the functions basis_list_creation works :
    id_col_list = obj_list_creation(N_rep = N_curves, id_col_obj)
  } else if(mode(id_col_obj) == "list" &&
            all(sapply(id_col_obj, function(x) inherits(x, "character")))){
    id_col_list = id_col_obj
  } else {
    stop("smoothPLS_multi() : id_col_obj have to be a character or
         a list of characters.")
  }

  if(mode(time_col_obj) == "character"){
    # Because the functions basis_list_creation works :
    time_col_list = obj_list_creation(N_rep = N_curves, time_col_obj)
  } else if(mode(time_col_obj) == "list" &&
            all(sapply(time_col_obj, function(x) inherits(x, "character")))){
    time_col_list = time_col_obj
  } else {
    stop("smoothPLS_multi() : time_col_obj have to be a character or
         a list of characters.")
  }

  # last check : number!
  if(length(basis_list) != N_curves){
    stop("smoothPLS_multi() : basis_obj should be only one basis or a list of
         length the number of curves.")
  }

  if(length(regul_time_list) != N_curves){
    stop("smoothPLS_multi() : regul_time_obj should be only one regul_time
         vector or a list of length the number of curves.")
  }

  if(length(curve_type_list) != N_curves){
    stop("smoothPLS_multi() : curve_type_obj should be only one curve_type
         'cat' or 'num', or a list of length the number of curves.")
  }

  if(length(id_col_list) != N_curves){
    stop("smoothPLS_multi() : id_col_obj should be only one character,
         or a list of length the number of curves.")
  }

  if(length(time_col_list) != N_curves){
    stop("smoothPLS_multi() : time_col_obj should be only one character,
         or a list of length the number of curves.")
  }

  assert_objs = list(N_curves, basis_list, regul_time_list, curve_type_list,
                     orth_list, id_col_list, time_col_list)
  names(assert_objs) = c("N_curves", "basis_list", "regul_time_list",
                         "curve_type_list", "orth_list",
                         "id_col_list", "time_col_list")

  return(assert_objs)
}

#' orthonormalize_basis_list
#'
#' This function orthonormalized a basis from a list if needed.
#'
#' @param basis_list a list of basis to orthonormalized
#' @param orth_list a list of boolean per basis if orthonormalization is needed
#' @param tol a float, tolerance parameter.
#'
#' @returns a list of list of fd object : the orthonormalized functions.
#' @export
#'
#' @author Francois Bassac
orthonormalize_basis_list <- function(basis_list, orth_list, tol=1e-9){
  # This function orthonormalize the basis of the basis list based on the
  # information given by orth_list.
  orth_basis_list <- vector("list", length(basis_list))

  for (i in 1:length(basis_list)) {
    if (orth_list[[i]]) {
      # Spectral method call
      orth_basis_list[[i]] <- gram_schmidt_orthonormalize(
        basis = basis_list[[i]],
        output_type = "fdlist",
        tol = tol
      )
    } else {
      orth_basis_list[[i]] <- from_basis_to_fdlist(basis_list[[i]])
    }
  }

  return(orth_basis_list)
}


#### Lambda #####

#' evaluate_lambda
#'
#' This function evaluates the Lambda matrix such as per column :
#' Lambda_i = int_0^T X(t) phi_i(t) dt.
#' The curve_type input is important function of the type of data you work with
#' 'cat' for Categorical Funcitonal Data
#' 'num' for Scalar Functional Data
#'
#' @param df dataframe X(t)
#' @param basis basis fd object
#' @param curve_type a character, 'cat' for Categorical FD, 'num' for Scalar FD
#' @param int_mode integration mode, 1 for integrate, 2 for pracma::trapz
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#' @param regul_time regul_time a vector of time regularization values default
#'  basis rangeval per 1
#'
#' @returns a matrix \eqn{\Lambda_i = \int_0^T X(t) \phi_i(t) dt}
#' @export
#'
#' @examples
#' df = generate_X_df(nind=100, start=0, end=100, curve_type = 'cat',
#' lambda_0=0.2, lambda_1=0.1, prob_start=0.5)
#' basis = create_bspline_basis(0, 100, 10, 4)
#' Lambda = evaluate_lambda(df, basis, curve_type = 'cat')
#'
#'
#' df = generate_X_df(nind=100, start=0, end=100, curve_type = 'num')
#' basis = create_bspline_basis(0, 100, 10, 4)
#' Lambda = evaluate_lambda(df, basis, curve_type = 'num', int_mode = 2)
#'
#' @author Francois Bassac
evaluate_lambda <- function(df, basis, curve_type = NULL, int_mode = 1,
                            id_col = 'id', time_col = 'time',
                            nb_pt = 10, subdivisions = 100,
                            regul_time = seq(basis$rangeval[1],
                                             basis$rangeval[2],
                                             1)){
  # This function evaluate the Lambda matrix depending of its curve_type

  if(is.null(curve_type)){
    stop("evaluate_lambda() : curve_type should be 'cat' or 'num'.")
  }else if(curve_type == 'cat'){

    Lambda = evaluate_lambda_CFD(df = df, basis = basis, int_mode = int_mode,
                                 id_col = id_col, time_col = time_col,
                                 nb_pt = nb_pt, subdivisions = subdivisions)

  }else if(curve_type == 'num'){
    if(int_mode==1){
      cat("evaluate_lambda_SFD() : int_mode to 2 for pracma::trapz for
        integration stability.\n")
      mode_int=2
    }else{
      mode_int = int_mode
    }

    Lambda = evaluate_lambda_SFD(df = df, basis = basis,
                                 regul_time = regul_time,
                                 int_mode = mode_int,
                                 id_col = id_col, time_col = time_col,
                                 subdivisions = subdivisions)

  }else{
    stop("evaluate_lambda() : curve_type should be 'cat' or 'num'.")
  }
  return(Lambda)
}

#' evaluate_lambda_CFD
#'
#' This function evaluates the Lambda matrix
#' Lambda_ij = int X_j(t) phi_i(t) dt
#'
#' @param df dataframe X(t)
#' @param basis basis fd object or a list of fd functions
#' @param int_mode integration mode, 1 for integrate, 2 for pracma::trapz
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#'
#' @returns a matrix of dimension nbasis columns and nind rows
#' @importFrom fda fd
#'
#' @author Francois Bassac
evaluate_lambda_CFD <- function(df, basis, int_mode = 1,
                                id_col = 'id', time_col = 'time',
                                nb_pt = 10, subdivisions = 100){
  # This function evaluate the matrix LAMBDA for step 1.
  # Lambda is a matrix of p columns (Lambda_i) and n rows (individuals)
  # Lambda_i = int X(t) phi_i(t) dt

  ids = unique(df[[id_col]])
  n_ind = length(ids)

  # Convert in fd list if necessary
  if (inherits(basis, "basisfd")) {
    nbasis <- basis$nbasis
    rangeval <- basis$rangeval
    coefs <- diag(nbasis)
    fd_list <- list()
    for (i in 1:nbasis) {
      fd_list[[i]] <- fda::fd(coefs[, i], basis)
    }
  } else if (is.list(basis) &&
             all(sapply(basis, function(x) inherits(x, "fd")))) {
    fd_list <- basis
    nbasis <- length(fd_list)
  }else {
    stop("evaluate_lambda_CFD() : basis have to be a basisobj from fda
         or a list of fd objects")
  }
  n_col = nbasis

  # Use other functions
  # fd_list = = from_basis_to_fdlist(basis)
  # n_col = basis_list[[1]]$basis$nbasis

  Lambda = data.frame(matrix(nrow=n_ind, ncol=n_col))
  # name creation
  lambda_names = c()
  for(i in 1:n_col){
    name = paste0("Lambda_",i)
    lambda_names = c(lambda_names, name)
  }
  colnames(Lambda) = lambda_names

  # Iteration per basis function per individual
  # LONG
  for(j in 1:n_col){
    # phi_i function definition for integration
    local_phi_fd <- fd_list[[j]]

    local_phi <- from_fd_to_func(local_phi_fd)

    for(i in 1:n_ind){
      df_id = df[df[[id_col]]==ids[i], ]
      int_df = evaluate_id_func_integral(id_df = df_id,
                                         func = local_phi,
                                         mode = int_mode,
                                         id_col = id_col,
                                         time_col = time_col,
                                         nb_pt = nb_pt,
                                         subdivisions = subdivisions)
      Lambda[i, j] = int_df[, 2]
    }
  }
  return(Lambda)
}

#' evaluate_lambda_SFD
#'
#' This function evaluates the Lambda matrix
#' Lambda_ij = int X_j(t) phi_i(t) dt for step > 1
#'
#' @param df dataframe X(t)
#' @param basis basis fd object
#' @param regul_time a vector of time regularization values default basis
#' rangeval per 1
#' @param int_mode int, integration mode, 1 for integrate, 2 for pracma::trapz
#' @param id_col default name of the id column
#' @param time_col default name of the time column
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#'
#' @returns a matrix of dimension nbasis columns and nind rows
#'
#' @importFrom fda eval.fd fd
#' @importFrom stats approxfun integrate
#' @importFrom pracma trapz
#'
#' @author Francois Bassac
evaluate_lambda_SFD <- function(df, basis,
                                regul_time = seq(basis$rangeval[1],
                                                 basis$rangeval[2],
                                                 1), int_mode = 1,
                                id_col = 'id', time_col = 'time',
                                subdivisions = 100){
  # This function evaluate the matrix LAMBDA for step 1.
  # Lambda is a matrix of p columns (Lambda_i) and n rows (individuals)
  # Lambda_i = int X(t) phi_i(t) dt

  ids = unique(df[[id_col]])
  n_ind = length(ids)

  basis_list = from_basis_to_fdlist(basis)

  n_col = basis_list[[1]]$basis$nbasis

  # if df has more than 3 columns, we remove state
  if(ncol(df) > 3){
    stop("evaluate_lambda_SFD() : df shoul have only 3 column, id, time, value")
  }

  value_col = setdiff(names(df), c(id_col, time_col))

  Lambda = data.frame(matrix(nrow=n_ind, ncol=n_col))
  # name creation
  lambda_names = c()
  for(i in 1:n_col){
    name = paste0("Lambda_",i)
    lambda_names = c(lambda_names, name)
  }
  colnames(Lambda) = lambda_names

  #eval_basis = matrix(nrow = basis$nbasis, ncol = length(regul_time))

  df_regul = regularize_time_series(df, time_seq = regul_time,
                                    curve_type = 'num',
                                    id_col = id_col, time_col = time_col)

  # Iteration per basis function per individual
  # LONG
  for(j in 1:n_col){

    # 1/ evaluate the j^th basis on the regul_time
    eval_basis = as.vector(
      fda::eval.fd(evalarg = regul_time, fdobj = basis_list[[j]])
    )

    for(k in 1:n_ind){

      df_id = df_regul[df_regul[[id_col]]==ids[k], ]
      # Order
      #df_id <- df_id[order(df_id[[time_col]]), ]

      # 2/ multiply by the X_1(t) values
      prod_values = eval_basis*df_id[[value_col]]

      # Integration
      if(int_mode == 1){
        # integrate
        local_fun = stats::approxfun(x = regul_time, y = prod_values)

        integral_sum = stats::integrate(f=local_fun,
                                 lower=regul_time[1],
                                 upper=regul_time[length(regul_time)],
                                 subdivisions = subdivisions)[1]$value

      } else if(int_mode == 2){
        # pracma::trapz
        integral_sum = pracma::trapz(x = regul_time,
                                     y = prod_values)

      }

      Lambda[k, j] = integral_sum
    }
  }
  return(Lambda)
}


#' build_new_data_list
#'
#' This function preprocess the different input in order to format them to the
#' right number of curve. Warning the "right" number of curve take into account
#' the number of different states for the CFDs.
#'
#' @param df_list a list of dataframes (id, time, value_or_state)
#' @param N_curves a integer, the number of curves
#' @param orth_basis_list a list of orthogonolized basis fd list
#' @param basis_list a list of basis fd functions
#' @param curve_type_list a list of the curve type of each curve
#' @param id_col_list a list of the id column name for each curve
#' @param time_col_list a list of the time column name for each curve
#' @param regul_time_list a list of the time regularization vector for each
#' curve
#'
#'
#' @returns a list
#'
#' @author Francois Bassac
build_new_data_list <- function(df_list, N_curves,
                                orth_basis_list = NULL,
                                basis_list = NULL,
                                curve_type_list,
                                id_col_list, time_col_list,
                                regul_time_list){

  df_processed_list = list()
  curves_names_list = list()
  new_curves_type_list = list()
  new_basis_list = list()
  new_orth_basis_list = list()
  new_id_col_list = list()
  new_time_col_list = list()
  new_regul_time_list = list()

  N_curves = length(curve_type_list)

  for(i in 1:N_curves){
    if(curve_type_list[[i]] == 'cat'){

      state_col = setdiff(names(df_list[[i]]),
                          c(id_col_list[[i]], time_col_list[[i]]))
      states = unique(df_list[[i]][state_col])

      if(all(states[[state_col]] %in% c(0, 1)) ||
         all(states[[state_col]] %in% c('0', '1')) ){
        # case : one state CFD in indicatrice form

        df_processed = list(df_list[[i]])
        names(df_processed) = paste0(state_col, "_1")
      } else {
        # case : multistates CFD
        df_processed = cat_data_to_indicatrice(df_list[[i]],
                                               id_col = id_col_list[[i]],
                                               time_col = time_col_list[[i]])
      }

      df_processed_list = append(df_processed_list, df_processed)

      name_df = paste0('CatFD_', i, "_", names(df_processed))
      curves_names_list = c(curves_names_list, name_df)

      nb_states = length(name_df)

      curve_types = obj_list_creation(nb_states, 'cat')
      new_curves_type_list = append(new_curves_type_list, curve_types)

      if(is.null(orth_basis_list)){
        new_orth_basis_list = append(new_orth_basis_list, list(NULL))
      }else{
        new_orth_basis_list = append(new_orth_basis_list,
                                   obj_list_creation(nb_states,
                                                     orth_basis_list[[i]]))
      }

      if(is.null(basis_list)){
        new_basis_list = append(new_basis_list, list(NULL))
      }else{
        new_basis_list = append(new_basis_list,
                                obj_list_creation(nb_states,
                                                  basis_list[[i]]))
      }

      new_id_col_list = append(new_id_col_list,
                               obj_list_creation(nb_states, id_col_list[[i]]))
      new_time_col_list = append(new_time_col_list,
                                 obj_list_creation(nb_states,
                                                   time_col_list[[i]]))
      new_regul_time_list = append(new_regul_time_list,
                                   obj_list_creation(nb_states,
                                                     regul_time_list[[i]]))

    }else if(curve_type_list[[i]] == 'num'){

      #if(N_curves == 1){
      #  df_processed_list = append(df_processed_list, list(df_list[[i]]))
      #}else{
      #  df_processed_list = append(df_processed_list, list(df_list[[i]]))
      #}
      df_processed_list = append(df_processed_list, list(df_list[[i]]))

      value_col = setdiff(names(df_list[[i]]),
                          c(id_col_list[[i]], time_col_list[[i]]))


      name_df = paste0('NumFD_', i, "_", value_col)
      curves_names_list = c(curves_names_list, name_df)

      new_curves_type_list = append(new_curves_type_list, list('num'))

      if(is.null(orth_basis_list)){
        new_orth_basis_list = append(new_orth_basis_list, list(NULL))
      }else{
        new_orth_basis_list = append(new_orth_basis_list,
                                     list(orth_basis_list[[i]]))
      }

      if(is.null(basis_list)){
        new_basis_list = append(new_basis_list, list(NULL))
      }else{
        new_basis_list = append(new_basis_list, list(basis_list[[i]]))
      }
      new_id_col_list = append(new_id_col_list, list(id_col_list[[i]]))
      new_time_col_list = append(new_time_col_list, list(time_col_list[[i]]))
      new_regul_time_list = append(new_regul_time_list,
                                   list(regul_time_list[[i]]))
    }
  }
  names(df_processed_list) = curves_names_list

  new_list_obj = list(df_processed_list,
                      curves_names_list,
                      new_curves_type_list,
                      new_basis_list,
                      new_orth_basis_list,
                      new_id_col_list,
                      new_time_col_list,
                      new_regul_time_list)
  names(new_list_obj) = c("df_processed_list", "curves_names_list",
                          "new_curves_type_list", "new_basis_list",
                          "new_orth_basis_list", "new_id_col_list",
                          "new_time_col_list", "new_regul_time_list")

  return(new_list_obj)
}

#' build_spls_functions
#'
#' This function build some intermediate functions of the smooth pls algorithm.
#' It also evaluates the different coefficients gamma_ij
#'
#'
#' @param curves_names_list a list of the names of the different curves
#' @param new_basis_list a list of the initial basis fd object for all curves
#' @param new_orth_basis_list a list of the orthogonalized basis as fd list for
#' all curves
#' @param d_i the pls coefficient such as X = d_i t_i :
#' plsr_model$loadings
#' @param u_i the pls coefficient such as t_i = X u_i :
#' plsr_model$loading.weigths
#'
#' @returns a list Lambda = sum d_i t_i
#' @export
#'
#' @author Francois Bassac
build_spls_functions <- function(curves_names_list,
                                 new_basis_list,
                                 new_orth_basis_list,
                                 d_i, u_i){

  N_curves_processed = length(curves_names_list)
  p_hij_list = p_w_building(coefficient = d_i,
                            N_curves_processed = N_curves_processed,
                            new_basis_list = new_basis_list,
                            new_orth_basis_list = new_orth_basis_list,
                            curves_names_list = curves_names_list)

  w_hij_list = p_w_building(coefficient = u_i,
                            N_curves_processed = N_curves_processed,
                            new_basis_list = new_basis_list,
                            new_orth_basis_list = new_orth_basis_list,
                            curves_names_list = curves_names_list)

  # gamma_ij
  gamma_i_list = list() #depend only on the curve
  # for all curves
  for(i in 1:N_curves_processed){
    g_ih = evaluate_gamma_ij(w_i_list = w_hij_list[[i]],
                             p_i_list = p_hij_list[[i]])
    gamma_i_list = append(gamma_i_list, list(g_ih))
  }
  names(gamma_i_list) = curves_names_list

  # v_ij list
  v_i_list = list()
  # for all curves
  for(i in 1:N_curves_processed){
    v_i = evaluate_V_i_function(w_i_list = w_hij_list[[i]],
                                gamma_ij_list = gamma_i_list[[1]])
    v_i_list = append(v_i_list, list(v_i))
  }
  names(v_i_list) = curves_names_list

  functions_gamma_list = list(p_hij_list, w_hij_list, gamma_i_list, v_i_list)
  names(functions_gamma_list) = c("p_hij_list", "w_hij_list",
                                  "gamma_i_list", "v_i_list")

  return(functions_gamma_list)
}

#' p_w_building
#'
#' This functions evalutates the p_i(t) (X_i regression functions) and w_i(t)
#' such as \eqn{t_i = \int_0^T w_i(t)X_{i-1}(t) dt}.
#' The evaluation is done for all the components.
#'
#' @param coefficient a table of the coefficient to use
#' @param N_curves_processed a integer, the real number of curves
#' (1 sfd = 1 curve, 1 cfd = 1curve per state)
#' @param new_basis_list a list of processed basis list
#' (1 basis per curve_processed)
#' @param new_orth_basis_list a list of orthonormalized basis list
#' (1 basis per curve_processed)
#' @param curves_names_list a list of curve name (1 name per curve_processed)
#'
#' @returns a list of fd objects
#'
#' @importFrom fda fd
#'
#' @author Francois Bassac
p_w_building <- function(coefficient,
                         N_curves_processed,
                         new_basis_list,
                         new_orth_basis_list,
                         curves_names_list){
  pw_hij_list = list() # coord 1 : Curve, coord 2 : pls step
  base = 0

  nbComp = dim(coefficient)[2]

  # for each curve
  for(i in 1:N_curves_processed){
    pw_ij_list = list()
    #cat(paste0("= FD Curve : ", i ,"\n"))
    n_i_basis = new_basis_list[[i]]$nbasis
    pw_i = fda::fd(coef = rep(0, n_i_basis), basisobj = new_basis_list[[i]])
    # for each pls step
    for(h in 1:nbComp){
      #cat(paste0("== PLS step : ", h ,"\n"))
      # for every basis function
      for(j in 1:n_i_basis){
        # cat(paste0("=== Global coordinate : ", c(base + j) ,"\n"))
        pw_i = pw_i + coefficient[c(base + j),h] * new_orth_basis_list[[i]][[j]]
      }
      pw_ij_list = append(pw_ij_list, list(pw_i))
    }
    base = base + n_i_basis
    names(pw_ij_list) = paste0("PLS_Step_", 1:nbComp)
    pw_hij_list = append(pw_hij_list, list(pw_ij_list))
  }

  pw_hij_names = curves_names_list
  names(pw_hij_list) = pw_hij_names

  return(pw_hij_list)
}

#' evaluate_gamma_ij
#'
#' This function evaluates int_0^T w_i_fd p_j_t dt for all j
#' return a list, i^th element length is (i-1)
#'
#' @param w_i_list a list of fd functions w_i(t)
#' @param p_i_list a list of fd functions p_i(t)
#'
#' @returns a list of list of the gamma_ij values
#'
#' @importFrom fda inprod
#' @export
#'
#' @author Francois Bassac
evaluate_gamma_ij <- function(w_i_list, p_i_list){
  # This function evaluates int_0^T w_i_fd p_j_t dt for all j
  # return a list, i^th element length is (i-1)

  gamma_ij_list = list()

  for(i in 1:length(w_i_list)){
    w_i_fd = w_i_list[[i]]

    if(i == 1){
      gamma_ij = 0
      gamma_ij_list = append(gamma_ij_list, list(gamma_ij))
    }else{
      gamma_ij = list()

      for(j in 1:(i-1)){
        # gamma_ij evaluation by fda::inprod

        p_j_t_fd = p_i_list[[j]]
        # if(int_mode = 2){
        #   w_eval = fda::eval.fd(regul_time, w_i_fd)
        #   p_eval = fda::eval.fd(regul_time, p_i_fd)
        #   gamma_ij_temp = w_eval * p_eval
        #}else{ }
        gamma_ij_temp = fda::inprod(w_i_fd, p_j_t_fd)

        gamma_ij = append(gamma_ij, list(gamma_ij_temp))
      }
      names(gamma_ij) = as.character(1:(i-1))

      gamma_ij_list = append(gamma_ij_list, list(gamma_ij))
    }

  }
  names(gamma_ij_list) = as.character(1:length(w_i_list))

  return(gamma_ij_list)
}

#' evaluate_V_i_function
#'
#' This function evaluates the different functions v_i_t bases on gamma_ij and
#' the functions w_i_fd. Recursive :
#' \eqn{v_i = W_i - \sum_{j=1}^{i-1}gamma_ij*v_j}
#'
#' @param w_i_list a list of fd functions w_i(t)
#' @param gamma_ij_list a list of list of gamma_ij values
#'
#' @returns a list of fd functions
#' @export
#'
#' @author Francois Bassac
evaluate_V_i_function <- function(w_i_list, gamma_ij_list){
  # This function evaluates the different functions v_i_t bases on gamma_ij and
  # the functions w_i_fd
  # recursive : v_i = W_i - \sum_{j=1}^{i-1}gamma_ij*v_j
  # return fd object

  V_list = list()
  for(i in 1:length(w_i_list)){
    # w_i_fd
    w_i_fd = w_i_list[[i]]
    if(i==1){
      V_i = w_i_fd
    } else {
      V_i = w_i_fd
      for(j in 1:(i-1)){
        V_i = V_i - gamma_ij_list[[i]][[j]] * V_list[[j]]
      }
    }
    V_list = append(V_list, list(V_i))
  }
  names(V_list) = as.character(1:length(w_i_list))
  return(V_list)
}

#' evaluate_delta_PLS
#'
#' This functions evalutates the regression function such as :
#' \eqn{Y = \int_0^T delta(t) X(t) dt}.
#'
#'
#' @param plsr_model a model from pls package
#' @param v_i_list a list of fd functions v_i(t) :
#' \eqn{t_i = \int_0^T v_i(t) X(t) dt}
#' @param nb_comp a value, number of components to take into account, default NULL
#' if NULL, then delta(t) will take every components available.
#'
#' @returns a fd function
#' @export
#'
#' @impotFrom fda fd
#'
#' @author Francois Bassac
evaluate_reg_curve_SPLS_uni <- function(plsr_model, v_i_list,
                                nb_comp=NULL){
  # nb_comp=length(v_i_list)

  delta = fda::fd(coef=rep(0,v_i_list[[1]]$basis$nbasis),
                  basisobj = v_i_list[[1]]$basis)

  if(is.null(nb_comp)){
    # take ALL the components >> all the v_i(t)
    nb_stop = length(v_i_list)
  }else if(nb_comp > length(v_i_list)){
    stop("evaluate_delta_PLS_V2() :
           nb_comp superior than the number of PLS steps!")
  }else{
    nb_stop = nb_comp
  }
  for(i in 1:nb_stop){
    delta = delta + plsr_model$Yloadings[i] * v_i_list[[i]]
  }

  return(delta)
}

#' build_reg_curve_spls
#'
#' This function builds the smooth PLS regression functions.
#'
#' @param plsr_model a pls model
#' @param curves_names_list a list of the names of the different curves
#' @param v_i_list a list of the v functions
#' @param nb_comp_pls_opt a integer, the number of component to take into
#' account. if null (default) all the components are used
#' @param print_steps a boolean to print steps
#'
#' @returns a list of fd object
#' @export
#'
#' @importFrom stats coef
#'
#' @author Francois Bassac
build_reg_curve_spls <- function(plsr_model, curves_names_list,
                             v_i_list, nb_comp_pls_opt = NULL,
                             print_steps = TRUE){

  if(is.null(nb_comp_pls_opt)){
    nb_comp_pls_opt = dim(plsr_model$scores)[2]
  }

  N_curves_processed = length(v_i_list)

  delta_list = list()
  for(i in 1:N_curves_processed){
    if(print_steps){
      cat(paste0("==> Build regression curve for : ",
                 curves_names_list[[i]], "\n"))
    }
    d_i = evaluate_reg_curve_SPLS_uni(plsr_model = plsr_model,
                                      v_i_list = v_i_list[[i]],
                                      nb_comp = nb_comp_pls_opt)
    delta_list = append(delta_list, list(d_i))
  }

  delta_0 = stats::coef(plsr_model,
                        ncomp = nb_comp_pls_opt,
                        intercept = TRUE)[1]


  delta_spls = list(delta_0)
  for(i in 1:length(delta_list)){
    delta_spls = append(delta_spls, list(delta_list[[i]]))
  }
  names(delta_spls) = c("Intercept", curves_names_list)
  return(delta_spls)
}

#### Prediction #####

#' smoothPLS_predict_uni
#'
#' This function make a prediction base on a dataframe and a list made of the
#' intercept and the regression curve. The input curve_type in needed to select
#' the good way of evaluate the integrals \eqn{\int_0^T delta(t) X(t) dt}.
#'
#' @param df_predict a dataframe ('id', 'time', 'state or value') to predict
#' from
#' @param delta_list a list of delta_spls : list(intercept, delta_fd)
#' @param curve_type a character, 'cat' for Categorical FD, 'num' for Scalar FD
#' @param int_mode a value of the integration mode, default 1
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#' @param regul_time a vector of time regularization values default delta_fd
#'  basis rangeval per 1, NEEDED for curve_type = 'num'!
#'
#'
#' @returns a vector of predicted values
#' \eqn{\hat{Y} = \delta_0 + \int_0^T X(t) \delta(t) dt}
#' @export
#'
#' @author Francois Bassac
smoothPLS_predict_uni <- function(df_predict, delta_list, curve_type = NULL,
                                  int_mode=1, id_col = 'id', time_col = 'time',
                                  nb_pt = 10, subdivisions = 100,
                                  regul_time = seq(
                                    delta_list[[2]]$basis$rangeval[1],
                                    delta_list[[2]]$rangeval[2],
                                    1)){

  if(!inherits(delta_list[[2]], "fd")){
    stop("smoothPLS_predict(): delta_list[[2]] have to be a
         fd object. delta_list[[1]] = delta_0")
  }

  # Time interval check
  basis_range <- delta_list[[2]]$basis$rangeval
  data_range <- range(df_predict[[time_col]])

  if(data_range[1] < basis_range[1] || data_range[2] > basis_range[2]) {
    stop(paste0("Time range mismatch: Data is in [", data_range[1], ", ",
                data_range[2],
                "], but basis is defined in [", basis_range[1], ", ",
                basis_range[2], "]."))
  }

  if(is.null(curve_type)){
    stop("smoothPLS_predict() : curve_type should be 'cat' or 'num'.")
  }else if(curve_type == 'cat'){

    y_hat = smoothPLS_CFD_predict(df_predict = df_predict,
                                  delta_spls = delta_list,
                                  int_mode = int_mode,
                                  id_col = id_col,
                                  time_col = time_col,
                                  nb_pt = nb_pt, subdivisions = subdivisions)

  }else if(curve_type == 'num'){

    if(int_mode == 1){
      cat("smoothPLS_SFD_predict() : Unstable integration with int_mode = 1,
          switch to int_mode = 2 to use pracma::trapz for integration stability.
          \n")
    }

    y_hat = smoothPLS_SFD_predict(df_predict = df_predict,
                                  delta_spls = delta_list,
                                  regul_time = regul_time,
                                  int_mode = 2,
                                  id_col = id_col, time_col = time_col,
                                  nb_pt = nb_pt, subdivisions = subdivisions)

  }else{
    stop("smoothPLS_predict() : curve_type should be 'cat' ror 'num'.")
  }
  return(y_hat)
}


#' smoothPLS_CFD_predict (Updated v0.1.2)
#'
#' @description
#' Predicts the response variable for Categorical Functional Data (CFD) by
#' integrating the regression coefficient function over active state intervals.
#'
#' @param df_predict Dataframe containing columns for id, time, and state.
#' @param delta_spls A list containing the scalar intercept and the functional
#' regression coefficient (fd object).
#' @param id_col Character, name of the id column, default 'id'.
#' @param time_col Character, name of the time column, default 'time'.
#' @param ... Additional parameters passed to evaluate_id_func_integral
#' (e.g., rel_tol, subdivisions).
#'
#' @return A numeric vector of predicted values for each individual.
#'
#' @author Francois Bassac
smoothPLS_CFD_predict <- function(df_predict, delta_spls, id_col = 'id',
                                  time_col = 'time', subdivisions = 100, ...) {

  delta_0 <- delta_spls[[1]]

  # One-time transformation of the fd object into an R function to optimize
  # performance within the loop.
  delta_1_func <- from_fd_to_func(delta_spls[[2]])

  ids_predict <- unique(df_predict[[id_col]])
  y_hat <- numeric(length(ids_predict))

  for(i in seq_along(ids_predict)) {
    df_id <- df_predict[df_predict[[id_col]] == ids_predict[i], ]

    # Call the exact segment-based integration function.
    # This respects the "Active Area" concept by only integrating over
    # intervals where the state is 1.
    res_int <- evaluate_id_func_integral(id_df = df_id,
                                         func = delta_1_func,
                                         id_col = id_col,
                                         time_col = time_col,
                                         subdivisions = 100,
                                         ...)

    y_hat[i] <- delta_0 + res_int$integral
  }

  return(y_hat)
}

#' smoothPLS_CFD_predict_Deprecated
#'
#' This function make prediction on new data from pls-cfd delta
#' Y = delta_0 + int delta_1(t) X(t) dt.
#' This function can be used for a classic functional PLS prediction.
#'
#' @param df_predict a dataframe ('id', 'time', 'state') to predict from
#' @param delta_spls a list of delta_spls
#' @param int_mode a value of the integration mode, default 1
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#'
#' @returns a vector of predicted values
#' @noRd
#'
#' @author Francois Bassac
smoothPLS_CFD_predict_Deprecated <- function(df_predict, delta_spls,
                                             int_mode = 1,
                                  id_col = 'id', time_col = 'time',
                                  nb_pt = 10, subdivisions = 100){
  # This function make prediction on new data from pls-cfd delta object
  # Y = delta_0 + int delta_1(t) X(t) dt

  delta_0 = delta_spls[[1]]
  delta_1_func = from_fd_to_func(delta_spls[[2]])

  y_hat = c()
  ids_predict = unique(df_predict[[id_col]])

  for(i in 1:length(ids_predict)){
    df_id = df_predict[df_predict[[id_col]]==ids_predict[i], ]

    y_hat_id = (delta_0 +
                  evaluate_id_func_integral(id_df=df_id,
                                            func = delta_1_func,
                                            mode = int_mode,
                                            id_col = id_col,
                                            time_col = time_col,
                                            nb_pt = nb_pt,
                                            subdivisions = subdivisions)[[2]])

    y_hat = c(y_hat, y_hat_id)
  }
  return(y_hat)
}


#' smoothPLS_SFD_predict
#'
#' @description
#' Predicts the response Y for Scalar Functional Data using the analytic L2 inner product.
#' This implementation follows the theory from Chapter 7.
#'
#' @param df_predict Dataframe with columns (id, time, value).
#' @param delta_spls List containing (intercept, delta_fd_object).
#' @param basis_obj Optional basis for signal reconstruction. If NULL, uses the basis from delta_fd
#' @param id_col Character, name of id column.
#' @param time_col Character, name of time column.
#' @param ... Additional arguments for Data2fd or inprod.
#'
#' @return A numeric vector of predicted values
#'
#' @importFrom fda Data2fd inprod
#' @author Francois Bassac
smoothPLS_SFD_predict <- function(df_predict, delta_spls, basis_obj = NULL,
                                  id_col = 'id', time_col = 'time', ...) {

  delta_0 <- delta_spls[[1]]
  delta_fd <- delta_spls[[2]]

  if(!inherits(delta_fd, "fd")) {
    stop("smoothPLS_SFD_predict() : delta_fd must be an 'fd' object.")
  }

  if(is.null(basis_obj)) {
    basis_obj <- delta_fd$basis
  }

  value_col <- setdiff(names(df_predict), c(id_col, time_col))
  ids <- unique(df_predict[[id_col]])
  y_hat <- numeric(length(ids))

  for(i in seq_along(ids)) {
    df_id <- df_predict[df_predict[[id_col]] == ids[i], ]

    # Functional representation of the test signal
    x_fd <- fda::Data2fd(argvals = df_id[[time_col]],
                         y = df_id[[value_col]],
                         basisobj = basis_obj)

    # Analytic inner product calculation <X, delta> [cite: 77, 1348]
    inner_prod <- fda::inprod(x_fd, delta_fd)

    y_hat[i] <- delta_0 + inner_prod
  }

  return(y_hat)
}



#' smoothPLS_SFD_predict_Deprecated
#'
#' This function make prediction on new data from pls-cfd delta
#' Y = delta_0 + int delta_1(t) X(t) dt.
#' This function can be used for a classic functional PLS prediction.
#'
#' @param df_predict a dataframe ('id', 'time', 'state') to predict from
#' @param delta_spls a list of delta_spls
#' @param regul_time a vector of time regularization values default delta_fd
#' basis rangeval per 1
#' @param int_mode a value of the integration mode, default 1
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#'
#' @returns a vector of predicted values
#' @noRd
#' @importFrom fda eval.fd
#' @importFrom stats approxfun integrate
#' @importFrom pracma trapz
#'
#' @author Francois Bassac
smoothPLS_SFD_predict_Deprecated <- function(df_predict, delta_spls,
                                  regul_time = seq(
                                    delta_spls[[2]]$basis$rangeval[1],
                                    delta_spls[[2]]$basis$rangeval[2], 1),
                                  int_mode=2,
                                  id_col ='id', time_col = 'time',
                                  nb_pt = 10, subdivisions = 100){

  # This function make prediction on new data from pls-cfd delta
  # Y = delta_0 + int delta_1(t) X(t) dt

  value_col = setdiff(names(df_predict), c(id_col, time_col))

  delta_0 = delta_spls[[1]]
  delta_1_func = from_fd_to_func(delta_spls[[2]])


  ids_predict = unique(df_predict[[id_col]])

  df_regul = regularize_time_series(df = df_predict, time_seq = regul_time,
                                    curve_type = 'num',
                                    id_col = id_col, time_col = time_col)

  # 1/ evaluate delta_fd on the regul_time
  eval_delta = as.vector(
    fda::eval.fd(evalarg = regul_time, fdobj = delta_spls[[2]])
  )

  y_hat = c()

  for(i in 1:length(ids_predict)){
    df_id = df_regul[df_regul[[id_col]]==ids_predict[i], ]

    prod_values = eval_delta*df_id[[value_col]]

    # Integration
    if(int_mode == 1){
      # integrate
      local_fun = stats::approxfun(x = regul_time, y = prod_values)

      y_hat_id = delta_0 + stats::integrate(f=local_fun,
                                     lower=regul_time[1],
                                     upper=regul_time[length(regul_time)],
                                     subdivisions = subdivisions)[1]$value

    } else if(int_mode == 2){
      # pracma::trapz
      y_hat_id = delta_0 + pracma::trapz(x = regul_time, y = prod_values)

    }

    y_hat = c(y_hat, y_hat_id)
  }
  return(y_hat)
}

#' smoothPLS_multi_predict
#'
#' This function use the list of regression functions to make a prediction
#' \eqn{\hat{Y} = \delta_0 + \int_0^T \delta(t) X(t) dt}
#'
#' @param df_predict_list a list of dataframe (id, time, value_or_state)
#' @param delta_list a list of regression object (intercept, delta_1_fd,
#' delta_2_fd, etc)
#' @param curve_type_obj a list of characters of the curve types 'cat' or 'num'
#' @param id_col_obj a list of character of the name of the id column,
#' default 'id'
#' @param time_col_obj a list of character of the name of the time column,
#' default 'time'
#' @param regul_time_obj a list of the time regularization values
#' @param int_mode a integer for the integration mode, 1 for integrate, 2 for
#' pracma::trapz
#' @param nb_pt a integrer, number of intermediate points for pracma::trapz,
#' default 10
#' @param subdivisions a integer, number of subdivision in integrate function,
#' default 100
#'
#' @returns a numeric vector of the prediction
#' @export
#'
#' @author Francois Bassac
smoothPLS_predict <- function(df_predict_list, delta_list,
                              curve_type_obj = NULL,
                              id_col_obj = 'id', time_col_obj = 'time',
                              regul_time_obj = NULL, int_mode = 1,
                              nb_pt = 10, subdivisions = 100){

  # assertions
  N_curves_processed = length(delta_list) - 1 # minus the intercept

  # N_curves
  if(length(df_predict_list) != 3){
    N_curves = length(df_predict_list)
  }else if(length(df_predict_list) == 3){
    # if 3, either 1 curve or 3
    if(mode(df_predict_list[[1]]) == "numeric" ||
       mode(df_predict_list[[1]]) == "character" ){
      N_curves = 1
    }else{
      N_curves = 3
    }
  }

  if(mode(df_predict_list) != "list" && N_curves != 1){
    stop("smoothPLS_multi_predict() : df_predict_list have to be a list of
         dataframe (id, time, value).")
  }else if(mode(df_predict_list) == "list" && N_curves != 1 &&
           !all(sapply(df_predict_list, function(x) ncol(x)==3 ))){
    stop("smoothPLS_multi_predict() : all df of df_predict_list should have
         3 columns, (id, time, value).")
  }else if(N_curves == 1 && ncol(df_predict_list) != 3){
    stop("smoothPLS_multi_predict() : df_predict_list have to be a list of
         dataframe (id, time, value).")
  }


  # curve_type_obj
  if(is.null(curve_type_obj)){
    stop("smoothPLS_multi_predict() : curve_type_obj have to be either :
    'cat' or 'num' to be applied to all the curves OR
     a list of ('cat', 'num', ....) to be used for each curve.")
  } else if(all(curve_type_obj == 'cat')){
    curve_type_list = obj_list_creation(N_rep = N_curves, obj = curve_type_obj)
  } else if(all(curve_type_obj == 'num')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(!(mode(curve_type_obj) == 'list' &&
              all(sapply(curve_type_obj, function(x) (x=='cat' | x=='num'))))){
    stop("smoothPLS_multi_predict() : curve_type_obj have to be either \'cat\'
         \'num\', for all curves, or a list  or a list specify 'cat' or 'num'
         or for each curve.")
  } else {
    curve_type_list = curve_type_obj
  }

  # regul_time_list
  if(is.null(regul_time_obj)  &&
     (curve_type_obj == 'cat' || curve_type_obj[[1]] == 'cat')){
    regul_time_obj = c(
      delta_list[[2]]$basis$rangeval[1]:delta_list[[2]]$basis$rangeval[2]
    )
  }
  if(mode(regul_time_obj) == "numeric"){
    # Because the functions basis_list_creation works :
    regul_time_list =  obj_list_creation(N_rep = N_curves, regul_time_obj)
  } else if(mode(regul_time_obj) == "list" &&
            all(sapply(regul_time_obj, function(x) inherits(x, "numeric")))){
    regul_time_list = regul_time_obj
  } else {
    stop("smoothPLS_multi_predict() : regul_time_obj have to be a numeric
    sequence or a list of numeric sequence.")
  }

  # id_col_obj & time_col_obj
  if(mode(id_col_obj) == "character"){
    # Because the functions basis_list_creation works :
    id_col_list = obj_list_creation(N_rep = N_curves, id_col_obj)
  } else if(mode(id_col_obj) == "list" &&
            all(sapply(id_col_obj, function(x) inherits(x, "character")))){
    id_col_list = id_col_obj
  } else {
    stop("smoothPLS_multi_predict() : id_col_obj have to be a character or
         a list of characters.")
  }

  if(mode(time_col_obj) == "character"){
    # Because the functions basis_list_creation works :
    time_col_list = obj_list_creation(N_rep = N_curves, time_col_obj)
  } else if(mode(time_col_obj) == "list" &&
            all(sapply(time_col_obj, function(x) inherits(x, "character")))){
    time_col_list = time_col_obj
  } else {
    stop("smoothPLS_multi_predict() : time_col_obj have to be a character or
         a list of characters.")
  }

  if(N_curves == 1){
    if((mode(df_predict_list[[1]]) == "numeric" ||
        mode(df_predict_list[[1]]) == "character") &&
       length(df_predict_list) !=3){
      stop("smoothPLS_multi_predict() : df_predict_list shoud have 3 columns :
            id, time, value_state.")
    }
  }
  if(mode(df_predict_list) != "list" && N_curves != 1){
    stop("smoothPLS_multi_predict() : df_predict_list have to be a list of
         dataframe (id, time, value).")
  }else if(mode(df_predict_list) == "list" && N_curves != 1 &&
           !all(sapply(df_predict_list, function(x) ncol(x)==3 ))){
    stop("smoothPLS_multi_predict() : all df of df_predict_list should have
         3 columns, (id, time, value).")
  }else if(N_curves == 1 && ncol(df_predict_list) != 3){
    stop("smoothPLS_multi_predict() : df_predict_list have to be a list of
         dataframe (id, time, value).")
  }

  if(N_curves == 1 &&
     mode(df_predict_list[[1]]) != 'list' &&
     ncol(df_predict_list) == 3){
    df_predict_list = list(df_predict_list)
  }

  # Data preparation
  new_list_obj = build_new_data_list(df_list = df_predict_list,
                                     N_curves = N_curves,
                                     orth_basis_list = NULL,
                                     basis_list = NULL,
                                     curve_type_list = curve_type_list,
                                     id_col_list = id_col_list,
                                     time_col_list = time_col_list,
                                     regul_time_list = regul_time_list)

  df_processed_list = new_list_obj$df_processed_list
  curves_names_list = new_list_obj$curves_names_list
  new_curves_type_list = new_list_obj$new_curves_type_list
  new_id_col_list = new_list_obj$new_id_col_list
  new_time_col_list = new_list_obj$new_time_col_list
  new_regul_time_list = new_list_obj$new_regul_time_list

  Y_predict = delta_list[[1]] # The intercept
  # for all N_curves
  for(i in 1:length(df_processed_list)){

    # print(paste0("Prediction for : ", curves_names_list[[i]]))

    # temp_df
    temp_df = df_processed_list[[i]]

    # temp_curve_type
    temp_curve_type = new_curves_type_list[[i]]

    # temp_delta_list
    temp_delta_list = list(0, delta_list[[i+1]])

    # temp_regul_time
    tem_regul_time = new_regul_time_list[[i]]

    temp_Y = smoothPLS_predict_uni(df_predict = temp_df,
                                   delta_list = temp_delta_list,
                                   curve_type = temp_curve_type,
                                   int_mode = int_mode,
                                   id_col = new_id_col_list[[i]],
                                   time_col = new_time_col_list[[i]],
                                   nb_pt = nb_pt,
                                   subdivisions = subdivisions,
                                   regul_time = tem_regul_time
    )

    Y_predict = Y_predict + temp_Y

  }
  return(Y_predict)
}


#' help_smoothPLS
#'
#' This function print some information on the package use.
#'
#' @export
#'
#' @author Francois Bassac
help_smoothPLS <- function(){
  cat("Here are some informations to use the SmoothPLS package :
      \n")
  cat("=> Performs univariate and multivariate PLS methods on functional data. \n")
  cat("==> The functional data can be scalar or categorial. \n")
  cat("==> 3 methods are implemented : smooth PLS, functional PLS and naive PLS.\n")
  cat("==> The functions are smoothPLS, funcPLS and naivePLS.
      \n")
  cat("=> The data format have to be a list of dataframe of 3 columns. \n")
  cat("==> The columns should reprensent the id, the time and the value_or_state. \n")
  cat("==> Per dataframe, all id must have the SAME start time and end time.
      \n")
  cat("=> In the functions, be carefull when you fill the 'curve_type_list' parameter.\n")
  cat("==> The 'curve_type_obj' parameter is a list with 'cat' or 'num'.\n")
  cat("==> For the arguments finishing in _obj, if you put only one value when
      you have more than 1 curve, the programm will build the list itself.
      \n")
}

