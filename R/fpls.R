
#### Multivariate ####

#' funcPLS
#'
#' This function performs the Multivariate Functional PLS as a matrix problem.
#'
#' @param df_list a list of dataframes (id, time, value_or_state)
#' @param Y a numeric vector of the response
#' @param basis_obj a basis fd obj or a list of basis fd obj. If basis fd obj,
#' the same basis is used for all the curves
#' @param regul_time_obj a vector of time regularization values or a list of vectors
#' @param curve_type_obj a character "cat" or 'num' or a list of those values
#' @param id_col_obj a character of the id column for all the curves or a list of id column character
#' @param time_col_obj a character of the time column for all the curves or a list of time column character
#' @param print_steps a boolean to cat the current step
#' @param plot_rmsep a boolean to plot the plsr RMSEP
#' @param print_nbComp a boolean to cat the optimal number of components
#' @param plot_reg_curves a boolean to directly plot the beta regression curves
#' @param jackknife a plsr input, defautl = TRUE
#' @param validation a plsr input, default = 'LOO'
#'
#' @returns a list ("curve_names", "alphas", "metric", "root_metric",
#' "trans_alphas", "mfpls_mfd", "nb_comp_pls_opt", "beta_0", "beta_pls_list")
#' @export
#'
#' @import pls
#' @importFrom mgcv mroot
#' @importFrom pls plsr RMSEP MSEP
#' @importFrom MASS ginv
#' @importFrom stats coef
#' @importFrom graphics plot title
#' @importFrom fda fd
#'
#' @author Francois Bassac
funcPLS <- function(df_list, Y, basis_obj, regul_time_obj,
                    curve_type_obj = NULL,
                    id_col_obj = 'id', time_col_obj = 'time',
                    print_steps = FALSE, plot_rmsep = TRUE,
                    print_nbComp = TRUE, plot_reg_curves = FALSE,
                    jackknife = TRUE, validation = 'LOO'){
  #multivariate_fpls
  if(print_nbComp){
    cat("### Functional PLS ### \n")
  }
  # 1 assertion

  if(print_steps){
    cat("=> Input format assertions.\n")
  }

  assert_objs = assert_funcPLS_inputs(df_list = df_list, Y = Y,
                                      basis_obj = basis_obj,
                                      regul_time_obj =  regul_time_obj,
                                      curve_type_obj = curve_type_obj,
                                      id_col_obj = id_col_obj,
                                      time_col_obj = time_col_obj)

  N_curves = assert_objs$N_curves
  basis_list = assert_objs$basis_list
  regul_time_list = assert_objs$regul_time_list
  curve_type_list = assert_objs$curve_type_list
  id_col_list = assert_objs$id_col_list
  time_col_list = assert_objs$time_col_list

  if(print_steps){
    cat("=> Input format assertions OK.\n")
  }

  if(N_curves == 1 && mode(df_list[[time_col_obj[[1]]]]) == 'numeric' &&
     ncol(df_list) == 3){
    df_list = list(df_list)
  }


  if(print_steps){
    cat("=> Building alpha matrix.\n")
    cat("=> Building curve names.\n")
  }
  # Alphas matrix and curve_names
  alpha_names_obj = multivariate_alpha_building(df_list = df_list,
                                                basis_list = basis_list,
                                                curve_type_list =curve_type_list,
                                                regul_time_list =regul_time_list,
                                                id_col_list = id_col_list,
                                                time_col_list = time_col_list,
                                                print_steps = print_steps)
  alphas = alpha_names_obj[[1]]
  curve_names = alpha_names_obj[[2]]
  new_basis_list = alpha_names_obj[[3]]

  # Metric
  if(print_steps){
    cat("=> Evaluate metrix and root_metric.\n")
  }
  # WARNING on MULTISTATES!!

  metric = multivariate_assemble_basis_metric(df_list = df_list,
                                              basis_list =  basis_list,
                                              curve_list = curve_type_list,
                                              id_col_list = id_col_list,
                                              time_col_list = time_col_list)

  root_metric = mgcv::mroot(metric, method="svd")

  # trans_alphas
  trans_alphas = alphas%*%root_metric #use the metric

  if(print_steps){
    cat("=> plsr(Y ~ trans_alphas).\n")
  }
  # FPLS
  mfpls_mfd = pls::plsr(Y ~ trans_alphas,
                        validation="LOO", jackknife = TRUE, center = TRUE,
                        intercept = TRUE)

  if(plot_rmsep){
    plot(pls::RMSEP(mfpls_mfd), main="Cross Validation",
         xlab="Components number")
  }

  nb_comp_pls_opt = which.min(pls::RMSEP(mfpls_mfd)$val[1, , ])-1
  if(print_nbComp){
    cat(paste("Optimal number of PLS components : ", nb_comp_pls_opt, ".\n"))
  }

  # Assert a non equal to zero optimal number of component

  if(nb_comp_pls_opt == 0){
    nb_total_cp = dim(trans_alphas)[2]
    cat("No Optimal number of component!\n")
    cat("Search the optimal number of component EXCLUDING the intercept.\n")
    nb_comp_pls_opt = which.min(pls::RMSEP(mfpls_mfd)$val[1,
                                                          , ][c(2:nb_total_cp)])
    cat(paste("New \'Optimal\' number of PLS components : ",
              nb_comp_pls_opt, ".\n"))
  }

  if(print_steps){
    cat("=> Build Intercept and regression curves for optimal number of components.\n")
  }
  # beta_t_i and beta_0_i
  beta_pls_list = reg_curve_funcPLS_evaluation(mfpls_mfd = mfpls_mfd,
                                               nb_comp_pls_opt =nb_comp_pls_opt,
                                               root_metric = root_metric,
                                               new_basis_list = new_basis_list,
                                               curve_names =  curve_names,
                                               print_steps = print_steps)

  if(plot_reg_curves){
    for(i in 1:length(beta_pls_list)){
      plot(beta_pls_list[[i]])
      title(curve_names[[i]])
    }
  }

  beta_0 = stats::coef(mfpls_mfd, ncomp = nb_comp_pls_opt, intercept = TRUE)[1]

  full_beta_list = append(list(beta_0), beta_pls_list)
  names(full_beta_list) = c("Intercept", curve_names)


  # Extracting coefficients with MASS::ginv for the inverse of the metric
  #coefs_alphas = stats::coef(mfpls_mfd,
  #                           ncomp = nb_comp_pls_opt,
  #                           intercept=FALSE)
  #betacoef = MASS::ginv(t(root_metric)) %*% coefs_alphas


  # return also the df_fd$curve_name?
  # fpls_mv_list = list(curve_names, alphas, metric, root_metric,
  #                    trans_alphas, mfpls_mfd, nb_comp_pls_opt,
  #                    full_beta_list)
  fpls_mv_list = list(metric, root_metric, alphas, trans_alphas,
                      mfpls_mfd, nb_comp_pls_opt, curve_names,
                      full_beta_list)
  names(fpls_mv_list) = c("metric", "root_metric", "alphas", "trans_alphas",
                          "plsr_model", "nbCP_opti", "curves_names", "reg_obj")

  return(fpls_mv_list)

}


#' assert_funcPLS_inputs
#'
#' This function checks the integrity of the input for funcPLS.
#' It returns a list of (basis_list, regul_time_list, curve_type_list,
#' id_col_list, time_col_list)
#'
#' @param df_list a list of dataframes (id, time, value_or_state)
#' @param Y a numeric vector of the response
#' @param basis_obj a list of basis object or a basis object
#' @param regul_time_obj a vector of time regularization values or a list of vectors
#' @param curve_type_obj a character "cat" or 'num' or a list of those values
#' @param id_col_obj a character of the id column for all the curves or a list of id column character
#' @param time_col_obj a character of the time column for all the curves or a list of time column character
#'
#' @returns a list of (basis_list, regul_time_list, curve_type_list,
#' id_col_list, time_col_list)
#'
#' @author Francois Bassac
assert_funcPLS_inputs <- function(df_list, Y, basis_obj,
                                            regul_time_obj,
                                            curve_type_obj = NULL,
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
    stop("smoothPLS_multi() : df_list have to be a list of dataframe (id, time, value).")
  }else if(mode(df_list) == "list" && N_curves != 1 &&
           !all(sapply(df_list, function(x) ncol(x)==3 ))){
    stop("smoothPLS_multi() : all df of df_list should have 3 columns, (id, time, value).")
  }else if(N_curves == 1 && mode(df_list[[1]]) != 'list' && ncol(df_list) != 3){
    stop("naivePLS() : df_list have to be a list of dataframe
         (id, time, value).")
  }

  # Y
  if((!(is.vector(Y) && mode(Y)=="numeric") &
      !(is.list(Y) && all(sapply(Y, function(x) inherits(x, "numeric")))))){
    stop("funcPLS() : Y should be a numeric vector.")
  }

  # basis_list
  if(inherits(basis_obj, "basisfd")){
    basis_list = obj_list_creation(N_rep = N_curves, obj = basis_obj)
  }else if(mode(basis_obj) == "list" &&
           all(sapply(basis_obj, function(x) inherits(x, "basisfd")))){
    basis_list = basis_obj
  } else {
    stop("funcPLS() : basis_obj have to be a basis fd object or a list
          of basis fd objects.")
  }

  # regul_time_list
  if(mode(regul_time_obj) == "numeric"){
    regul_time_list =  obj_list_creation(N_rep = N_curves, regul_time_obj)
  } else if(mode(regul_time_obj) == "list" &&
            all(sapply(regul_time_obj, function(x) inherits(x, "numeric")))){
    regul_time_list = regul_time_obj
  } else {
    stop("funcPLS() : regul_time_obj have to be a numeric sequence or
         a list of numeric sequence.")
  }

  if(is.null(curve_type_obj)){
    stop("funcPLS() : curve_type_obj missing and have to be either :
    'cat' or 'num' to be applied to all the curves OR
     a list of ('cat', 'num', ....) to be used for each curve.")
  } else if(all(curve_type_obj == 'cat')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(all(curve_type_obj == 'num')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(!(mode(curve_type_obj) == 'list' &&
             all(sapply(curve_type_obj, function(x) (x=='cat' | x=='num'))))){
    stop("funcPLS() : curve_type_obj have to be either \'cat\' or
         \'num\', for all curves, or a list  or a list specify 'cat' or 'num'
         for each curve.")
  } else {
    curve_type_list = curve_type_obj
  }

  # id_col_obj & time_col_obj
  if(mode(id_col_obj) == "character"){
    id_col_list = obj_list_creation(N_rep = N_curves, id_col_obj)
  } else if(mode(id_col_obj) == "list" &&
            all(sapply(id_col_obj, function(x) inherits(x, "character")))){
    id_col_list = id_col_obj
  } else {
    stop("funcPLS() : id_col_obj have to be a character or
         a list of characters.")
  }

  if(mode(time_col_obj) == "character"){
    time_col_list = obj_list_creation(N_rep = N_curves, time_col_obj)
  } else if(mode(time_col_obj) == "list" &&
            all(sapply(time_col_obj, function(x) inherits(x, "character")))){
    time_col_list = time_col_obj
  } else {
    stop("funcPLS() : time_col_obj have to be a character or
         a list of characters.")
  }

  # last check : number!
  if(length(basis_list) != N_curves){
    stop("funcPLS() : basis_obj should be only one basis or a list of
         length the number of curves.")
  }

  if(length(regul_time_list) != N_curves){
    stop("funcPLS() : regul_time_obj should be only one regul_time
         vector or a list of length the number of curves.")
  }

  if(length(curve_type_list) != N_curves){
    stop("funcPLS() : curve_type_obj should be only one curve_type
         'cat' or 'num', or a list of length the number of curves.")
  }

  if(length(id_col_list) != N_curves){
    stop("funcPLS() : id_col_obj should be only one character,
         or a list of length the number of curves.")
  }

  if(length(time_col_list) != N_curves){
    stop("funcPLS() : time_col_obj should be only one character,
         or a list of length the number of curves.")
  }

  assert_objs = list(N_curves, basis_list, regul_time_list, curve_type_list,
                     id_col_list, time_col_list)
  names(assert_objs) = c("N_curves", "basis_list", "regul_time_list", "curve_type_list",
                        "id_col_list", "time_col_list")

  return(assert_objs)
}

#' multivariate_alpha_building
#'
#' This function builds the alphas matrix by columns binding the alpha of each
#' curve.
#' It returne the alphas matrix and a vector containing all the curves names.
#'
#' @param df_list a list of data frame (id, time, state_or_value)
#' @param basis_list a list of basis
#' @param curve_type_list a list of curves
#' @param regul_time_list a list of regul_time
#' @param id_col_list a list if the characters for the id column
#' @param time_col_list a list of the characters for the time column
#' @param print_steps a boolean to print the current step, default FALSE
#'
#' @returns a list of the alphas matrix and the curve_names vector
#' and the new_basis_list
#'
#' @author Francois Bassac
multivariate_alpha_building <- function(df_list, basis_list,
                                        curve_type_list, regul_time_list,
                                        id_col_list, time_col_list,
                                        print_steps = FALSE){

  N_curves = length(df_list)

  alphas = c()
  curve_names = c()
  new_basis_list = list()

  for(i in 1:N_curves){

    if(curve_type_list[[i]] =='cat'){
      state_col = setdiff(names(df_list[[i]]),
                        c(id_col_list[[i]], time_col_list[[i]]))
      states = unique(df_list[[i]][[state_col]])
      states_names = states[order(states)]

      if(all(states_names %in% c(0, 1)) ||
         all(states_names %in% c('0', '1')) ){
        # case : one state CFD in indicatrice form
        states_names = paste0(state_col, "_1")
      }

      curve_name_loc = determine_curve_name(curve_number = i,
                                            curve_type = curve_type_list[[i]],
                                            states_names = states_names)

      basis_list_local = obj_list_creation(N_rep = length(states_names),
                                           obj = basis_list[[i]])
      new_basis_list = append(new_basis_list, basis_list_local)

    } else if(curve_type_list[[i]] =='num'){

      curve_name_loc = determine_curve_name(curve_number = i,
                                            curve_type = curve_type_list[[i]])

      new_basis_list = append(new_basis_list, list(basis_list[[i]]))
    }

    if(print_steps){
      cat(paste0("==> Evaluating alpha for : ", curve_name_loc, ".\n"))
    }


    alphas_loc = univariate_alpha_building(df = df_list[[i]],
                                           basis = basis_list[[i]],
                                           curve_type = curve_type_list[[i]],
                                           regul_time = regul_time_list[[i]],
                                           id_col = id_col_list[[i]],
                                           time_col = time_col_list[[i]])


    alphas = cbind(alphas, alphas_loc)
    curve_names = c(curve_names, curve_name_loc)

  }
  return(list(alphas, curve_names, new_basis_list))
}

#' univariate_alpha_building
#'
#' This function build the alphas matrix for a dataframe of a curve_type curve
#'  on a basis on a regul_time vector.
#'
#' @param df a dataframe (id, time, value_or_state)
#' @param basis a basis fd object
#' @param curve_type a character, 'cat' or 'num'
#' @param regul_time a numeric vector for time regularization
#' @param id_col a character for the id column name
#' @param time_col a character for the time column name
#'
#' @returns a matrix
#'
#' @importFrom fda Data2fd
#' @importFrom stats coef
#'
#' @author Francois Bassac
univariate_alpha_building <- function(df, basis, curve_type = NULL, regul_time,
                                      id_col = 'id', time_col = 'time'){

  if(is.null(curve_type)){
    stop("univariate_alpha_building() : curve_type should be 'cat' ror 'num'.")
  }else if(curve_type == 'num'){

    df_new = regularize_time_series(df, time_seq = regul_time,
                                    curve_type = curve_type,
                                    id_col = id_col, time_col = time_col)

    df_new_wide = convert_to_wide_format(df_new, id_col = id_col,
                                         time_col = time_col)

    # Data2fd : compute the coefs
    df_fd = fda::Data2fd(argvals = regul_time,
                         y = t(df_new_wide[, -c(1)]), basis)

    alphas = t(stats::coef(df_fd))

  }else if(curve_type == 'cat'){
    alphas = c()

    state_col = setdiff(names(df), c(id_col, time_col))
    states = unique(df[state_col])

    if(all(states[[state_col]] %in% c(0, 1)) ||
       all(states[[state_col]] %in% c('0', '1')) ){
      # case : one state CFD in indicatrice form

      df_processed = list(df)
      names(df_processed) = paste0(state_col, "_1")
      N_states = 1

    } else {
      # case : multistates CFD
      N_states = length(unique(df[[state_col]]))

      df_processed = cat_data_to_indicatrice(df,
                                             id_col = id_col,
                                             time_col = time_col)
    }

    for(i in 1:N_states){
      #print(names(df_processed)[i])

      df_temp = df_processed[[i]]

      df_new = regularize_time_series(df = df_temp, time_seq = regul_time,
                                      curve_type = curve_type,
                                      id_col = id_col, time_col = time_col)

      df_new_wide = convert_to_wide_format(df_new,
                                           id_col = id_col,
                                           time_col = time_col)

      # Data2fd : compute the coefs
      df_fd = fda::Data2fd(argvals = regul_time,
                           y = t(df_new_wide[, -c(1)]),
                           basis)

      #plot(df_fd)
      #title(names(df_processed)[i])

      # alpha
      alphas_loc = t(stats::coef(df_fd))

      alphas = cbind(alphas, alphas_loc)

    }

    # return aussi la liste ordonnée des états,
    # return aussi le nom de la courbe, SFD_1, CFD_2_state_1, CFD_2_state_2, etc
  }

  return(alphas)
}

#' determine_curve_name
#'
#' This function determines the curves names bases on its place in the data and
#' its states if it is a categorical functional data.
#'
#' @param curve_number a int, the place of the vurve in the data
#' @param curve_type a character, 'cat' or 'num'
#' @param states_names a character or list of character with the states of the CFD
#'
#' @returns a vector of names.
#'
#' @author Francois Bassac
determine_curve_name <- function(curve_number, curve_type=NULL,
                                 states_names=NULL){
  # return the name of the curve
  # is a character vector

  if(!(curve_number > 0 & curve_number %% 1 == 0)){
    stop("determine_curve_name() : curve_number should be an integer > 0.")
  }

  if(is.null(curve_type) || (curve_type != 'cat' & curve_type != 'num')){
    stop("determine_curve_name() : curve_type should be 'cat' or 'num'.")
  } else if(curve_type == 'num'){
    curve_name = c(paste0("NumFD_", curve_number))
  } else if(curve_type == 'cat'){
    if(is.null(states_names) |
       !(is.vector(states_names) &&
        (all(sapply(states_names, function(x) inherits(x, "character"))) ||
         all(sapply(states_names, function(x) inherits(x, "numeric")))))
       ){
          stop("determine_curve_name() : states_names should be a character or
               a vector of characters or numeric only.")
    }else{
            states_names_ordered = states_names[order(states_names)]
            curve_name = c(paste0("CatFD_", curve_number,
                                  "_", states_names_ordered))
          }
  }
  return(curve_name)
}

#' multivariate_assemble_basis_metric
#'
#' This function assemble the metrics of all the basis of the basis list for
#' multivariate use case.
#' It take care on the categorical case by de-multiplying the basis per state.
#'
#' @param basis_list a list of basis fd object
#' @param df_list a list of dataframe (id, time, value_or_state)
#' @param curve_list a list of the curves type 'cat' or 'num'
#' @param id_col_list a list of the character of the id columns
#' @param time_col_list a list of the character of the time columns
#'
#' @returns a matrix of the metric to consider
#' @export
#'
#'
#' @author Francois Bassac
multivariate_assemble_basis_metric <- function(df_list, basis_list, curve_list,
                                               id_col_list, time_col_list){
  # This function assemble the metrics of all the basis of the basis list.

  N_curves = length(df_list)

  for(i in 1:N_curves){
    if(curve_list[[i]] == 'num'){
      if(i==1){
        metric = evaluate_metric(basis = basis_list[[i]])
      }else{
        metric0 = evaluate_metric(basis = basis_list[[i]])
        metric <- block_diag(metric, metric0)
      }
    } else if(curve_list[[i]] == 'cat'){

      state_col = setdiff(names(df_list[[i]]),
                          c(id_col_list[[i]], time_col_list[[i]]))

      states = unique(df_list[[i]][[state_col]])
      states_names = states[order(states)]

      if(all(states_names %in% c(0, 1)) ||
         all(states_names %in% c('0', '1')) ){
        # case : one state CFD in indicatrice form
        N_states = 1

      } else{
        N_states = length(unique(df_list[[i]][[state_col]]))

      }

      # list of basis for all the states of the i^th 'cat' curve
      list_basis = obj_list_creation(N_rep = N_states,
                                       obj = basis_list[[i]])
      if(i==1){
        metric = assemble_basis_metric(basis_list = list_basis)
      }else{
        metric0 = assemble_basis_metric(basis_list = list_basis)
        metric <- block_diag(metric, metric0)
      }
    }
  }
  return(metric)
}

#' reg_curve_funcPLS_evaluation
#'
#' This function evaluate the beta regression functions of the multivariate
#' functional PLS.
#'
#' @param mfpls_mfd a MFPLS model
#' @param nb_comp_pls_opt the optimal number of components
#' @param root_metric the root metric
#' @param curve_names a list of curves to keep
#' @param new_basis_list a list of basis fd objects
#' @param print_steps a boolean to cat or not the steps
#'
#' @returns a list of the betas per state
#' @export
#'
#' @importFrom fda fd
#' @importFrom MASS ginv
#' @importFrom stats coef
#'
#' @author Francois Bassac
reg_curve_funcPLS_evaluation <- function(mfpls_mfd, nb_comp_pls_opt,
                                              root_metric, new_basis_list,
                                              curve_names,
                                              print_steps = FALSE){
  # This function evaluates the beta regression curves for FPLS for multivariate
  # FD.
  # It returns a list of beta.

  coefs_alphas = stats::coef(mfpls_mfd,
                             ncomp = nb_comp_pls_opt,
                             intercept=FALSE)
  betacoef = MASS::ginv(t(root_metric))%*%coefs_alphas

  beta_mfpls = list()
  start_indice = 1
  end_indice = new_basis_list[[1]]$nbasis
  for(i in 1:length(curve_names)){
    name = curve_names[i]
    if(print_steps){
      cat(paste0(" ==> Build  regression curve for : ", name, "\n"))
    }
    beta_mfpls[[name]] = fda::fd(betacoef[start_indice:end_indice],
                                 new_basis_list[[i]])
    start_indice = start_indice + new_basis_list[[i]]$nbasis
    if(i == length(curve_names)){
      break
    } else {
      end_indice = end_indice + new_basis_list[[i+1]]$nbasis
    }
  }

  return(beta_mfpls)
}


#' funcPLS_predict
#'
#' This function performs the prediction on a df_predict_ms using the delta_list
#' \eqn{\hat(Y) = \delta_0 + \sum_{i=1}^K \int_0^T \delta_i(t) X_i(t) dt}.
#'
#'
#' @param df_predict_list a list of dataframe (id, time, state_or_value)
#' @param delta_list a list of regression objects (Intercept, fd, etc)
#' @param curve_type_obj a list of the curves types 'cat' or 'num'
#' @param regul_time_obj a list of time regularization values
#' @param id_col_obj a list of characters of the names of the id columns
#' @param time_col_obj a list of characters of the names of the time columns
#' @param int_mode a integer for integration mode, 1 for integrate,
#' 2 for pracma::trapz, default 1
#' @param nb_pt a integer, the number of intermediate points for integration
#' mode 2, default 10
#' @param subdivisions a integer, the number of subdivisions for integration
#' mode 1, default 100
#'
#' @returns a numeric vector
#' @export
#'
#' @author Francois Bassac
funcPLS_predict <- function(df_predict_list,
                            delta_list,
                            curve_type_obj = NULL,
                            regul_time_obj = NULL,
                            id_col_obj = 'id',
                            time_col_obj = 'time',
                            int_mode = 1,
                            nb_pt = 10, subdivisions = 100){

  Y_hat = smoothPLS_predict(df_predict_list = df_predict_list,
                            delta_list = delta_list,
                            curve_type_obj = curve_type_obj,
                            regul_time_obj = regul_time_obj,
                            id_col_obj = id_col_obj,
                            time_col_obj = time_col_obj,
                            int_mode = int_mode,
                            nb_pt = nb_pt,
                            subdivisions = subdivisions)
  return(Y_hat)
}

