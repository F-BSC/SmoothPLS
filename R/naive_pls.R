# Contain pls non functional (naive). the idea is to be able to compare the approaches.

#' naivePLS
#'
#' This function performs the naive PLS method for Categorical functional data,
#' Scalar functional data and multivariate data.
#'
#' @param df_list a list of dataframe (id, time, value_or_state)
#' @param Y a numeric vector for the scalar response
#' @param regul_time_obj a list of time regularization values
#' @param curve_type_obj a list of the curve types 'cat' or 'num'
#' @param id_col_obj a list of character of the names of the id columns
#' @param time_col_obj a list of character of the names of the time columns
#' @param print_steps a boolean to print the different steps, default FALSE
#' @param plot_rmsep a boolean to plot the RMSEP, default TRUE
#' @param print_nbComp a boolean to print the optimal number or components,
#' default TRUE
#' @param plot_reg_curves a boolean to plot the regression curves, default
#' FALSE
#' @param validation a character, pls::plsr input, default 'LOO'
#' @param jackknife a boolean, pls::plsr input, default TRUE
#'
#' @returns a list of ("plsr_model", "nbCP_opti", "curves_names",
#' "opti_reg_coef", "reg_obj")
#' @export
#'
#' @import pls
#' @importFrom pls plsr RMSEP MSEP
#' @importFrom stats coef predict residuals lm.influence
#'
#' @author Francois Bassac
naivePLS <- function(df_list, Y,
                     regul_time_obj = NULL, curve_type_obj = NULL,
                     id_col_obj = 'id', time_col_obj = 'time',
                     print_steps = FALSE,
                     plot_rmsep = TRUE, print_nbComp = TRUE,
                     plot_reg_curves = FALSE,
                     validation ='LOO', jackknife = TRUE){
  # performs the naive pls for all cases.
  cat("### Naive PLS ### \n")

  # N_curves
  if(print_steps){
    cat("=> Input format assertions.\n")
  }

  assert_obj = assert_multivariate_naivePLS_inputs(
    df_list = df_list, Y = Y,
    regul_time_obj = regul_time_obj,
    curve_type_obj = curve_type_obj,
    id_col_obj = id_col_obj,
    time_col_obj = time_col_obj)

  N_curves = assert_obj$N_curves
  regul_time_list = assert_obj$regul_time_list
  curve_type_list = assert_obj$curve_type_list
  id_col_list = assert_obj$id_col_list
  time_col_list = assert_obj$time_col_list

  if(print_steps){
    cat("=> Input format assertions OK.\n")
  }

  if(N_curves == 1 && mode(df_list[[1]]) != 'list' && ncol(df_list) == 3){
    df_list = list(df_list)
  }

  curves_names = c()
  curve_type_updated = c()
  regul_time_updated = list()

  if(print_steps){
    cat("=> Data formatting.\n")
  }

  for(i in 1:N_curves){

    temp_df = df_list[[i]]
    # For all curves
    if(curve_type_list[[i]] == 'cat'){

      state_col = setdiff(names(df_list[[i]]),
                          c(id_col_list[[i]], time_col_list[[i]]))
      states = unique(df_list[[i]][state_col])

      if(all(states[[state_col]] %in% c(0, 1)) ||
         all(states[[state_col]] %in% c('0', '1')) ){
        # case : one state CFD in indicator form

        df_processed = list(df_list[[i]])
        names(df_processed) = paste0(state_col, "_1")
      } else {
        # case : multistates CFD
        df_processed = cat_data_to_indicator(df_list[[i]],
                                               id_col = id_col_list[[i]],
                                               time_col = time_col_list[[i]])
      }

      for(j in 1:length(df_processed)){

        temp_df_new = regularize_time_series(df_processed[[j]],
                                             time_seq = regul_time_list[[i]],
                                             curve_type = 'cat')
        temp_df_new_wide = convert_to_wide_format(temp_df_new,
                                                  id_col = id_col_list[[i]],
                                                  time_col = time_col_list[[i]])
        if(i == 1 && j == 1){
          df_mod_wide = temp_df_new_wide[, -c(1)]
        } else {
          df_mod_wide = cbind(df_mod_wide, temp_df_new_wide[, -c(1)])
        }
        curves_names = c(curves_names, paste0("Cat_", i,
                                              "_", names(df_processed[j])))

        curve_type_updated = c(curve_type_updated, 'cat')

        regul_time_updated = append(regul_time_updated,
                                    list(regul_time_list[[i]]))

      }

    }else if(curve_type_list[[i]] == 'num'){
      # Data manipulation
      temp_df_new = regularize_time_series(temp_df,
                                           time_seq = regul_time_list[[i]],
                                           curve_type = 'num')

      temp_df_new_wide = convert_to_wide_format(temp_df_new,
                                                id_col = id_col_list[[i]],
                                                time_col = time_col_list[[i]])
      if(i == 1){
        df_mod_wide = temp_df_new_wide[, -c(1)]
      } else {
        df_mod_wide = cbind(df_mod_wide, temp_df_new_wide[, -c(1)])
      }

      curves_names = c(curves_names, paste0("Num_", i))
      curve_type_updated = c(curve_type_updated, 'num')
      regul_time_updated = append(regul_time_updated,
                                  list(regul_time_list[[i]]))
    }
  }

  # PLS
  if(print_steps){
    cat("=> PLS model.\n")
  }
  mpls = pls::plsr(Y ~ as.matrix(df_mod_wide),
                   validation = validation, jackknife = jackknife)
  if(plot_rmsep){
    plot(pls::RMSEP(mpls),
         main="Cross Validation", xlab="Components number")
  }

  opt_n_comp = which.min(pls::RMSEP(mpls)$val[1, , ])-1
  if(print_nbComp==TRUE){
    print(paste("Optimal number of PLS components : ", opt_n_comp))
  }

  if(opt_n_comp == 0){
    nb_total_cp = dim(df_mod_wide)[2]
    cat("No Optimal number of component!\n")
    cat("Search the optimal number of component EXCLUDING the intercept.\n")
    opt_n_comp = which.min(pls::RMSEP(mpls)$val[1, , ][c(2:nb_total_cp)])
    cat(paste("New \'Optimal\' number of PLS components : ", opt_n_comp, "\n"))
  }

  opt_coef = stats::coef(mpls, ncomp = opt_n_comp)

  otp_coef_list = list(stats::coef(mpls, ncomp = opt_n_comp, intercept=TRUE)[1])
  start = 0

  for(i in 1:length(curves_names)){
    end = length(regul_time_updated[[i]])
    temp_list = data.frame(time = regul_time_updated[[i]],
                           coef = opt_coef[(start + 1) : (start + end)])
    start = start + end
    otp_coef_list = append(otp_coef_list, list(temp_list))
  }
  names(otp_coef_list) = c("Intercept", curves_names)

  if(plot_reg_curves){
    for(i in 2:length(otp_coef_list)){
      plot(otp_coef_list[[i]])
      lines(otp_coef_list[[i]], col = 'blue')
      title(curves_names[i-1])
    }
  }

  naive_pls_obj = list(mpls, opt_n_comp, curves_names, opt_coef, otp_coef_list)
  names(naive_pls_obj) = c("plsr_model", "nbCP_opti", "curves_names",
                            "opti_reg_coef", "reg_obj")

  return(naive_pls_obj)

}


#' assert_multivariate_naivePLS_inputs
#'
#' This function checks the input of naivePLS function.
#'
#' @param df_list a list of dataframe (id, time, value_or_state)
#' @param Y a numeric vector
#' @param regul_time_obj a list of time regularisation values
#' @param curve_type_obj a list of curve type 'cat' or 'num'
#' @param id_col_obj a list of the names of the id columns
#' @param time_col_obj a list of the names of the time columns
#'
#' @returns a list
#'
#' @author Francois Bassac
assert_multivariate_naivePLS_inputs <- function(df_list, Y,
                                                regul_time_obj = NULL,
                                                curve_type_obj = NULL,
                                                id_col_obj = 'id',
                                                time_col_obj = 'time'){

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
    stop("naivePLS() : df_list have to be a list of dataframe
         (id, time, value).")
  }else if(mode(df_list) == "list" && N_curves != 1 &&
           !all(sapply(df_list, function(x) ncol(x)==3 ))){
    stop("naivePLS() : all df of df_list should have 3 columns
         (id, time, value).")
  }else if(N_curves == 1 && mode(df_list[[1]]) != 'list' && ncol(df_list) != 3){
    stop("naivePLS() : df_list have to be a list of dataframe
         (id, time, value).")
  }

  # Y
  if((!(is.vector(Y) && mode(Y)=="numeric") &
      !(is.list(Y) && all(sapply(Y, function(x) inherits(x, "numeric")))))){
    stop("naivePLS() : Y should be a numeric vector.")
  }

  # regul_time_list
  if(is.null(regul_time_obj)){
    stop("naivePLS() : regul_time_obj missing!")
  }
  if(mode(regul_time_obj) == "numeric"){
    regul_time_list =  obj_list_creation(N_rep = N_curves, regul_time_obj)
  } else if(mode(regul_time_obj) == "list" &&
            all(sapply(regul_time_obj, function(x) inherits(x, "numeric")))){
    regul_time_list = regul_time_obj
  } else {
    stop("naivePLS() : regul_time_obj have to be a numeric sequence or
         a list of numeric sequence.")
  }

  if(is.null(curve_type_obj)){
    stop("naivePLS() : curve_type_obj have to be either :
    'cat' or 'num' to be applied to all the curves OR
     a list of ('cat', 'num', ....) to be used for each curve.")
  } else if(all(curve_type_obj == 'cat')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(all(curve_type_obj == 'num')){
    curve_type_list = obj_list_creation(N_rep = N_curves, curve_type_obj)
  } else if(!(mode(curve_type_obj) == 'list' &&
              all(sapply(curve_type_obj, function(x) (x=='cat' | x=='num'))))){
    stop("naivePLS() : curve_type_obj have to be either \'cat\' or
         \'num\', for all curves, or a list  or a list specify 'cat' or 'num'
         for each curve.")
  } else {
    curve_type_list = curve_type_obj
  }

  # id_col_obj & time_col_obj
  if(mode(id_col_obj) == "character"){
    # Because the functions basis_list_creation works :
    id_col_list = obj_list_creation(N_rep = N_curves, id_col_obj)
  } else if(mode(id_col_obj) == "list" &&
            all(sapply(id_col_obj, function(x) inherits(x, "character")))){
    id_col_list = id_col_obj
  } else {
    stop("naivePLS() : id_col_obj have to be a character or
         a list of characters.")
  }

  if(mode(time_col_obj) == "character"){
    # Because the functions basis_list_creation works :
    time_col_list = obj_list_creation(N_rep = N_curves, time_col_obj)
  } else if(mode(time_col_obj) == "list" &&
            all(sapply(time_col_obj, function(x) inherits(x, "character")))){
    time_col_list = time_col_obj
  } else {
    stop("naivePLS () : time_col_obj have to be a character or
         a list of characters.")
  }

  assert_objs = list(N_curves, regul_time_list, curve_type_list,
                     id_col_list, time_col_list)
  names(assert_objs) = c("N_curves", "regul_time_list", "curve_type_list",
                         "id_col_list", "time_col_list")
  return(assert_objs)
}

#' naivePLS_formatting
#'
#' This function format the list of given dataframe into a time regularized
#' matrix usable for the naivePLS function or prediction.
#'
#' @param df_list a list of dataframe of 3 columns (id, time, state_or_value)
#' @param regul_time_obj a list of vector of time regularization values
#' @param curve_type_obj a list of character of the type of curves
#' @param id_col_obj a list of character for the id column names
#' @param time_col_obj a list of character for the time column names
#'
#' @returns a list
#' @export
#'
#' @author Francois Bassacs
naivePLS_formatting <- function(df_list,
                              regul_time_obj = NULL,
                              curve_type_obj = NULL,
                              id_col_obj = 'id',
                              time_col_obj = 'time'){

  # Assertions
  assert_obj = assert_multivariate_naivePLS_inputs(
    df_list = df_list, Y = c(0),
    regul_time_obj = regul_time_obj,
    curve_type_obj = curve_type_obj,
    id_col_obj = id_col_obj,
    time_col_obj = time_col_obj)

  N_curves = assert_obj$N_curves
  regul_time_list = assert_obj$regul_time_list
  curve_type_list = assert_obj$curve_type_list
  id_col_list = assert_obj$id_col_list
  time_col_list = assert_obj$time_col_list


  if(N_curves == 1 && ncol(df_list) == 3){
    df_list = list(df_list)
  }

  # Loop
  for(i in 1:N_curves){

    temp_df = df_list[[i]]
    # For all curves
    if(curve_type_list[[i]] == 'cat'){

      state_col = setdiff(names(df_list[[i]]),
                          c(id_col_list[[i]], time_col_list[[i]]))
      states = unique(df_list[[i]][state_col])

      if(all(states[[state_col]] %in% c(0, 1)) ||
         all(states[[state_col]] %in% c('0', '1')) ){
        # case : one state CFD in indicator form

        df_processed = list(df_list[[i]])
        names(df_processed) = paste0(state_col, "_1")

      } else {
        # case : multistates CFD
        df_processed = cat_data_to_indicator(df_list[[i]],
                                               id_col = id_col_list[[i]],
                                               time_col = time_col_list[[i]])
      }
      for(j in 1:length(df_processed)){

        temp_df_new = regularize_time_series(df_processed[[j]],
                                             time_seq = regul_time_list[[i]],
                                             curve_type = 'cat')
        temp_df_new_wide = convert_to_wide_format(temp_df_new,
                                                  id_col = id_col_list[[i]],
                                                  time_col = time_col_list[[i]])
        if(i == 1 && j == 1){
          df_mod_wide = temp_df_new_wide[, -c(1)]
        } else {
          df_mod_wide = cbind(df_mod_wide, temp_df_new_wide[, -c(1)])
        }
      }

    }else if(curve_type_list[[i]] == 'num'){
      # Data manipulation
      temp_df_new = regularize_time_series(temp_df,
                                           time_seq = regul_time_list[[i]],
                                           curve_type = 'num')

      temp_df_new_wide = convert_to_wide_format(temp_df_new,
                                                id_col = id_col_list[[i]],
                                                time_col = time_col_list[[i]])
      if(i == 1){
        df_mod_wide = temp_df_new_wide[, -c(1)]
      } else {
        df_mod_wide = cbind(df_mod_wide, temp_df_new_wide[, -c(1)])
      }

    }
  }
  return(df_mod_wide)

}


#' naivePLS_predict
#'
#' This function use the df_pred to make a prediction for a naivePLS object.
#'
#' @param naive_pls_obj a list of a naivePLS object
#' @param df_predict_list a list of dataframe (id, time, state_or_value)
#' @param regul_time_obj a list of time regularization values
#' @param curve_type_obj a list of curves types 'cat' or 'num'
#' @param id_col_obj a list of id column names
#' @param time_col_obj a list of time column names
#'
#' @returns a numeric vector
#' @export
#'
#' @author Francois Bassac
naivePLS_predict <- function(naive_pls_obj, df_predict_list,
                             regul_time_obj = NULL,
                             curve_type_obj = NULL,
                             id_col_obj = 'id',
                             time_col_obj = 'time'){

  df_mod_wide = naivePLS_formatting(df_list = df_predict_list,
                                   regul_time_obj = regul_time_obj,
                                   curve_type_obj = curve_type_obj,
                                   id_col_obj = id_col_obj,
                                   time_col_obj = time_col_obj)

  y_hat = predict(naive_pls_obj$plsr_model,
                  ncomp = naive_pls_obj$nbCP_opti,
                  newdata = as.matrix(df_mod_wide))

  return(y_hat)

}
