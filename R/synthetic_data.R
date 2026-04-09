#' beta_1_real_func
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default NULL
#'
#' @returns a value
#' @export
#'
#' @examples
#' beta_1_real_func(0)
#' beta_1_real_func(10)
#' beta_1_real_func(10:90)
#' plot(x=0:100, y=beta_1_real_func(0:100, 100), type='l', main="Beta_1")
#'
#' @author Francois Bassac
beta_1_real_func <- function(t, end_time=100, drop = NULL){
  return(beta_real = sin(t*2*pi/end_time + pi) * exp (1.5*t/end_time))
}

#' beta_2_real_func
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default 3*100/5
#'
#' @returns a value
#' @export
#'
#' @importFrom stats plogis
#'
#' @examples
#' beta_2_real_func(0)
#' beta_2_real_func(10)
#' beta_2_real_func(10:90)
#' plot(x=0:100, y=beta_2_real_func(0:100, 100), type='l', main="Beta_2")
#'
#' @author Francois Bassac
beta_2_real_func <- function(t, end_time=100, drop=3*100/5){
  return(beta_2_real = 10*t*(plogis(t-drop)-0.5)/end_time)
}

#' beta_3_real_func
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default 27
#'
#' @returns a value
#' @export
#'
#' @examples
#' beta_3_real_func(0)
#' beta_3_real_func(10)
#' beta_3_real_func(10:90)
#' plot(x=0:100, y=beta_3_real_func(0:100, 100), type='l', main="Beta_3")
#'
#' @author Francois Bassac
beta_3_real_func <- function(t, end_time=100, drop=27*100/100){
  return(beta_3_real = -(t-drop)^2/3000 + 0.2)
}

#' beta_4_real_func
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default NULL
#'
#' @returns a value
#' @export
#'
#' @examples
#' beta_4_real_func(0)
#' beta_4_real_func(10)
#' beta_4_real_func(10:90)
#' plot(x=0:100, y=beta_4_real_func(0:100, 100), type='l', main="Beta_4")
#'
#' @author Francois Bassac
beta_4_real_func <- function(t, end_time=100, drop = NULL){
  return(beta_4_real = sin(t*2*pi/end_time + 2*pi) * exp (1.5*t/end_time))
}

#' beta_5_real_func
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default 3*100/5
#'
#' @returns a value
#' @export
#'
#' @importFrom stats plogis
#'
#' @examples
#' beta_5_real_func(0)
#' beta_5_real_func(10)
#' beta_5_real_func(10:90)
#' plot(x=0:100, y=beta_5_real_func(0:100, 100), type='l', main="Beta_5")
#'
#' @author Francois Bassac
beta_5_real_func <- function(t, end_time=100, drop=3*100/5){
  return(beta_5_real = 10*(end_time-t)*(stats::plogis((end_time-t)
                                                      -drop)-0.5)/end_time)
}

#' beta_6_real_func
#'
#' Constant function = 1
#' Can adjuste the constant value by drop input
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default 1
#'
#' @returns a value
#' @export
#'
#'
#' @examples
#' beta_5_real_func(0)
#' beta_5_real_func(10)
#' beta_5_real_func(10:90)
#' plot(x=0:100, y=beta_5_real_func(0:100, 100), type='l', main="Beta_5")
#'
#' @author Francois Bassac
beta_6_real_func <- function(t, end_time=100, drop=1){
  return(beta_5_real = drop*t/t)
}

#' beta_7_real_func
#'
#' Triangular function with angle at (x=drop, y=drop) with slope of 1 and -1
#'
#' @param t evaluation time
#' @param end_time end time; default 100
#' @param drop particular point of the curve, default 3*end_time/5
#'
#' @returns a value
#' @export
#'
#'
#' @examples
#' beta_5_real_func(0)
#' beta_5_real_func(10)
#' beta_5_real_func(10:90)
#' plot(x=0:100, y=beta_5_real_func(0:100, 100), type='l', main="Beta_5")
#'
#' @author Francois Bassac
beta_7_real_func <- function(t, end_time=100, drop=3*end_time/5){
  ifelse(t < drop, t, 2*drop - t)
}


#' beta_list_generation
#'
#' @param N_states a int of the number of states wanted, default 3
#'
#' @returns a list of functions
#' @export
#'
#' @examples
#'
#' beta_list = beta_list_generation()
#' beta_list_2 = beta_list_generation(6)
#'
#' @author Francois Bassac
beta_list_generation <- function(N_states=3){
  # This function returns a list of N_states functions beta based on
  # beta_1_real_func, beta_2_real_func, beta_3_real_func

  if(N_states > 7){
    stop('beta_list_generation () can not generate more than 7 differents functions')
  }

  beta_list = list(beta_1_real_func, beta_2_real_func, beta_3_real_func,
                   beta_4_real_func, beta_5_real_func, beta_6_real_func,
                   beta_7_real_func)
  #, -beta_1_real_func, -beta_2_real_func,
  #  -beta_3_real_func, -beta_4_real_func, -beta_5_real_func,
  #  -beta_6_real_func, -beta_7_real_func)

  beta_list = beta_list[1:N_states]

  return(beta_list)
}



#### Synthetic data - single state ####

#' generate_X_df
#'
#' This function generate synthetic data of nind X(t).
#' For 'cat' curve_type it is in state 0 or 1 by two
#' different exponential laws.
#' For curve_type = 'num' it is noised cosinus.
#'
#' @param nind number of individuals, default 100
#' @param start first time, default 0
#' @param end last time, default 100
#' @param curve_type character, type of wanted synthetic data, default 'cat', need to be \'cat\' or \'num\'.
#' @param lambda_0 lambda parameter for exponential law for state 0, default 0.2
#' @param lambda_1 lambda parameter for exponential law for state 1, default 0.1
#' @param prob_start Start state 1 probability, binomial law, default 0.5
#' @param noise_sd noise added to the signal, default 0.1
#' @param seed seed for reproductability, default 123
#'
#' @returns the dataframe of the individuals
#' @export
#'
#'
#' @examples
#' generate_X_df()
#'
#' @author Francois Bassac
generate_X_df <- function(nind=100, start=0, end=100, curve_type = 'cat',
                          lambda_0=0.2, lambda_1=0.1, prob_start=0.5,
                          noise_sd = 0.1, seed = 123){
  # This function generate synthetic data of nind X(t) in state 0 or 1 by
  # 2 different exponential laws.
  # INPUTS :
  # nind  number of individuals
  # start : first time
  # end : end time
  # lambda_0 / _1 : exponential law parameters of 0 or 1
  # prob_start : probability of binomial law on start state.

  set.seed(seed = seed)

  if(is.null(curve_type)){
    stop('generate_X_df() : curve_type is missing! need to be \'cat\'
         or \'num\' ')
  }

  if(curve_type == 'num'){
    df = generate_X_df_SFD (nind = nind, start = start, end = end,
                            noise_sd = noise_sd, seed = seed)
  }else if(curve_type == 'cat'){
    df = generate_X_df_CFD(nind=nind, start=start, end=end,
                           lambda_0=lambda_0, lambda_1=lambda_1,
                           prob_start=prob_start)
  }else{
    stop('generate_X_df() : curve_type has to be either \'cat\' or \'num\'!')
  }

  return(df)
}


#' generate_Y_df
#'
#' This function generates Y_df bases on df, beta_func with the following link
#' Y = beta_0 + int(X(t)*beta(t))dt
#' It generates also the noised values of Y.
#' Here the NotS_ratio is the Noise over total Signal ratio meaning that a
#' value of 0.2 means that the noise represents 20% of the TOTAL variance.
#'
#' @param df the X(t) dataframe to evaluate Y on
#' @param curve_type a character, the type of data, default NULL, need to be \'cat\' or \'num\'.
#' @param beta_real_func_or_list a function or a list of functions beta(t), function used for the Y evaluation
#' @param beta_0_real the intercept, default 5.4321
#' @param NotS_ratio the Noise over total Signal ratio, default 0.2
#' @param seed a integer value for the seed to be reproducible 123
#' @param id_col a character of the id column, default 'id'
#' @param time_col a character of the time column, default 'time'
#' @param int_mode integration mode, 1 for integrate, 2 for pracma::trapz
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#' @param parallel a boolean to use parallelization, default FALSE
#'
#' @returns a dataframe of Y real and noised values
#' @export
#'
#'
#' @examples
#' df = generate_X_df(nind=100, curve_type = 'cat')
#' beta_real_func<-function(t, end_time=100, drop = NULL){
#' return(beta_real = sin(t*2*pi/end_time + pi) * exp (1.5*t/end_time))
#' }
#' beta_0_real=5.4321
#' Y_df = generate_Y_df(df, curve_type = 'cat',
#' beta_real_func, beta_0_real, NotS_ratio=0.2)
#'
#' @author Francois Bassac
generate_Y_df <- function(df, curve_type = NULL,
                          beta_real_func_or_list,
                          beta_0_real=5.4321,
                          NotS_ratio=0.2, seed = 123,
                          id_col = 'id', time_col = 'time',
                          int_mode = 1,
                          nb_pt = 10, subdivisions = 100,
                          parallel = FALSE){
  # This function generates Y_df bases on df, beta_func with the following link
  # Y = beta_0 + int(X(t)*beta_real_func(t))dt
  # It generates also the noised values of Y.
  # INPUT :
  # df : df to build Y for.
  # beta_real_func : link between X and Y
  # beta_0_real : intercept
  # NotS_ratio =  Noise over total Signal ratio
  # (0.10 means 10% of the TOTAL signal is noise.)

  set.seed(seed = seed)

  if(is.null(curve_type)){
    stop('generate_Y_df() : curve_type is missing! need to be \'cat\'
         or \'num\' ')
  }

  if(curve_type == 'num'){
    Y_df = generate_Y_df_SFD(df = df, beta_real_func = beta_real_func_or_list,
                      beta_0_real = beta_0_real, NotS_ratio = NotS_ratio,
                      id_col = id_col, time_col = time_col, seed = seed)
  }else if(curve_type == 'cat'){
    Y_df = generate_Y_df_CFD(df = df,
                             beta_real_func_or_list = beta_real_func_or_list,
                             beta_0_real=beta_0_real, NotS_ratio=NotS_ratio,
                             id_col = id_col, time_col = time_col,
                             int_mode = int_mode, seed = seed,
                             nb_pt = nb_pt, subdivisions = subdivisions,
                             parallel = parallel)
  }else{
    stop('generate_Y_df() : curve_type has to be either \'cat\' or \'num\'!')
  }

  return(Y_df)
}

#' number_of_test_id
#'
#' This function evaluates the number of individuals for a test set.
#'
#' @param TTRatio Train Test ratio, default 0.2
#' @param nind number of individuals, default 100
#'
#' @returns int, the number of id for test set
#' @export
#'
#' @examples
#' number_of_test_id(TTRatio = 0.2, nind=100)
#'
#' @author Francois Bassac
number_of_test_id <- function(TTRatio = 0.2, nind=100){
  TTr2 = TTRatio/(1-TTRatio)
  nind_test = floor(nind*TTr2)
  return(nind_test)
}

##### CFD data #####
#' generate_X_df_CFD
#'
#' This function generate synthetic data of nind X(t) in state 0 or 1 by two
#' different exponential laws.
#'
#' @param nind number of individuals, default 500
#' @param start first time, default 0
#' @param end last time, default 100
#' @param lambda_0 lambda parameter for exponential law for state 0, default 0.2
#' @param lambda_1 lambda parameter for exponential law for state 1, default 0.1
#' @param prob_start Start state 1 probability, binomial law, default 0.5
#' @param seed a integer, rabdom seed 123
#'
#' @returns the dataframe of the individuals
#' @export
#'
#' @importFrom stats rbinom rexp
#'
#' @examples
#' generate_X_df_CFD()
#' generate_X_df_CFD(10, 13, 60, 0.21, 0.13, 0.7)
#'
#' @author Francois Bassac
generate_X_df_CFD <- function(nind=500, start=0, end=100,
                              lambda_0=0.2, lambda_1=0.1, prob_start=0.5,
                              seed = 123){
  # This function generate synthetic data of nind X(t) in state 0 or 1 by
  # 2 different exponential laws.
  # INPUTS :
  # nind  number of individuals
  # start : first time
  # end : end time
  # lambda_0 / _1 : exponential law parameters of 0 or 1
  # prob_start : probability of binomial law on start state.

  set.seed(seed = seed)

  df = data.frame()

  df_list <- vector("list", nind)
  for(i in c(1:nind)){
    id = i
    #Initialisation
    time = c(start)
    T_c = start
    X_0 = stats::rbinom(n=1, prob = prob_start, size = 1)
    X_cs = X_0 # X_Current_state
    state = c(X_cs)

    while(T_c < end){
      if(X_cs == 0){
        lambda = lambda_0
      }else {
        lambda = lambda_1
      }
      # Exp law
      duration = stats::rexp(n=1, rate=lambda)
      # I stay in X_cs for duration
      # Update
      T_c = T_c + duration
      if(T_c > end){
        break
      }
      time = c(time, T_c)
      X_cs = (X_cs + 1)%%2
      state = c(state, X_cs)
    }
    # Add last state
    state = c(state, X_cs)
    time = c(time, end)
    df_list[[i]] <-data.frame(id = id, time=time, state=state)
  }
  df <- do.call(rbind, df_list)
  return(df)
}

#' generate_Y_df_CFD
#'
#' This function generates Y_df bases on df, beta_func with the following link
#' Y = beta_0 + int(X(t)*beta(t))dt
#' It generates also the noised values of Y.
#' Here the NotS_ratio is the Noise over total Signal ratio meaning that a
#' value of 0.2 means that the noise represents 20% of the TOTAL variance.
#'
#' @param df the X(t) dataframe to evaluate Y on
#' @param beta_real_func_or_list a function beta(t), or a list of function used for the Y evaluation
#' @param beta_0_real the intercept, default 5.4321
#' @param NotS_ratio the Noise over total Signal ratio, default 0.2
#' @param id_col a character of the id column, default 'id'
#' @param time_col a character of the time column, default 'time'
#' @param int_mode integration mode, 1 for integrate, 2 for pracma::trapz
#' @param nb_pt number of points for the integration, default value : 10
#' @param subdivisions default parameter of R function integrate;
#' default value : 100
#' @param seed a integer, random seed
#' @param parallel a boolean to enable parallel processing, default FALSE
#'
#' @returns a dataframe of Y real and noised values
#' @export
#'
#' @importFrom stats rnorm var
#' @importFrom dplyr full_join
#' @importFrom future availableCores plan multisession sequential
#' @importFrom future.apply future_lapply future_sapply
#'
#' @examples
#' df = generate_X_df(nind=100, curve_type='cat')
#' beta_real_func<-function(t, end_time=100, drop = NULL){
#' return(beta_real = sin(t*2*pi/end_time + pi) * exp (1.5*t/end_time))
#' }
#' beta_0_real=5.4321
#' Y_df = generate_Y_df_CFD(df, beta_real_func, beta_0_real)
#'
#' @author Francois Bassac
generate_Y_df_CFD <- function(df, beta_real_func_or_list,
                              beta_0_real=5.4321, NotS_ratio=0.2,
                              id_col = 'id', time_col = 'time',
                              int_mode = 1,
                              nb_pt = 10, subdivisions = 100, seed = 123,
                              parallel = FALSE){

  set.seed(seed = seed)

  if (parallel) {
    nb_cores <- max(future::availableCores() - 2, 1)
    old_plan <- future::plan(future::multisession, workers = nb_cores)
    on.exit(future::plan(old_plan), add = TRUE)
  } else {
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  Y_df = data.frame()
  state_col = setdiff(names(df), c(id_col, time_col))

  if(mode(beta_real_func_or_list) == 'list'){
    state_names = names(df)
    for(i in 1:length(state_names)){
      temp_df = data.frame(df[[state_names[i]]])
      temp_id_list <- split(temp_df, temp_df[[id_col]])
      temp_beta_f = beta_real_func_or_list[[i]]

      temp_Y_df <- do.call(rbind,
                           future.apply::future_lapply(temp_id_list,
                                                       evaluate_id_func_integral,
                                                       func=temp_beta_f,
                                                       mode = int_mode,
                                                       id_col = id_col,
                                                       time_col = time_col,
                                                       nb_pt = nb_pt,
                                                       subdivisions = subdivisions,
                                                       future.seed = TRUE))
      names(temp_Y_df) = c('id', paste0('Y_beta', i))

      if(i == 1){
        Y_df = temp_Y_df
      }else{
        Y_df = dplyr::full_join(Y_df, temp_Y_df, by = id_col)
      }
    }
    Y_df$Y_real <- rowSums(Y_df[, setdiff(names(Y_df), c(id_col))])
    Y_df$Y_real = Y_df$Y_real + beta_0_real

    if(NotS_ratio == 0 | is.null(NotS_ratio)){
      Y_df$Y_noised = Y_df$Y_real
    }else{
      noise_var = stats::var(Y_df$Y_real) * (NotS_ratio / (1-NotS_ratio))
      Y_df$Y_noised = Y_df$Y_real + stats::rnorm(n=dim(Y_df)[1],
                                                 mean=0, sd=sqrt(noise_var))
    }

  } else if(mode(beta_real_func_or_list)=='function'){
    id_list <- split(df, df[[id_col]])

    Y_df <- do.call(rbind, future.apply::future_lapply(id_list, evaluate_id_func_integral,
                                                       func=beta_real_func_or_list,
                                                       mode = int_mode,
                                                       id_col = id_col,
                                                       time_col = time_col,
                                                       nb_pt = nb_pt,
                                                       subdivisions = subdivisions,
                                                       future.seed = TRUE))

    names(Y_df) = c('id', 'Y_real')
    Y_df$Y_real = Y_df$Y_real + beta_0_real
    noise_var = stats::var(Y_df$Y_real) * (NotS_ratio / (1-NotS_ratio))
    Y_df$Y_noised = Y_df$Y_real + stats::rnorm(n=length(unique(df[, id_col])),
                                               mean=0, sd=sqrt(noise_var))
  }
  return(Y_df)
}

#' generate_X_df_test
#'
#' @param TTRatio Train Test ratio, default 0.2
#' @param nind number of individuals, default 500
#' @param start first time, default 0
#' @param end last time, default 100
#' @param curve_type  character, type of wanted synthetic data, default 'cat',
#'  need to be \'cat\' or \'num\'.
#' @param lambda_0 lambda parameter for exponential law for state 0, default 0.2
#' @param lambda_1 lambda parameter for exponential law for state 1, default 0.1
#' @param prob_start Start state 1 probability, binomial law, default 0.5
#' @param seed a integer, random seed 123
#'
#' @returns a dataframe of the test data
#' @export
#'
#' @examples
#' generate_X_df_test()
#' generate_X_df_test(TTRatio = 0.4, nind=8, start=0, end=10, curve_type = 'num')
#'
#' @author Francois Bassac
generate_X_df_test <- function(TTRatio = 0.2, nind=500, start=0, end=100,
                               curve_type = 'cat',
                               lambda_0=0.2, lambda_1=0.1, prob_start=0.5,
                               seed = 123){
  # This function generate a test set in the same way of generate_X_df.
  # TTRatio : is the train test ratio in order to have TTRatio nind for *
  # test set based on the nind of the train set.
  set.seed(seed = seed)

  TTr2 = TTRatio/(1-TTRatio)

  df_test = data.frame()

  nind_test = floor(nind*TTr2)
  df_test = generate_X_df(nind_test, start, end, curve_type,
                          lambda_0, lambda_1, prob_start, seed)

  return(df_test)
}




##### SFD data #####

#' generate_X_df_SFD
#'
#' This function generate synthetic data of nind X which are noised cosinuses.
#'
#' @param nind number of individuals, default 500
#' @param start first time, default 0
#' @param end last time, default 100
#' @param noise_sd noise added to the signal, default 0.1
#' @param seed seed for reproductability, default 123
#'
#' @returns a dataframe
#' @export
#'
#' @importFrom stats rnorm runif
#'
#' @examples
#' generate_X_df_SFD(nind = 100, start = 0, end = 100,
#' noise_sd = 0.1, seed = 123)
#'
#' @author Francois Bassac
generate_X_df_SFD <- function(nind = 100, start = 0, end = 100,
                              noise_sd = 0.1, seed = 123){
  set.seed(seed)

  time_vec <- seq(start, end, length.out = 101)

  df <- data.frame()

  df_list <- vector("list", nind)
  for (i in 1:nind) {
    # Frequency, amplitude and phase
    freq <- stats::runif(1, 0.02, 0.1)
    freq <- 1
    amplitude <- stats::runif(1, 0.5, 1.5)
    amplitude <- 2
    phase <- stats::runif(1, 0, pi/2)
    phase <- 1

    # cosinus + noise
    values <- amplitude * cos(2 * pi * freq * time_vec/max(time_vec) + phase) +
      stats::rnorm(length(time_vec), mean = 0, sd = noise_sd)

    df_list[[i]] <-data.frame(
      id = i,
      time = time_vec,
      value = values
    )
  }
  df <- do.call(rbind, df_list)
  return(df)
}

#' generate_X_df_SFD_data
#'
#' This function generate test data of nind X which are noised cosinuses.
#'
#' @param nind number of individuals, default 500
#' @param start first time, default 0
#' @param end last time, default 100
#' @param noise_sd noise added to the signal, default 0.1
#' @param seed seed for reproductability, default 123
#'
#' @returns a dataframe
#' @export
#'
#' @importFrom stats rnorm
#'
#' @examples
#' generate_X_df_SFD(nind = 100, start = 0, end = 100,
#' noise_sd = 0.1, seed = 123)
#'
#' @author Francois Bassac
generate_X_df_SFD_data <- function(nind = 100, start = 0, end = 100,
                              noise_sd = 0.1, seed = 123){
  set.seed(seed)

  time_vec <- seq(start, end, length.out = 101)

  df <- data.frame()

  df_list <- vector("list", nind)
  for (i in 1:nind) {
    # Frequency, amplitude and phase
    freq <- 1
    amplitude <- 2
    phase <- 1

    # cosinus + noise
    values <- amplitude * cos(2 * pi * freq * time_vec/max(time_vec) + phase) +
      stats::rnorm(length(time_vec), mean = 0, sd = noise_sd)

    df_list[[i]] <-data.frame(
      id = i,
      time = time_vec,
      value = values
    )
  }
  df <- do.call(rbind, df_list)
  return(df)
}

#' generate_Y_df_SFD
#'
#' This function generates Y_df bases on df, beta_func with the following link
#' Y = beta_0 + int(X(t)*beta(t))dt
#' It generates also the noised values of Y.
#' Here the NotS_ratio is the Noise over total Signal ratio meaning that a
#' value of 0.2 means that the noise represents 20% of the TOTAL variance.
#'
#' @param df the X(t) dataframe to evaluate Y on
#' @param beta_real_func a function beta(t), function used for the Y evaluation
#' @param beta_0_real the intercept, default 5.4321
#' @param NotS_ratio the Noise over total Signal ratio, default 0.2
#' @param id_col a character of the id column, default 'id'
#' @param time_col a character of the time column, default 'time'
#' @param seed a integer, random seed, 123
#'
#' @returns a dataframe
#'
#' @importFrom pracma trapz
#' @importFrom stats rnorm
#'
#' @author Francois Bassac
generate_Y_df_SFD <- function(df, beta_real_func, beta_0_real=5.4321,
                              NotS_ratio=0.2, id_col = 'id', time_col = 'time',
                              seed = 123){
  # This function generates Y_df bases on df, beta_func with the following link
  # Y = beta_0 + int(X(t)*beta_real_func(t))dt
  # It generates also the noised values of Y.
  # INPUT :
  # df : df to build Y for.
  # beta_real_func : link between X and Y
  # beta_0_real : intercept
  # NotS_ratio =  Noise over total Signal ratio
  # (0.10 means 10% of the TOTAL signal is noise.)

  value_col = setdiff(names(df), c(id_col, time_col))

  Y_df = data.frame()

  # Mode for a single function
  id_list <- split(df, df$id)
  for(i in 1:length(id_list)){

    func_values = beta_real_func(id_list[[i]][[time_col]])
    product_values = func_values * id_list[[i]][[value_col]]

    Y_temp = pracma::trapz(x = id_list[[i]][[time_col]], y = product_values)

    Y_df = rbind(Y_df, data.frame(id = i, Y_real = Y_temp))
  }

  names(Y_df) = c('id', 'Y_real')

  # Add beta_0_real
  Y_df$Y_real = Y_df$Y_real + beta_0_real

  if(NotS_ratio == 0 | is.null(NotS_ratio)){
    Y_df$Y_noised = Y_df$Y_real
  }else{
    noise_var = var(Y_df$Y_real) * (NotS_ratio / (1-NotS_ratio))
    Y_df$Y_noised = Y_df$Y_real + stats::rnorm(n=length(unique(df[, id_col])),
                                        mean=0, sd=sqrt(noise_var))
  }

  return(Y_df)
}

#### Data regularisation ####
#' regularize_time_series
#'
#' This function regularize the data on a new time interval time_seq given
#' the curve_type.
#' For curve_type = 'cat' the output state stay the same between 2 times.
#' For curve_type = 'num' the intermediate values are interpolated linearly.
#'
#' @param df dataframe with one or more different ids
#' @param time_seq New time sequence where we want to regularize
#' @param curve_type A string giving the type of the curve, 'cat' for a
#' categorical functional data, 'num' for a scalar functional data, default : NULL
#' @param id_col col_name of df for the id
#' @param time_col col_name of df for the time
#'
#' @returns a dataframe with the regularized data on time_seq
#' @export
#'
#'
#' @examples
#' id_df = data.frame(id=rep(1,5), time=seq(0, 40, 10), state=c(0, 1, 1, 0, 1))
#' regularize_time_series(id_df, time_seq = seq(0, 40, 2), curve_type = 'cat')
#' regularize_time_series(id_df, time_seq = seq(0, 40, 2), curve_type = 'num')
#'
#' @author Francois Bassac
regularize_time_series <- function(df, time_seq = 0:100, curve_type = NULL,
                                   id_col='id', time_col='time') {
  # this function regularize the data on a new time interval time_seq

  if(is.null(curve_type)){
    stop('regularize_time_series() : curve_type missing! need to be \'cat\'
         or \'num\' ')
  }

  if(curve_type == 'cat'){
    result_df = regularize_time_series_CFD(df=df, time_seq = time_seq,
                                           id_col=id_col, time_col=time_col)
  } else if(curve_type == 'num'){
    result_df = regularize_time_series_SFD(df=df, time_seq = time_seq,
                                           id_col=id_col, time_col=time_col)
  }else{
    stop('regularize_time_series() : curve_type has to be either \'cat\' or \'num\'!')
  }

  return(result_df)
}


#' regularize_time_series_CFD
#'
#' This function regularize the data CFD on a new time interval time_seq
#'
#' @param df dataframe with one or more different ids
#' @param time_seq New time sequence where we want to regularize
#' @param id_col col_name of df for the id
#' @param time_col col_name of df for the time
#'
#' @returns a dataframe with the regularized data on time_seq
#'
#' @importFrom utils tail
#'
#' @author Francois Bassac
regularize_time_series_CFD <- function(df, time_seq = 0:100,
                                       id_col='id', time_col='time') {
  # this function regularize the data on a new time interval time_seq

  state_col = setdiff(names(df), c(id_col, time_col))

  unique_ids <- unique(df[[id_col]])
  result_df <- data.frame()

  for (curr_id in unique_ids) {
    id_data <- df[df[[id_col]] == curr_id, ]
    # Order
    id_data <- id_data[order(id_data[[time_col]]), ]

    new_times <- data.frame(curr_id, time_seq)
    colnames(new_times) = c(id_col, time_col)

    new_times[[state_col]] = NA

    # Infer state for the new time_seq
    for (i in 1:nrow(new_times)) {
      t <- new_times$time[i]
      # Case 1: Before first point
      if (t < min(id_data[[time_col]])) {
        new_times[[state_col]][i] <- id_data[[state_col]][1]# First known state
      }
      # Case 2: After last knowed point
      else if (t > max(id_data[[time_col]])) {
        # Last knowed state
        new_times[[state_col]][i] <- utils::tail(id_data[[state_col]], 1)
      }
      # Case 3: On an existing time
      else if (t %in% id_data[[time_col]]) {
        new_times[[state_col]][i] <- id_data[id_data[[time_col]]== t, state_col]
      }
      # Case 4: between 2 points
      else {
        # Find times before and after
        prev_idx <- max(which(id_data[[time_col]] <= t))
        new_times[[state_col]][i] <- id_data[[state_col]][prev_idx]
      }
    }

    # Add
    result_df <- rbind(result_df, new_times)
  }

  #keeps <- setdiff(names(result_df), c(state_col))
  #result_df = result_df[keeps]

  return(result_df)
}


#' regularize_time_series_SFD
#'
#' This function regularizz a Scalar Functional Data SFD on the time_seq input.
#' This function uses linear interpolation.
#'
#' @param df dataframe with one or more different ids
#' @param time_seq New time sequence where we want to regularize
#' @param id_col col_name of df for the id
#' @param time_col col_name of df for the time
#'
#' @returns a regularized dataframe
#' @importFrom stats approx
#'
#'@author Francois Bassac
regularize_time_series_SFD <- function(df, time_seq = 0:100, id_col='id',
                                       time_col='time'){
  # This function regularize the SFD on the time_seq.

  value_col = setdiff(names(df), c(id_col, time_col))

  unique_ids <- unique(df[[id_col]])
  result_df <- data.frame()

  for (curr_id in unique_ids) {
    id_data <- df[df[[id_col]] == curr_id, ]
    # Order
    id_data <- id_data[order(id_data[[time_col]]), ]

    new_times <- data.frame(curr_id, time_seq)
    colnames(new_times) = c(id_col, time_col)

    new_times[[value_col]] = stats::approx(x = id_data[[time_col]],
                                           y = id_data[[value_col]],
                                           xout = time_seq)$y

    # Add
    result_df <- rbind(result_df, new_times)
  }

  #keeps <- setdiff(names(result_df), c(state_col))
  #result_df = result_df[keeps]

  return(result_df)

}


#' convert_to_wide_format
#'
#' This function takes the df_new output from regularize_time_series and
#' convert it into another format
#'
#' @param df_new a regularized dataframe
#' @param id_col col_name of df_new for the id
#' @param time_col col_name of df_new for the id
#'
#' @returns the dataframe in wide format
#' @export
#'
#' @importFrom tidyr pivot_wider all_of
#' @importFrom magrittr %>%
#'
#' @examples
#' id_df = data.frame(id=rep(1,5), time=seq(0, 40, 10), state=c(0, 1, 1, 0, 1))
#' id_df_new = regularize_time_series(id_df, time_seq = seq(0, 40, 2), curve_type = 'cat')
#' convert_to_wide_format(id_df_new)
#'
#' @author Francois Bassac
convert_to_wide_format <- function(df_new, id_col='id', time_col='time') {
  # This function takes the df_new output from regularize_time_series and
  # convert it into another format

  value_col = setdiff(names(df_new), c(id_col, time_col))

  # Check if all IDs have the same time
  time_counts <- table(df_new[[id_col]])
  if (length(unique(time_counts)) != 1) {
    warning("Warning! All IDs does not have the same shared time!")
  }
  df_new$time_col <- paste0("t", df_new[[time_col]])
  wide_df <- df_new %>%
    tidyr::pivot_wider(
      id_cols = tidyr::all_of(id_col),
      names_from = tidyr::all_of(time_col),
      values_from = tidyr::all_of(value_col)
    )
  return(wide_df)
}



#### Synthetic data - multi state ####

#' lambda_determination
#'
#' This function determines a nomber of lambda parameters for exponantial laws.
#'
#' @param N_states a int the number of states
#' @param lambda_values a vector of min and max values authorized for lambda
#'
#' @returns a numeric vector with the lambda values
#' @export
#'
#' @importFrom stats runif
#'
#' @examples
#' lambda_determination(3)
#' lambda_determination(7)
#'
#' @author Francois Bassac
lambda_determination <- function(N_states, lambda_values = c(0.05, 0.25)){
  # This function give some lambda values for exponential laws.
  return(stats::runif(n=N_states,
                      min=min(lambda_values),
                      max=max(lambda_values)))
}


#' generate_probabilities
#'
#' This function generates probabilities whose sum is 1.
#'
#' @param N_proba a int the number of values requested
#'
#' @returns a vector of the probabilities
#' @export
#'
#' @importFrom stats rgamma
#'
#' @examples
#' generate_probabilities(3)
#' generate_probabilities(5)
#'
#' @author Francois Bassac
generate_probabilities <- function(N_proba) {
  # This function generate probabilities whose sum is 1,
  # following Dirichlet law.
  x <- stats::rgamma(N_proba, shape = 1)
  return(x / sum(x))
}


#' transfert_probabilities
#'
#' This function gives transfer probabilities between states. row -> columns.
#'
#' @param N_states a int the number a states considered.
#'
#' @returns a dataframe containing the transition probabilities
#' @export
#'
#' @examples
#' transfert_probabilities(3)
#' transfert_probabilities(5)
#'
#' @author Francois Bassac
transfert_probabilities <- function(N_states){
  # This function gives transfer probabilities between states.

  names = c()
  for(i in 1:(N_states)){
    name = paste0("dir", i)
    names = c(names, name)
  }

  proba_df = data.frame()
  for(i in 1:N_states){
    proba_values = generate_probabilities(N_states - 1)

    proba_row = c()
    i_eq_j = FALSE
    for(j in 1:N_states){
      if(i == j){
        proba_row = c(proba_row, 0)
        i_eq_j = TRUE
      } else if(i_eq_j){
        proba_row = c(proba_row, proba_values[j-1])
      }else{
        proba_row = c(proba_row, proba_values[j])
      }
    }

    proba_df = rbind(proba_df, proba_row)
  }

  colnames(proba_df) = names
  rownames(proba_df) = paste0("state_",seq(1, N_states))
  return(proba_df)
}

#' initial_state_determination
#'
#' This functions determine randomly the initial state of a categorical
#' functional data, uniform probability between states.
#'
#' @param N_states a int the number of considered states
#'
#' @returns a value of the designated state
#' @export
#'
#' @examples
#' initial_state_determination(2)
#' initial_state_determination(7)
#'
#' @author Francois Bassac
initial_state_determination <- function(N_states){
  # This functions determine randomly the initial state of a categorical
  # functional data
  # uniform probability
  val = floor(N_states * runif(N_states, min=0, max=1) + 1)
  # return the random value
  val_choosen = sample(val, size = 1)
  return(val_choosen)
}

#' determine_next_state
#'
#' This function determines the next state base on the current state and
#' the transition matrix.
#'
#' @param current_state a value of the current state
#' @param transition_df a dataframe of the transistion matrix
#'
#' @returns a value for the next state
#' @export
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#' determine_next_state(1, transition_df)
#'
#' @author Francois Bassac
determine_next_state <- function(current_state, transition_df){
  # This function determines the next state base on the current state and
  # the transition matrix.

  # Select the right transition probabilities
  probs = transition_df[current_state, ]

  # Determine the possible next states
  col = seq(1, dim(transition_df)[1], 1)
  avaliable_states = col[col != current_state]

  if(length(avaliable_states)==1){
    next_state = avaliable_states
  }else{
    next_state = sample(avaliable_states, size = 1,
                        prob = probs[, avaliable_states])
  }


  return(next_state)
}


#' generate_X_df_multistates
#'
#' This function generates a multistates CFD.
#'
#' @param nind a int of the number of the individuals, default 100
#' @param N_states a int of the number of states wanted, default 3
#' @param start a value of starting time, default 0
#' @param end a value of ending time, default 100
#' @param lambdas a vector of N_states lambda values from lambda_determination(N_states)
#' @param transition_df a dataframe with the transition matrix from tranfert_probabilities(N_states)
#' @param seed a integer, random seed
#'
#' @returns a dataframe of a multistates Categorical Funcitonal Data
#' @export
#'
#' @importFrom stats rexp
#'
#' @examples
#' N_states = 4
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas,  transition_df)
#'
#' @author Francois Bassac
generate_X_df_multistates <- function(nind = 100, N_states=3, start=0, end=100,
                                      lambdas, transition_df, seed = 123){
  # This function generates a multistates CFD.

  set.seed(seed = seed)

  df = data.frame()
  df_list <- vector("list", nind)
  for(i in c(1:nind)){
    id = i
    #Initialisation
    time = c(start)
    T_c = start
    X_0 = initial_state_determination(N_states)
    X_cs = X_0 # X_Current_state
    state = c(X_cs)

    while(T_c < end){
      # Lambda selection
      lambda = lambdas[X_cs]

      # Exp law
      duration = stats::rexp(n=1, rate=lambda)
      # I stay in X_cs for duration
      # Update
      T_c = T_c + duration
      if(T_c > end){
        break
      }
      time = c(time, T_c)

      # New state determination
      X_cs = determine_next_state(X_cs, transition_df)

      state = c(state, X_cs)
    }
    # Add last state
    state = c(state, X_cs)
    time = c(time, end)
    df_list[[i]] <-data.frame(id = id, time=time, state=state)
  }
  df <- do.call(rbind, df_list)
  return(df)
}

#### Data transformation ####

#' state_indicatrices
#'
#' This function takes functional categorical curve as input and tranform it
#' into as many indicatrices curves as the number of state input
#' return DATAFRAME
#' Works even on dataframe without time condition respected
#' (same start and end)
#' This function sort the states by ascending order (if numeric) and put the
#' name 'state_X' as the column of the output concerning the 'X' state.
#' This function will also work with character states.
#' Now for the different lists, the ith element of a list concern the ith states
#' ordered.
#'
#' @param data a multistates dataframe ('id', 'time', 'states')
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#'
#' @returns a dataframe with columns ('id', 'time', list of states_XX)
#' @export
#' @importFrom stats setNames
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas, transition_df)
#'
#' si_df = state_indicatrices(df, id_col='id',
#' time_col='time')
#'
#' @author Francois Bassac
state_indicatrices <- function(data, id_col='id', time_col = 'time'){
  state_col = setdiff(names(data), c(id_col, time_col))
  states = sort(unique(data[[state_col]])) # Tri direct

  indicator_matrix <- sapply(states, function(s) as.numeric(data[[state_col]] == s))
  colnames(indicator_matrix) <- paste0("state_", states)

  ind <- data.frame(data[, c(id_col, time_col)], indicator_matrix)

  return(ind[order(ind[[id_col]], ind[[time_col]]), ])
}

#' state_indicatrices_old
#'
#' This function takes functional categorical curve as input and tranform it
#' into as many indicatrices curves as the number of state input
#' return DATAFRAME
#' Works even on dataframe without time condition respected
#' (same start and end)
#' This function sort the states by ascending order (if numeric) and put the
#' name 'state_X' as the column of the output concerning the 'X' state.
#' This function will also work with character states.
#' Now for the different lists, the ith element of a list concern the ith states
#' ordered.
#'
#' @param data a multistates dataframe ('id', 'time', 'states')
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#'
#' @returns a dataframe with columns ('id', 'time', list of states_XX)
#' @export
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas, transition_df)
#'
#' si_df = state_indicatrices_old(df, id_col='id',
#' time_col='time')
#'
#' @author Francois Bassac
state_indicatrices_old <- function(data, id_col='id', time_col = 'time'){
  # This function takes functional categorical curve as input and tranform it
  # into as many indicatrices curves as the number of state input
  # return DATAFRAME
  # Works even on dataframe without time condition respected
  # (same start and end).

  if(ncol(data) != 3){
    stop("state_indicatrices() : The dataframe should have 3 columns : id, time, state.")
  }else{
    state_col = setdiff(names(data), c(id_col, time_col))
  }

  states = unique(data[[state_col]])

  # I ordered the states
  states_ordered = states[order(states)]

  indicator_matrix <- sapply(states_ordered, function(state) {
    as.numeric(data[[state_col]] == state)
  })
  colnames(indicator_matrix) <- paste0("state_", states_ordered)


  ind = data.frame(id= data[[id_col]],
                   time =  data[[time_col]],
                   indicator_matrix)

  colnames(ind)[1] <- id_col
  colnames(ind)[2] <- time_col


  # I ordered by id_col and time_col.
  ind = ind[order(ind[[id_col]], ind[[time_col]] ), ]

  return(ind)
}



#' split_in_state_df
#'
#' This function transform a categorical functional data with its indicatrices
#' into a dedicated list of all the state (one per different state)
#' This function will also work with character states.
#'
#' @param data a dataframe cotaining the indicatrices, output of state_indicatrices()
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#'
#' @returns a list containing the dataframe of the indicatrice of each state.
#' @export
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas, transition_df)
#'
#' si_df = state_indicatrices(df, id_col='id',
#' time_col='time')
#'
#' split_df = split_in_state_df(si_df, id_col='id', time_col='time')
#'
#' @author Francois Bassac
split_in_state_df <- function(data, id_col='id', time_col='time'){
  # This function transform a categorical functional data with its indicatrices
  # into a dedicated list of all the state (one per different state)
  #
  # This function will also work with character states.

  col_to_use = setdiff(colnames(data), c(id_col, time_col))
  dfs = list()
  for(name in col_to_use){
    temp = data.frame(matrix(nrow=nrow(data), ncol=3))
    colnames(temp) = c(id_col, time_col, 'state')
    temp[[id_col]] = data[[id_col]]
    temp[[time_col]] = data[[time_col]]
    temp[['state']] = data[[name]] # HERE make state_col_name = state
    dfs[[name]] = temp
  }

  return(dfs)
}


#' remove_duplicate_states
#'
#' This function removes the duplicated states and keep the last line.
#' this function works both with (0,1) or categorical states on the state_col!
#'
#' @param data a single- or multi-state dataframe.
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#'
#' @returns a dataframe with no duplicated states.
#' @export
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas, transition_df)
#' remove_duplicate_states(df)
#'
#' @author Francois Bassac
remove_duplicate_states <- function(data, id_col='id',
                                    time_col='time') {
  # This function removes the duplicated states and keep the last line.
  # this function works both with (0,1) or categorical states on the state_col!

  unique_ids = unique(data[[id_col]])

  state_col = setdiff(names(data), c(id_col, time_col))

  filtered_data = data.frame()

  for(id in unique_ids){
    subset = data[data[[id_col]]==id, ]
    # Identify the changes
    changes <- c(TRUE,
                 subset[[state_col]][-1] != subset[[state_col]][-nrow(subset)])
    # Always including the last line
    changes[length(changes)] <- TRUE
    # Filtered dataset
    temp_df <- subset[changes, ]

    filtered_data = rbind(filtered_data, temp_df)
  }

  filtered_data = filtered_data[order(filtered_data[[id_col]],
                                      filtered_data[[time_col]] ), ]

  return(filtered_data)
}

#' build_df_per_state
#'
#' @param data_list a list containing the dataframe of the indicatrice of each state.
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#'
#' @returns a list of the ordered states in the indicatrice form.
#' @export
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas, transition_df)
#'
#' si_df = state_indicatrices(df, id_col='id',
#' time_col='time')
#'
#' split_df = split_in_state_df(si_df, id_col='id', time_col='time')
#'
#' @author Francois Bassac
build_df_per_state <- function(data_list, id_col='id', time_col='time'){
  # This function takes the data_list with one df per state indicatrice
  # and remove the duplicated state of each state indicatrice.
  results_df = list()
  for(name in names(data_list)){
    results_df[[name]] = remove_duplicate_states(data_list[[name]],
                                                 id_col=id_col,
                                                 time_col=time_col
    )
  }
  return(results_df)
}

#' cat_data_to_indicatrice
#'
#' This function apply all functions to go from a categorical functional data
#' with different states to a list of one dataframe per state indicatrice
#' (in the ascending order) whose duplicated states where removed.
#'
#' @param data a multistates dataframe ('id', 'time', 'states')
#' @param id_col a character for the id column, default 'id'
#' @param time_col a character for the time column, default 'time'
#'
#' @returns a list of the ordered states in the indicatrice form.
#' @export
#'
#' @examples
#' N_states = 3
#' lambdas = lambda_determination(N_states)
#' transition_df = transfert_probabilities(N_states)
#'
#' df = generate_X_df_multistates(nind = 100, N_states, start=0, end=100,
#' lambdas, transition_df)
#' df_list = cat_data_to_indicatrice(df)
#'
#' @author Francois Bassac
cat_data_to_indicatrice <- function(data, id_col='id',
                                    time_col = 'time'){
  # This function apply all functions to go from a categorical functional data
  # with different states to a list of one dataframe per state indicatrice
  # (in the ascending order) whose duplicated states where removed.

  if(ncol(data) != 3){
    stop("cat_data_to_indicatrice() : The dataframe should have 3 columns : id, time, state.")
  }else{
    state_col = setdiff(names(data), c(id_col, time_col))
  }

  temp_1 = state_indicatrices(data, id_col=id_col,
                              time_col=time_col)
  temp_2 = split_in_state_df(temp_1, id_col=id_col, time_col=time_col)
  temp_3 = build_df_per_state(temp_2, id_col=id_col, time_col=time_col)
  return(temp_3)
}

