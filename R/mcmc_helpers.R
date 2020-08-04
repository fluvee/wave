# ----------------------------------------------------------
#' Likelihood function for maximum likelihood (ML) method
#'
#' The ML method is based on calculating the contribution to the likelihood function of each study participant.
#' See XX for more details.
#' @param x a list with the information necessary to calculate the likelihood function
#' @param alpha alpha parameter
#' @param theta_0 theta_0 parameter
#' @param phi phi parameter
#' @return the negative sum of the loglikelihood function
#' @keywords wave
#' @export
likelihood_xxx <- function(x, alpha, theta_0, phi){
  pi_0u <- pi_0v <- pi_1u <- pi_1v <- c(1,rep(0,x$n_days-1))
  psi_0u <- psi_0v <- psi_1u <- psi_1v <- pi_0u

  #initialise period
  period <- 1
  period_start_days <- seq(1, x$n_days, by = x$n_days_period)
  # loop over days
  for (d in 2:x$n_days){
    if(d %in% period_start_days){period <- period + 1}
    #print(period)
    lambda <- phi - 1
    theta_d <- theta_0 + lambda * period
    alpha_d <- alpha * x$prev[d]
    # conditional probabilities: pi_ju & pi_jv, where
    #   j = infection status,
    #   u = unvaccinated,
    #   v = vaccinated
    pi_0u[d] <- ifelse(1 - alpha_d < 0, 0.0001, 1 - alpha_d)
    pi_0v[d] <- ifelse(1 - alpha_d * theta_d < 0, 0.0001, 1 - alpha_d * theta_d)
    pi_1u[d] <- alpha_d
    pi_1v[d] <- ifelse(alpha_d * theta_d > 1, 1, alpha_d * theta_d)
    # unconditional probabilities: psi_ju & psi_jv, where
    #   j = infection status,
    #   u = unvaccinated,
    #   v = vaccinated
    if (d == 1){
      psi_0u[d] <- pi_0u[d]
      psi_0v[d] <- pi_0v[d]
      psi_1u[d] <- pi_1u[d]
      psi_1v[d] <- pi_1v[d]
    } else {
      psi_0u[d] <- pi_0u[d] * psi_0u[d-1]          # not infected, unvaccinated
      psi_0v[d] <- pi_0v[d] * psi_0v[d-1]          # not infected, vaccinated
      psi_1u[d] <- pi_1u[d] * psi_0u[d-1]          # infected, unvaccinated
      psi_1v[d] <- pi_1v[d] * psi_0v[d-1]          # infected, vaccinated
    }
  }
  # personal contribution to the likelihood ---------------------------------
  lik <- numeric(x$n)

  for (i in 1:x$n){
    if( x$dinf[i] == 999 ){
      if( x$v[i] == 0 ){ lik[i] <- psi_0u[x$n_days] }
      if( x$v[i] == 1 ){ lik[i] <- psi_0v[x$n_days] }
    }
    if ( x$dinf[i] != 999 ){
      if( x$v[i] == 0 ){ lik[i] <- psi_1u[x$dinf[i]] }
      if( x$v[i] == 1 ){ lik[i] <- psi_1v[x$dinf[i]] }
    }
  }
  rtn <- sum(log(lik))
  return(rtn)
}

# ----------------------------------------------------------
#' Likelihood function for maximum likelihood (ML) method
#'
#' The ML method is based on calculating the contribution to the likelihood function of each study participant.
#' See XX for more details.
#' @param parTab The parameter table controlling information such as bounds, initial values etc.
#' @param data The data frame of data to be fitted.
#' @param PRIOR_FUNC Prior distribution function. If no priors are specified (PRIOR_FUNC = NULL), then the uniform prior of the sampling restrictions will apply from parTab.
#' @param ... other function arguments
#' @return the negative sum of the loglikelihood function
#' @keywords wave
#' @export
## You MUST specify the arguments to this function as parTab, data then PRIOR_FUNC.
## Use the `...` to pass additional arguments through.
my_creation_function <- function(parTab, data, PRIOR_FUNC, ...){
  ##############################
  ## This is where you would manipulate all
  ## of your model arguments. For example,
  ## unpackaging the stuff from `...`,
  ## or transforming model parameters
  ##############################

  ## Somewhat redundant example
  parameter_names <- parTab$names

  N <- length(unique(data$ID))
  n_days <- parTab[parTab$names == "n_days",]$values
  n_periods <- parTab[parTab$names == "n_periods",]$values
  n_days_period <- parTab[parTab$names == "n_days_period",]$values
  latent_period <- parTab[parTab$names == "latent_period",]$values
  infectious_period <- parTab[parTab$names == "infectious_period",]$values

  prev <- numeric(n_days)

  for (d in 1:n_days){
    # calculate which days individuals got infected to be infectious on day d
    possible_day_of_infection <- (d  - latent_period - infectious_period):(d - latent_period)
    prev[d] <- length(which(data$DINF_new %in% possible_day_of_infection))/N
  }
  prev <- ifelse(prev == 0, 0.001, prev)

  x <- list(n = N,
            n_days = n_days,
            n_days_period = n_days_period,
            prev = prev,
            dinf = data$DINF_new,
            v = data$V
  )

  ##############################
  ## This is where you would put your own model code
  ##############################
  likelihood_func <- function(pars){
    #names(pars) <- parameter_names
    alpha = pars[1]
    theta_0 = pars[2]
    phi = pars[3]

    ## Note the use of closures to capture the `data` argument
    loglik <- likelihood_xxx(x, alpha, theta_0, phi) #sum(dnorm(data, mu, sd, TRUE))
    if(!is.null(PRIOR_FUNC)) loglik <- loglik + PRIOR_FUNC(pars)
    loglik
  }
  return(likelihood_func)
}
