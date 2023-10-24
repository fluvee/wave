# ----------------------------------------------------------
#' Likelihood function for maximum likelihood (ML) method
#'
#' The ML method is based on calculating the contribution to the likelihood function of each study participant.
#' See XX for more details.
#' @param x a list with the information necessary to calculate the likelihood function
#' @param pars parameters to estimate
#' @return the negative sum of the loglikelihood function
#' @keywords wave
#' @export
loglik <- function(x, pars){
  #names(pars) <- parameter_names
  alpha = pars[1]    #pars["alpha"]
  theta_0 = pars[2]  #pars["theta_0"]
  phi = pars[3]      #pars["phi"]
  # for debugging
  #print(pars)

  # initialise unconditional probabilities
  pi_0u <- c(1,rep(0,x$n_days-1))
  pi_0v <- c(1,rep(0,x$n_days-1))
  pi_1u <- c(1,rep(0,x$n_days-1))
  pi_1v <- c(1,rep(0,x$n_days-1))

  # initialise conditional probabilities
  psi_0u <- c(1,rep(0,x$n_days-1))
  psi_0v <- c(1,rep(0,x$n_days-1))
  psi_1u <- c(1,rep(0,x$n_days-1))
  psi_1v <- c(1,rep(0,x$n_days-1))
  pi_0u  <- c(1,rep(0,x$n_days-1))

  #initialise period
  period <- 1
  period_start_days <- seq(1, x$n_days, by = x$n_days_period)
  # loop over days
  for (d in 1:x$n_days){
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
  Li <- numeric(x$n)

  for (i in 1:x$n){
    if( x$dinf[i] == 999 ){
      if( x$v[i] == 0 ){ Li[i] <- psi_0u[x$n_days] }
      if( x$v[i] == 1 ){ Li[i] <- psi_0v[x$n_days] }
    }
    if ( x$dinf[i] != 999 ){
      if( x$v[i] == 0 ){ Li[i] <- psi_1u[x$dinf[i]] }
      if( x$v[i] == 1 ){ Li[i] <- psi_1v[x$dinf[i]] }
    }
  }

  # for debugging
  #print(Li)

  # output
  return(-sum(log(Li)))
}

# ----------------------------------------------------------
#' Estimate VE using a maximum likelihood (ML) method
#'
#' This method is similar to that of Ainslie et al. 2017 (SIM). The ML method is based on calculating the
#' contribution to the likelihood function of each study participant. See XX for more details.
#' @param dat data set
#' @param n_days number of days in the study
#' @param n_periods number of periods in the study
#' @param n_days_period number of days per period
#' @param latent_period length of latent period
#' @param infectious_period length of infectious period
#' @return list with a tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period) and the maximum likelihood estimates of the parameters.
#' @keywords wave
#' @import dplyr
#' @import tidyr
#' @export
ml_ve <- function(dat, n_days, n_periods, n_days_period, latent_period = 1,
                  infectious_period = 4){

  N <- length(unique(dat$ID))
  prev <- numeric(n_days)

  for (d in 1:n_days){
    # calculate which days individuals got infected to be infectious on day d
    possible_day_of_infection <- (d  - latent_period - infectious_period):(d - latent_period)
    prev[d] <- length(which(dat$DINF_new %in% possible_day_of_infection))/N
  }
  prev <- ifelse(prev == 0, 0.0001, prev)
  x <- list(n = N,
            n_days = n_days,
            n_days_period = n_days_period,
            prev = prev,
            dinf = dat$DINF_new,
            v = dat$V
  )


  # use DE optim to get initial values
  # initial <- DEoptim(fn=logLik,
  #                    x = x,
  #                    lower = c(0.0001, 0.0001, 0.0001),
  #                    upper = c(1, 1, 2),
  #                    control = list(itermax = 100, trace = FALSE)
  # )
  # print(initial$optim$bestmem)
  # maximum likelihood estimates ----------------------------------------------
  #tryCatch({
  mle <- stats::optim(par = c(0.3, 0.4, 1),
               fn = loglik,
               x = x,
               method = "L-BFGS-B",
               lower = c(0.0001, 0.0001, 1),
               upper = c(1, 1, 2),
               hessian = TRUE
               # control = list(trace = 3,
               #                maxit = 1000,
               #                ndeps = 1e-4)
  )
  #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  se <- sqrt(diag(solve(mle$hessian)))

  param_est <- tibble(param = c("alpha", "theta_0", "lambda"), mle = mle$par - c(0,0,1), se = se,
                      lower = mle - 1.96 * se, upper = mle + 1.96 * se)

  # output
  rtn <- param_est

  return(rtn)
}
# ----------------------------------------------------------
