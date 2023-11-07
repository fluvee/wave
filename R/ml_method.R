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
  # alpha = pars[1]    #pars["alpha"]
  theta_0 = pars[1]  # baseline vaccine efficacy
  lambda = pars[2]   # 1- waning rate
  # for debugging
  #print(pars)

  # initialise betas
  beta_d0 <- c(1,rep(0,x$n_days))

  # initialise unconditional probabilities
  pi_id00 <- c(1,rep(0,x$n_days-1))
  pi_id01 <- c(1,rep(0,x$n_days-1))
  pi_id02 <- c(1,rep(0,x$n_days-1))

  pi_id10 <- c(1,rep(0,x$n_days-1))
  pi_id11 <- c(1,rep(0,x$n_days-1))
  pi_id12 <- c(1,rep(0,x$n_days-1))

  # initialise conditional probabilities
  psi_id00 <- c(1,rep(0,x$n_days-1))
  psi_id01 <- c(1,rep(0,x$n_days-1))
  psi_id10 <- c(1,rep(0,x$n_days-1))
  psi_id11 <- c(1,rep(0,x$n_days-1))


  #initialise period
  period <- 1
  period_start_days <- seq(1, x$n_days, by = x$n_days_period)
  # loop over days
  for (d in 1:x$n_days){
    if(d %in% period_start_days){period <- period + 1}
    #print(period)
    eta <- lambda - 1
    # estimate hazard of infection in unvaccinated for each day
    #   the ratio of the number of unvaccinated susceptible persons who became
    #   infected on day d to the total number of unvaccinated susceptible
    #   persons at the beginning of that day
    beta_d0[d] <- length(which(x$dinf == d & x$v == 0)) / length(which(x$dinf == 999 & x$v == 0))
    theta_d <- theta_0 + eta * d #period

    # conditional probabilities: pi_idvj = P(Y_idv = j|Y_i(d-1)v = 0)
    #.  i = person
    #   d = day
    #   v = vaccination status (v = 0,1)
    #   j = infection status (j = 0,1,2)

    # v = 0
    pi_id00[d] <- 1 - beta_d0[d]
    pi_id01[d] <- beta_d0[d]
    pi_id02[d] <- 0

    # v = 1
    pi_id10[d] <- 1 - theta_d * beta_d0[d]
    pi_id11[d] <- theta_d * beta_d0[d]
    pi_id12[d] <- 0

    # unconditional probabilities: psi_idvj
    if (d == 1){
      psi_id00[d] <- pi_id00[d]
      psi_id01[d] <- pi_id01[d]
      psi_id10[d] <- pi_id10[d]
      psi_id11[d] <- pi_id11[d]
    } else {
      psi_id00[d] <- pi_id00[d] * psi_id00[d-1]          # not infected, unvaccinated
      psi_id10[d] <- pi_id10[d] * psi_id10[d-1]          # not infected, vaccinated
      psi_id01[d] <- pi_id01[d] * psi_id00[d-1]          # infected, unvaccinated
      psi_id11[d] <- pi_id11[d] * psi_id10[d-1]          # infected, vaccinated
    }
  }
  # personal contribution to the likelihood ---------------------------------
  Li <- numeric(x$n)

  for (i in 1:x$n){
    if( x$dinf[i] == 999 ){
      if( x$v[i] == 0 ){ Li[i] <- psi_id00[x$n_days] }
      if( x$v[i] == 1 ){ Li[i] <- psi_id10[x$n_days] }
    }
    if ( x$dinf[i] != 999 ){
      if( x$v[i] == 0 ){ Li[i] <- psi_id01[x$dinf[i]] }
      if( x$v[i] == 1 ){ Li[i] <- psi_id11[x$dinf[i]] }
    }
  }

  Li <- ifelse(Li <= 0, 0.00001, Li)
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
  mle <- stats::optim(par = c(0.4, 1),
               fn = loglik,
               x = x,
               method = "L-BFGS-B",
               lower = c(0.0001, 1.5),
               upper = c(1, 1),
               hessian = TRUE
               # control = list(trace = 3,
               #                maxit = 1000,
               #                ndeps = 1e-4)
  )
  #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  se <- sqrt(diag(solve(mle$hessian)))

  param_est <- tibble(param = c("theta_0", "eta"), mle = 1 - mle$par, se = se,
                      lower = mle - 1.96 * se, upper = mle + 1.96 * se)

  # output
  rtn <- param_est

  return(rtn)
}
# ----------------------------------------------------------
