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
  names(pars) <- parameter_names
  theta_0 <- pars["theta_0"]
  eta <- pars["eta"]
  D <- x$n_days

  # define empty vectors to store probabilities
  # conditional probabilities
  pi_0u <-  c(1,rep(0,x$n_days-1))
  pi_0v <- pi_0u
  pi_1u <- pi_0u
  pi_1v <- pi_0u

  #unconditional probabilities
  psi_0u <- pi_0u
  psi_0v <- pi_0u
  psi_1u <- pi_0u
  psi_1v <- pi_0u

  # loop over days
  #for (d in 2:x$n_days){
    # estimate beta_d0s
    beta_d0s <- x$num_inf_d_unvac/x$num_susc_d_unvac
    # define theta_d
    theta_d <- theta_0 + eta * d
    # calculate beta_d1s
    beta_d1s<- theta_d * beta_d0s

    # conditional probabilities: pi_udj & pi_vdj, where
    #   j = infection status,
    #   d = day,
    #   u = unvaccinated,
    #   v = vaccinated
    pi_ud0 <- ifelse(1 - beta_d0s <= 0, 0.0001, 1 - beta_d0s)  # bound probability above 0
    pi_vd0 <- ifelse(1 - beta_d1s <= 0, 0.0001, 1 - beta_d1s)  # bound probability above 0
    pi_ud1 <- beta_d0s
    pi_vd1 <- ifelse(beta_d1s > 1, 1, beta_d1s)                # bound probability below 1
    # unconditional probabilities: psi_ju & psi_jv, where
    #   j = infection status,
    #   d = day,
    #   u = unvaccinated,
    #   v = vaccinated
    # if (d == 1){
      psi_ud0 <- pi_udo
      psi_vd0 <- pi_vd0
      psi_ud1 <- pi_ud1
      psi_vd1 <- pi_vd1
    # } else {
      psi_ud0[2:D] <- pi_ud0[2:D] * psi_ud0[1:(D-1)]          # not infected, unvaccinated
      psi_vd0[2:D] <- pi_vd0[2:D] * psi_vd0[1:(D-1)]          # not infected, vaccinated
      psi_ud1[2:D] <- pi_ud1[2:D] * psi_ud0[1:(D-1)]          # infected, unvaccinated
      psi_vd1[2:D] <- pi_vd1[2:D] * psi_vd0[1:(D-1)]          # infected, vaccinated
    # }
  #}
  # personal contribution to the likelihood ---------------------------------
  Li <- numeric(x$n)

  # Li <- ifelse(x$dinf == 999 & x$v == 0, psi_ud0[D],
  #              ifelse(x$dinf == 999 & x$v == 1, psi_vd0[D],
  #                     ifelse(x$dinf != 999 & x$v == 0, psi_ud1)))

  for (i in 1:x$n){
    if( x$dinf[i] == 999 ){
      if( x$v[i] == 0 ){ Li[i] <- psi_ud0[D] }
      if( x$v[i] == 1 ){ Li[i] <- psi_vd0[D] }
    }
    if ( x$dinf[i] != 999 ){
      if( x$v[i] == 0 ){ Li[i] <- psi_ud1[x$dinf[i]] }
      if( x$v[i] == 1 ){ Li[i] <- psi_vd1[x$dinf[i]] }
    }
  }

  return(-sum(log(Li)))
}

# ----------------------------------------------------------
#' Estimate VE using a maximum likelihood (ML) method
#'
#' This method is similar to that of Ainslie et al. 2017 (SIM). The ML method is based on calculating the
#' contribution to the likelihood function of each study participant. See XX for more details.
#' @param dat data set
#' @param n_days number of days in the study
#' @param latent_period length of latent period
#' @return list with a tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period) and the maximum likelihood estimates of the parameters.
#' @keywords wave
#' @import dplyr
#' @import tidyr
#' @export
ml_ve <- function(dat,
                  n_days,
                  latent_period = 1#,
                  #infectious_period = 4
                  ){

  N <- length(unique(dat$ID))


  for (d in 1:n_days){
    # calculate which days individuals got infected to be infectious on day d
    possible_day_of_infection <- (d  - latent_period - infectious_period):(d - latent_period)
    prev[d] <- length(which(dat$DINF_new %in% possible_day_of_infection))/N
  }
  prev <- ifelse(prev == 0, 0.001, prev)
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
               lower = c(0.0001, 0.0001, 0.0001),
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

  periods <- rep(1:n_periods, each = n_days_period)
  ve_dat <- tibble(day = 1:n_days, period = periods, ve = 1-(mle$par[2] + (mle$par[3] - 1) * .data$day)) %>%
    select(-.data$day) %>%
    group_by(.data$period) %>%
    summarise_all(.funs = mean)

  # output
  rtn <- list(param_est = param_est, ve_dat = ve_dat)

  return(rtn)
}
# ----------------------------------------------------------
