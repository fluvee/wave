# ----------------------------------------------------------
#' Estimate VE using method from Durham et al. 1988 (AJE)
#'
#' This method is based on smoothing scaled residuals from the Cox proportional hazard regression model.
#' It consists of four steps. First, an ordinary proportional hazard model is fitted using a partial likelihood
#' function. Second, Schoenfeld residuals are calculated. These residuals are used to test the independence
#' between residuals and time. Third, the residuals are scaled and added to the coefficient from the ordinary
#' proportional hazard regression model. Fourth, after smoothing we can get the estimated hazard ratios as
#' function of time by recovering the time-varying regression coefficient. This allows testing the hypothesis of
#' no VE waning and the estimation of TVE at each time point.
#'
#' @param x survival object
#' @param df degrees of freedom
#' @param n_days number of days in the study
#' @param n_periods number of periods in the study
#' @param n_days_period number of days per period
#' @param var name of variable to assess
#' @return tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period.
#' @keywords wave
#' @import splines
#' @import dplyr
#' @import tidyr
#' @export
durham_ve <- function(x, df = 2, n_days, n_periods, n_days_period, var){
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  if (length(n_days) == 1) {
    pred.x <- seq(from = min(xx), to = max(xx), length = n_days)
  } else {pred.x <- n_days}
  temp <- c(pred.x, xx)
  lmat <- ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:length(pred.x), ]
  xmat <- lmat[-(1:length(pred.x)), ]
  qmat <- qr(xmat)
  if (x$transform!="identity")
    stop("please re-fit the Cox model with the identity transform")
  if (qmat$rank < df)
    stop("Spline fit is singular, try a smaller degrees of freedom")
  # se
  bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
  xtx <- bk %*% t(bk)
  seval <- ((pmat %*% xtx) * pmat) %*% rep(1, df)

  # ve estimate
  yhat.beta <- pmat %*% qr.coef(qmat, yy)
  yhat <- 1-exp(yhat.beta) # VE estimate

  # se
  temp <- 2 * sqrt(x$var[1,1] * seval)
  ylow.beta <- yhat.beta - temp
  yup.beta <- yhat.beta + temp
  ylow <- 1-exp(yup.beta)
  yup <- 1-exp(ylow.beta)
  yr <- range(yhat, yup, ylow)

  # output
  periods <- rep(1:n_periods, each = n_days_period)
  ve_dat <- tibble(day = round(pred.x), period = periods, ve = as.numeric(yhat),
                   ve_lower = ylow, ve_upper = yup) %>%
    select(-.data$day) %>%
    group_by(.data$period) %>%
    summarise_all(.funs = mean)

  return(ve_dat)
}
# ----------------------------------------------------------
