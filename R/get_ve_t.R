#' function to obtain VE(t) from the data by bootstrapping
#' @param n
#' @param data
#' @param indices
#' @return yhat
#' @keywords wave
#' @import splines
#' @import dplyr
#' @import tidyr
#' @export

get_ve_t <- function(n, data, indices) {
  d <- data[indices,] # allows boot to select sample
  coxmod <- coxph(Surv(DINF_new,influenza) ~ V, data = d) # fit ordinary Cox propotional hazards model (just for IIV first)
  zph <- cox.zph(fit = coxmod, transform = "identity") # test the proportional hazards assumption and compute the Schoenfeld residuals ($y)
  xx <- zph$x # time points
  yy <- zph$y # Schoenfeld residuals
  #d <- nrow(yy)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = n)
  #seq(from = (n_days_period/2), to = n_days - (n_days_period/2), length = n_periods)
  temp <- c(pred.x, xx)
  lmat <- ns(temp, df = 2, intercept = TRUE)
  pmat <- lmat[1:length(pred.x), ]
  xmat <- lmat[-(1:length(pred.x)), ]
  qmat <- qr(xmat)

  # ve estimate
  yhat.beta <- pmat %*% qr.coef(qmat, yy)
  yhat <- 1-exp(yhat.beta) # VE estimate

  return(yhat)
}
