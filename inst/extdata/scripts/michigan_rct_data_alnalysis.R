# -------------------------------------------------------------------------------------
# Michigan clinical trial analysis script #
# -------------------------------------------------------------------------------------

# load necessary packages -------------------------------------------------------------
library(readr)
library(tidyr)
library(dplyr)
library(survival)
library(timereg)
library(lazymcmc)
library(foreign)
library(splines)
library(ggplot2)

# load data ---------------------------------------------------------------------------
mi_data <- read_csv("inst/extdata/data/Michigan_0708.csv")

# data wrangling ----------------------------------------------------------------------
mi_data1 <- mi_data %>%
  mutate(V = ifelse(vaccine_type != "PLACEBO" , 1, 0),
         DINF = onset_date - min(mi_data$vaccination_date),
         DINF_new = ifelse(is.na(DINF), 999, DINF))

# filter just IIV and placebo
mi_data_iiv <- mi_data1 %>%
  filter(vaccine_type != "LAIV")
# filter just LAIV and placebo
mi_data_laiv <- mi_data1 %>%
  filter(vaccine_type != "IIV")

# method from Durham et al. 1988 ------------------------------------------------------
# This method is based on smoothing scaled residuals from the Cox proportional hazard
# regression model. It consists of four steps. First, an ordinary proportional hazard
# model is fitted using a partial likelihood function. Second, Schoenfeld residuals are
# calculated. These residuals are used to test the independence between residuals and
# time. Third, the residuals are scaled and added to the coefficient from the ordinary
# proportional hazard regression model. Fourth, after smoothing we can get the estimated
# hazard ratios as function of time by recovering the time-varying regression coefficient.
# This allows testing the hypothesis of no VE waning and the estimation of TVE at each
# time point.

# fit ordinary Cox propotional hazards model (just for IIV first)
flu_coxmod <- coxph(Surv(DINF_new,influenza) ~ V, data = mi_data_iiv)

# test the proportional hazards assumption and compute the Schoenfeld residuals
flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
# hypothesis test of whether beta varies over time
flu_zph$table[1,3]

# calculate VE and CI using method from Petrie et al 2016 (JID) ------------------------
# absolute VE = 1 - exp(beta)
abs_ve <- 1 - exp(flu_coxmod$coefficients[1])
# VE(t)
time_points <- flu_zph$x # time points
schoenfeld <- flu_zph$y # Schoenfeld residuals

# fit loess to Schoenfeld residuals + beta
# this approach is from Petrie et al. 2016 (JID)
for_loess <- data.frame(x = time_points,
                        V = schoenfeld[,1] + flu_coxmod$coefficients[1])

l <- loess(V~x, data = for_loess, degree = 1)
ve_t <- 1 - exp(l$fitted) # VE(t)
# get confidence intervals
pred <- predict(l, for_loess, se = TRUE)
ve_t_lower <- 1- exp(pred$fit - 1.96*pred$se.fit)
ve_t_upper <- 1 - exp(pred$fit + 1.96*pred$se.fit)
# integrate loess function and divide by days between first and last case
# of symptomatic influenza to get average VE(t)
f <- function(x) predict(l, newdata = x)
int <- integrate(f, min(time_points), max(time_points))
ve_t_avg <- 1 - exp(int$value / (max(time_points) - min(time_points)))
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
# output (save for each vaccine)
# IIV
first_date <- mi_data_iiv[which(mi_data_iiv$DINF_new == min(time_points)),]$onset_date
pred_dates <- first_date + (time_points - min(time_points))
ve_dat_iiv <- tibble(time_point = time_points,
                     date = pred_dates,
                     ve = ve_t,
                     lower = ve_t_upper,
                     upper = ve_t_lower,
                     vac_type = "IIV")

# LAIV
first_date <- mi_data_laiv[which(mi_data_laiv$DINF_new == min(time_points)),]$onset_date
pred_dates <- first_date + (time_points - min(time_points))
ve_dat_laiv <- tibble(time_point = time_points,
                      date = pred_dates,
                      ve = ve_t,
                      lower = ve_t_upper,
                      upper = ve_t_lower,
                      vac_type = "LAIV")
# for output
ve_output <- bind_rows(ve_dat_iiv, ve_dat_laiv)
write_csv(ve_output, file = "VE_est_durham.csv")

# combine data sets for plotting and output
ve_dat_plot <- ve_output %>%
  filter(vac_type == "IIV")
  #mutate(lower = ifelse(lower < -0.25, -0.25, lower))

p_petrie <- ggplot(data = ve_dat_plot, aes(x = date, y = ve)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
  labs(y = "VE(t)", x = "Date", title = "Loess") +
  scale_y_continuous(limits = c(-5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() #+
  #facet_grid(. ~ vac_type)
p_petrie

# ---------------------------------------------------------------------------------------
# calculate VE and CIs using bootstrap --------------------------------------------------
library(boot)
# function to obtain VE(t) from the data
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
# bootstrapping with 1000 replications
# n_days <- max(time_points) - min(time_points)
# n_days_period <- 7
# n_period <- round(n_days/n_days_period)
results <- boot(data=mi_data_iiv, statistic=get_ve_t,
                R=1000, n = 63
                # n_days_period = n_days_period,
                # n_days = n_days, n_period = n_period
                )

# view results
#results
plot(results)

# get 95% confidence interval
# boot.ci(results, index = 1)
first_date <- mi_data_iiv[which(mi_data_iiv$DINF_new == min(time_points)),]$onset_date
pred_dates <- seq.Date(from = first_date, to = first_date + (max(time_points) - min(time_points)), length.out = 63)

boot_dat <- data.frame(date = pred_dates,
                       V = results$t0,
                       lower = apply(results$t, 2, quantile, probs = 0.025),
                       upper = apply(results$t, 2, quantile, probs = 0.975)
                       )
fit <- lm(V ~ date, boot_dat)

# plots
plot(V ~ date, data = boot_dat)
abline(fit)

p_durham <- ggplot(data = boot_dat, aes(x = date, y = V)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.1) +
  labs(y = "VE(t)", x = "Date") + # , title = "Natural Spline"
  #scale_y_continuous(limits = c(-5, 1)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_grid(. ~ vac_type) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
p_durham

ggsave(filename = "inst/extdata/output/figure2.jpg", plot = p_durham)
p_both <- plot_grid(p_durham, p_petrie)
p_both
