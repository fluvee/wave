### ML method script ###
# ------------------------------------------------------------------------------
# We will use the ML method to estimate the rate of waning, first from simulated
# data and then from RCT data from the Michigan trial.
# ------------------------------------------------------------------------------

# Load required packages
library(dplyr)
library(foreign)
library(survival)
library(splines)
#library(timereg)
library(ggplot2)
#library(lazymcmc)
library(openxlsx)
# Load wave package
#   make sure you're in the wave root directory!
devtools::load_all()
# ------------------------------------------------------------------------------

# Simulated data ---------------------------------------------------------------
# Read parameters from input files
#   you can specify the folder and file names of the input file within the ""
#   if no path is specified a window will pop up and allow you to choose a file
#   from your computer

params00 <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_00_input.csv")
params03 <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_03_input.csv")
params06 <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_06_input.csv")

# Run simulation (or alternatively read in previously output sim outcomes file)
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
# outcomes_dat <- run_simvee(params, path = "/inst/extdata/output/")

# Read in outcomes files
#   you can specify the file name/path of the output file inside ""
outcomes_dat00 <- read.csv("./inst/extdata/output/Outcomes_MI_RCT_06_04_00.csv")
outcomes_dat03 <- read.csv("./inst/extdata/output/Outcomes_MI_RCT_06_04_03.csv")
outcomes_dat06 <- read.csv("./inst/extdata/output/Outcomes_MI_RCT_06_04_06.csv")

# select which input/putput files to use for estimation
my_params <- params00
my_dat <- outcomes_dat00
lambda_true <- 0

# add FARI indicator variable
dat1 <- my_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                       DINF_new = ifelse(DINF == 0, 999, DINF))

# initialise
reject_h0 <- 0
# ------------------------------------------------------------------------------

# Estimate waning VE for each simulation ---------------------------------------
# loop through simulations
for (i in 1:params$sim){
  print(i)

  # subset data for each simulation
  dat2 <- dat1 %>% filter(.data$Sim == i)

  # estimate ve waning using ML method
  #   we'll use the version without MCMC
  #   the parameters estimated are:
  #   alpha = pars["alpha"]
  #   theta_0 = pars["theta_0"]
  #   phi = pars["phi"]
  param_estimates <- ml_ve(dat = dat2, n_days = params$ND, n_periods = params$NJ,
    n_days_period = params$NDJ, latent_period = 1, infectious_period = 4)

  # store parameter estimates for each simulation
  # mle parameter estimates
  mle_est <- param_estimates %>%
    mutate(Sim = i,Method = "ML") %>%
    select(Sim, param, mle, Method)

  # estimate VE from MLE parameters for each day
  ve_dat <- tibble(day = 1:n_days,
                   period = rep(1:n_periods, each = n_days_period),
                   ve = 1-(mle$par[2] + (mle$par[3] - 1) * .data$day)
                   )

  # store daily VE estimates for each simulation
  daily_ve <- ve_dat %>%
    mutate(Sim = i, Method = "ML") %>%
    select(Sim, day, period, ve, Method)

  # amend results for each successive simulation
  if (i > 1){
    df_mle_est <- bind_rows(df_mle_est, mle_est)
    df_ve_est <- bind_rows(df_ve_est, daily_ve)
  } else {
    df_mle_est <- mle_est
    df_ve_est <- daily_ve
    }

  # calculate number of sims where null hypothesis is rejected
  # that is, if lambda is significantly different from its true value
  # the null hypothesis is ACCEPTED if (lower, upper) CONTAINS the true value
  # reject_h0 <- reject_h0 +
  #   ifelse(lambda_true %in% param_estimates$lower[3]:param_estimates$lower[3], 0, 1)
}

# Post process results after MLE estimation for each simulation ----------------
# get mean MLE and 95% bounds for each parameter
mean_mle_est <- df_mle_est %>%
  group_by(param) %>%
  summarise(mean = mean(mle),
            q025 = quantile(mle, probs = c(0.025)),
            q975 = quantile(mle, probs = c(0.975))
            )

# get mean and 95% confidence bounds for daily VE
mean_ve_est <- df_ve_est %>%
  group_by(day) %>%
  summarise(mean = mean(ve),
            q025 = quantile(ve, probs = c(0.025)),
            q975 = quantile(ve, probs = c(0.975))
  )

ve_plot <- ggplot(data = mean_ve_est, aes(x = day, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1) +
  labs(y = "VE(t)", x = "Day",) +
  #scale_y_continuous(limits = c(-0.5, 1)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
