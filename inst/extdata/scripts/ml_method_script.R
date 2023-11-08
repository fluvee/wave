### ML method script ###
# ------------------------------------------------------------------------------
# We will use the ML method to estimate the rate of waning, first from simulated
# data and then from RCT data from the Michigan trial.
# ------------------------------------------------------------------------------

# Load required packages
library(dplyr)
library(foreign)
#library(survival)
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
params_test <- readParams("./inst/extdata/input/test_input_new.csv")

params00 <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_00_input.csv")
params03 <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_03_input.csv")
params06 <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_06_input.csv")

# Run simulation (or alternatively read in previously output sim outcomes file)
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
outcomes_test <- run_simvee(params_test, path = "./inst/extdata/output")

outcomes_dat00 <- run_simvee(params00, path = "./inst/extdata/output")
outcomes_dat03 <- run_simvee(params03, path = "./inst/extdata/output")
outcomes_dat06 <- run_simvee(params06, path = "./inst/extdata/output")

# Read in outcomes files
#   you can specify the file name/path of the output file inside ""
# outcomes_dat00 <- read.csv("./inst/extdata/output/Outcomes_sim_test_n_5000.csv")
# outcomes_dat03 <- read.csv(
#   paste0("./inst/extdata/output/Outcomes_",params03$title,".csv")
#   )
# outcomes_dat06 <- read.csv("./inst/extdata/output/Outcomes_MI_RCT_06_04_06.csv")

# select which input/putput files to use for estimation
my_params <- params_test
my_dat <- outcomes_test

# add FARI indicator variable
dat1 <- my_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                       DINF_new = ifelse(DINF == 0, 999, DINF))

# ------------------------------------------------------------------------------

# Estimate waning VE for each simulation ---------------------------------------
# loop through simulations
for (i in 1:my_params$sim){
  print(i)

  # subset data for each simulation
  dat2 <- dat1 %>% filter(.data$Sim == i)

  # estimate ve waning using ML method
  #   we'll use the version without MCMC
  #   the parameters estimated are:
  #   alpha = pars["alpha"]
  #   theta_0 = pars["theta_0"]
  #   phi = pars["phi"]
  param_estimates <- ml_ve(dat = dat2,
                           n_days = my_params$ND,
                           n_periods = my_params$NJ,
                           n_days_period = my_params$NDJ#,
                           # latent_period = 1,
                           # infectious_period = 4
                           )

  # store parameter estimates for each simulation
  # mle parameter estimates
  mle_est <- param_estimates %>%
    mutate(Sim = i,Method = "ML") %>%
    select(Sim, param, mle, Method)

  # estimate VE from MLE parameters for each day
  ve_dat <- tibble(day = 1:my_params$ND,
                   period = rep(1:my_params$NJ, each = my_params$NDJ),
                   ve = mle_est$mle[1] + (mle_est$mle[2] * .data$day)
                   )

  # store daily VE estimates for each simulation
  daily_ve <- ve_dat %>%
    mutate(Sim = i, Method = "ML") %>%
    select(Sim, day, period, ve, Method) %>%
    # if VE is negative, change to 0
    mutate(ve = if_else(ve < 0, 0, ve))

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
            median = median(mle),
            q025 = quantile(mle, probs = c(0.025)),
            q975 = quantile(mle, probs = c(0.975))
            )
mean_mle_est
# get mean and 95% confidence bounds for daily VE
mean_ve_est <- df_ve_est %>%
  group_by(day) %>%
  summarise(mean = mean(ve),
            median = median(ve),
            q025 = quantile(ve, probs = c(0.025)),
            q05  = quantile(ve, probs = c(0.05)),
            q25  = quantile(ve, probs = c(0.25)),
            q75  = quantile(ve, probs = c(0.75)),
            q95  = quantile(ve, probs = c(0.90)),
            q975 = quantile(ve, probs = c(0.975))
  )

# plot of mean VE and confidence bounds
ve_plot <- ggplot(data = mean_ve_est, aes(x = day, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1) +
  labs(y = "VE(t)", x = "Day",) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(yintercept = 1 - theta_0_true, linetype = "dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ve_plot
ggsave("./inst/extdata/output/ve_sim_plot_00.jpg", ve_plot)
# # line plot of all VE estimates
# ve_plot_all <- ggplot(data = df_ve_est, aes(x = day, y = ve, color = as.factor(Sim))) +
#   geom_line() +
#   geom_hline(yintercept = 1 - theta_0_true, linetype = "dashed") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "none")
# ve_plot_all


### Debugging
# look at the simulations where mle of lambda is positive (for case where
# waning rate = 0)
# find simulation numbers where lambda_hat > 0
lambda_gt_0 <- df_mle_est %>% filter(param == "lambda", mle > 0)

# compare sims where lambda_hat = 0 and lambda_hat > 0
sim1 <- dat1 %>% filter(Sim == 1)    # lambda_hat = 0
sim2 <- dat1 %>% filter(Sim == 194)    # lambda_hat > 0

hist(sim1$DINF[sim1$DINF > 0])
hist(sim2$DINF[sim2$DINF > 0])
