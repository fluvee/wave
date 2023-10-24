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
params <- params00
dat <- outcomes_dat00

# add FARI indicator variable
dat1 <- dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                       DINF_new = ifelse(DINF == 0, 999, DINF))

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
  temp <- ml_ve(dat = dat2, n_days = params$ND, n_periods = params$NJ,
                n_days_period = params$NDJ,
                latent_period = 1, infectious_period = 4)
temp3a <- temp3$ve_dat %>%
  mutate(Sim = i, Method = "ML")

ve_est <- bind_rows(ve_est,temp3a)

# proportion of sims where null hypothesis is rejected
reject_h0_ml <- reject_h0_ml + ifelse(temp3$param_est$lambda[2] > 0, 1, 0)

# mle parameter estimates
temp3b <- temp3$param_est %>%
  mutate(Sim = i,
         Method = "ML",
         value = c("mle", "quantile_2.5", "quantile_97.5"))


if (i > 1){
  mle_param_est <- bind_rows(mle_param_est,temp3b)
} else {mle_param_est <- temp3b}
}
