# SIMVEE minimal script
# update: 11 May, 2020
# time: 18:02
# install.packages("timereg")
# install.packages("DEoptim")
# library(dplyr)
# library(foreign)
# library(survival)
# library(splines)
# library(timereg)
# library(ggplot2)
# library(DEoptim)

### Source functions from other files
# source('R/simvee.R')
# source('R/readParams.R')
# source('R/ve_methods.R')

### Load wave package

# make sure you're in the wave root directory!
devtools::load_all()
### Read parameters from input files
#   you can specify the folder and file names of the input file within the ""
#   if no path is specified a window will pop up and allow you to choose a file from your computer
params <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_00_input.csv")

### run simulation
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
outcomes_dat <- run_simvee(params)

### read in outcomes file
# you can specify the file name/path of the output file inside ""
# outcomes_dat <- read.csv(file.choose())

# add FARI indicator variable
outcomes_dat <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                                        DINF_new = ifelse(DINF == 0, 999, DINF))

# define input for ML method with MCMC
parTab <- data.frame(values=c(params$alpha_0, params$theta_d[1], params$theta_d[1] - params$theta_d[2] + 1,
                              params$ND, params$NJ, params$NDJ, 1, 4),
                     names=c("alpha","theta_0","phi", "n_days", "n_periods", "n_days_period",
                             "latent_period", "infectious_period"),
                     fixed=c(0,0,0,rep(1,5)),
                     steps=c(rep(0.01,8)),
                     lower_bound=c(rep(0.0001,3),rep(0,5)),
                     upper_bound=c(1,1,2, rep(1000, 5)),
                     stringsAsFactors=FALSE)

## Update the proposal step size every 1000 iterations (opt_freq) for the first 5000 iterations
## (adaptive_period) to aim for an acceptance rate of 0.44 (popt). After the adaptive period,
## run for a further 10000 (iterations) steps. Save every nth rows, where n is "thin" (ie. 1 here).
## Write to disk every 100 iterations (post thinning). Note that the adaptive period is also saved
## to disk
mcmcPars <- c("iterations"=10000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=5000,"save_block"=1000)

# apply VE estimation methods
ve_estimates <- estimate_ve(dat = outcomes_dat, params, par_tab = parTab, mcmc_pars = mcmcPars)

# print proportion of null hypotheses rejected
ve_estimates$prop_reject_h0

# print mean mle parameter estimates
ve_estimates$mean_mle_params




