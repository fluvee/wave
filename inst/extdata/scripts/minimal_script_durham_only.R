# SIMVEE minimal script
# update: 15 April, 2021
# time: 12:27
# install.packages("timereg")
# install.packages("DEoptim")
library(dplyr)
# library(foreign)
# library(survival)
# library(splines)
# library(timereg)
library(ggplot2)
# library(DEoptim)
library(readr)

### Source functions from other files
# source('R/simvee.R')
# source('R/readParams.R')
# source('R/ve_methods.R')

### Load wave package
# makse sure you're in the wave root directory!
devtools::load_all()
### Read parameters from input files
#   you can specify the folder and file names of the input file within the ""
#   if no path is specified a window will pop up and allow you to choose a file from your computer
params <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_03_input.csv")
# params <- readParams("inst/extdata/input/Input_ban_400.csv")

### run simulation
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
outcomes_dat <- run_simvee(params)

### read in outcomes file
# you can specify the file name/path of the output file inside ""
# outcomes_dat <- read.csv(file.choose())

# add FARI indicator variable
outcomes_dat1 <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                                         DINF_new = ifelse(DINF == 0, 999, DINF))

# only estimate VE using Durham method
# this is so we can troubleshoot this method

reject_h0_durham <- 0 # set h0 rejection counter to 0
rtn <- list()         # empty list for output

# loop over simulations
for (i in 1:params$sim){

  # subset data for current simulation (i)
  dat <- outcomes_dat1 %>%
    filter(Sim == i)
  # fit ordinary Cox propotional hazards model
  flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data = dat)

  # test the proportional hazards assumption and compute the Schoenfeld residuals ($y)
  flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
  reject_h0_durham <- reject_h0_durham + ifelse(flu_zph$table[1,3] < 0.05, 1, 0)

  # calculate VE
  rtn[[i]] <- durham_ve(flu_zph, df = 4, n_days = params$ND, n_periods = params$NJ,
                    n_days_period = params$NDJ, var = "V") %>%
    mutate(Sim = i, Method = "Durham")
}

# bind the outputs of each simulation
results <- bind_rows(rtn) %>%
  select(Sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         lower_new = ifelse(lower < -0.5, -0.5, lower),
         upper_new = ifelse(upper > 1.5, 1.5, upper))

# plot
p <- ggplot(data = results, aes(x = period, y = mean)) + # , color = waning_rate
  geom_line() +
  geom_ribbon(aes(ymin=lower_new, ymax=upper_new), alpha = 0.1, linetype = "dashed") + # , fill = waning_rate
  #facet_wrap(~waning_rate, nrow = 3) +
  theme_bw()
p







