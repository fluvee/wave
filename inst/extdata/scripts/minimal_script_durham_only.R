# SIMVEE minimal script --------------------------------------------------------
# update: 16 August, 2021
# time: 12:27

# install required packages ----------------------------------------------------
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

### Source functions from other files ------------------------------------------
# source('R/simvee.R')
# source('R/readParams.R')
# source('R/ve_methods.R')

### Load wave package ----------------------------------------------------------
# make sure you're in the wave root directory!
devtools::load_all()
### Read parameters from input files -------------------------------------------
#   you can specify the folder and file names of the input file within the ""
#   if no path is specified a window will pop up and allow you to choose a file
#   from your computer
params <- readParams("./inst/extdata/input/SimVEE_MI_RCT_06_04_06_input.csv")
params$N <- 5000
### run simulation -------------------------------------------------------------
#   there is an optional path argument for run_simvee(params, path = )
#   if no path is specified, it will default to current working directory
params$sim <- 100 # change number of sims (if necessary)
outcomes_dat <- run_simvee(params, path = "inst/extdata/output/")

### read in outcomes file ------------------------------------------------------
# you can specify the file name/path of the output file inside ""
file_path <- "inst/extdata/output/Outcomes_MI_RCT_06_04_06.csv"
outcomes_dat <- read_csv(file_path)

# add FARI indicator variable
outcomes_dat1 <- outcomes_dat %>% mutate(FARI = ifelse(DINF == 0, 0, 1),
                                         DINF_new = ifelse(DINF == 0, 999, DINF))

### estimate VE using Durham method --------------------------------------------
# this is so we can troubleshoot this method

reject_h0_durham <- 0 # set h0 rejection counter to 0
rtn <- list()         # empty list for output
slopes <- data.frame(sim = 1:params$sim, slope = c(rep(NA, params$sim)))
# loop over simulations
for (i in 1:params$sim){

  # subset data for current simulation (i)
  dat <- outcomes_dat1 %>%
    filter(Sim == i)
  # fit ordinary Cox propotional hazards model
  flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data = dat)

  # test the proportional hazards assumption and compute the Schoenfeld
  # residuals ($y)
  flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
  reject_h0_durham <- reject_h0_durham + ifelse(flu_zph$table[1,3] < 0.05, 1, 0)

  # calculate VE
  end_point <- params$ND - (params$NDJ/2)
  if (max(flu_zph$x) <= end_point){
    rtn[[i]] <- NA
  } else {
    rtn[[i]] <- durham_ve(flu_zph, df = 4, n_days = params$ND,
                          n_periods = params$NJ, n_days_period = params$NDJ,
                          var = "V") %>%
      mutate(Sim = i, Method = "Durham")
    fit <- lm(V ~ period, rtn[[i]])
    slopes [i,2] <- fit$coefficients[2]
  }
}
# ------------------------------------------------------------------------------

# results ----------------------------------------------------------------------

# what proportion of the time do we reject H0?
prop_reject <- reject_h0_durham/params$sim
prop_reject
# bind the outputs of each simulation
# VE estimates
find_nas <- is.na(rtn)
length(which(find_nas == TRUE))
rtn1 <- rtn[!find_nas]
results <- list(ve = bind_rows(rtn1),
                waning_rate = slopes
)

# waning rate estimates
mean(slopes$slope, na.rm = TRUE)
sd(slopes$slope, na.rm = TRUE)

# save output
saveRDS(results, file = paste0("inst/extdata/output/sim_results_", params$title, ".rds"))

# plot all sims
p <- ggplot(data = results$ve, aes(x = period, y = V, color = as.factor(Sim))) + #
  geom_line() +
  #geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.1, linetype = "dashed") + # , fill = waning_rate
  #facet_wrap(~waning_rate, nrow = 3) +
  theme(legend.position = "none")
p

# summarize results
summary_results <- results$ve %>%
  select(Sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "median", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         # lower_new = ifelse(lower < -0.5, -0.5, lower),
         # upper_new = ifelse(upper > 1.5, 1.5, upper)
         )

# plot summarized results
p <- ggplot(data = summary_results, aes(x = period, y = mean)) + # , color = waning_rate
  geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.1, linetype = "dashed") + # , fill = waning_rate
  #facet_wrap(~waning_rate, nrow = 3) +
  theme_bw()
p







