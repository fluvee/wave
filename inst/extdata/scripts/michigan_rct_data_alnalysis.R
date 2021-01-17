# Michigan clinical trial analysis script

# load necessary packages
library(readr)
library(tidyr)
library(dplyr)
library(survival)
library(timereg)
library(lazymcmc)
library(foreign)
library(splines)

# load data
mi_data <- read_csv("inst/extdata/data/Michigan_0708.csv")

# data wrangling
mi_data1 <- mi_data %>%
  mutate(V = ifelse(vaccine_type != "PLACEBO" , 1, 0),
         FARI = ifelse(is.na(onset_date),0,1),
         DINF = onset_date - min(mi_data$vaccination_date),
         DINF_new = ifelse(is.na(DINF), 999, DINF)) %>%
  select(studyid, V, DINF_new, FARI)


# method from Durham et al. 1988 -------------------------------------------------------------------------

# fit ordinary Cox propotional hazards model -------------------------------------------------------------
flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data = mi_data1)

# test the proportional hazards assumption and compute the Schoenfeld residuals ($y) ---------------------
flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
reject_h0_durham <- ifelse(flu_zph$table[1,3] < 0.05, 1, 0)

# calculate VE -------------------------------------------------------------------------------------------
# the nsmo argument indicates the number of time points to calculate VE at
durham_est <- durham_ve(flu_zph, n_days = 197, n_periods = 7,
                  n_days_period = 30, var = "V")

# method from Tian et al. 2005 ---------------------------------------------------------------------------

# calculate VE
# n_timepoint_breaks argument specifies the number of time points to calculate VE for
temp2 <- tian_ve(dat1, n_days = params$ND, n_periods = params$NJ,
                 n_days_period = params$NDJ)
#temp2a <- temp2$output %>% mutate(Sim = i, Method = "Tian")
#ve_est <- bind_rows(ve_est,temp2a)
# proportion of sims where null hypothesis is rejected
reject_h0_tian <- reject_h0_tian + temp2 #temp2$reject_h0


# method from Ainslie et al. 2017 ------------------------------------------------------------------------

temp3 <- ml_ve2(dat1, params, par_tab = par_tab, mcmc_pars = mcmc_pars, file_name = params$title)
temp3a <- temp3$ve_dat %>% mutate(Sim = i, Method = "ML")
ve_est <- bind_rows(ve_est,temp3a)

# proportion of sims where null hypothesis is rejected
reject_h0_ml <- reject_h0_ml + ifelse(temp3$param_est$lambda[2] > 0, 1, 0)
