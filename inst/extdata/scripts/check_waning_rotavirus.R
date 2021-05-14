# --------------------------------------------------
# Check waning rate using rotavirus input files
# --------------------------------------------------

# load required packages
library(readr)
library(dplyr)
devtools::load_all()

# read in one parameter to get number of sims
params <- readParams("./inst/extdata/input/Rotavirus_Input_ban_406.csv")

# read in outcomes files
no_waning <- read_csv("inst/extdata/output/Rotavirus_Outcomes_ban_400.csv") %>%
  mutate(FARI = ifelse(DINF == 0, 0, 1),
         DINF_new = ifelse(DINF == 0, 999, DINF))
waning_03 <- read_csv("inst/extdata/output/Rotavirus_Outcomes_ban_403.csv") %>%
  mutate(FARI = ifelse(DINF == 0, 0, 1),
         DINF_new = ifelse(DINF == 0, 999, DINF))
waning_06 <- read_csv("inst/extdata/output/Rotavirus_Outcomes_ban_406.csv") %>%
  mutate(FARI = ifelse(DINF == 0, 0, 1),
         DINF_new = ifelse(DINF == 0, 999, DINF))

# use Durham method to calculate VE over time

# loop over data sets
for (j in 1:3){
  print(j)
  reject_h0_durham <- 0 # set h0 rejection counter to 0
  rtn <- list()         # empty list for output

  if (j == 1){ outcomes_dat1 <- no_waning
  } else if (j == 2) { outcomes_dat1 <- waning_03
  } else { outcomes_dat1 <- waning_06 }

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
  print(reject_h0_durham/params$sim)
  if (j == 1){ rtn1 <- rtn
  } else if (j == 2) { rtn2 <- rtn
  } else { rtn3 <- rtn }
}
# bind the outputs of each simulation
results0 <- bind_rows(rtn1) %>%
  group_by(Sim, period)

results03 <- bind_rows(rtn2) %>%
  group_by(Sim, period)

results06 <- bind_rows(rtn3) %>%
  group_by(Sim, period)
