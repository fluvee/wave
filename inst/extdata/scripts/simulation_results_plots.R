# ------------------------------------------------------------------------------
# plot VE estimates over time for simulations using MI RCT params --------------
# ------------------------------------------------------------------------------

# load required packages
devtools::load_all() # need to be in wave root directory!
library(dplyr)
library(ggplot2)
#library(readr)

# load outcomes files
results_00 <- readRDS("inst/extdata/output/sim_results_MI_RCT_06_04_00.rds") %>%
  mutate(waning_rate = "0%")
results_03 <- readRDS("inst/extdata/output/sim_results_MI_RCT_06_04_03.rds") %>%
  mutate(waning_rate = "3%")
results_06 <- readRDS("inst/extdata/output/sim_results_MI_RCT_06_04_06.rds") %>%
  mutate(waning_rate = "6%")

