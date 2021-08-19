# ------------------------------------------------------------------------------
# plot VE estimates over time for simulations using MI RCT params --------------
# ------------------------------------------------------------------------------

# load required packages -------------------------------------------------------
devtools::load_all() # need to be in wave root directory!
library(dplyr)
library(ggplot2)
library(cowplot)
#library(readr)

# load outcomes files ----------------------------------------------------------
results_00 <- readRDS("inst/extdata/output/sim_results_MI_RCT_06_04_00.rds")
results_03 <- readRDS("inst/extdata/output/sim_results_MI_RCT_06_04_03.rds")
results_06 <- readRDS("inst/extdata/output/sim_results_MI_RCT_06_04_06.rds")

# summarise results files for ve plot ------------------------------------------
summary_res_00 <- results_00$ve %>%
  select(Sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "median", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         waning_rate = "0%",
         true_ve = c(rep(0.6, 7))
  )

summary_res_03 <- results_03$ve %>%
  select(Sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "median", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         waning_rate = "3%",
         true_ve = c(0.6, 0.57, 0.54, 0.51, 0.48, 0.45, 0.42)
  )

summary_res_06 <- results_06$ve %>%
  select(Sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "median", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         waning_rate = "6%",
         true_ve = c(0.6, 0.54, 0.48, 0.42, 0.36, 0.30, 0.24)
  )

ve_all <- bind_rows(summary_res_00, summary_res_03, summary_res_06)

# plot ve ----------------------------------------------------------------------
p <- ggplot(data = ve_all, aes(x = period, y = mean, color = waning_rate)) + #
  geom_line() +
  geom_line(aes(x = period, y = true_ve), color = "black", linetype = "dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = waning_rate), alpha = 0.1,
              linetype = "blank") +
  facet_wrap(~waning_rate, nrow = 1) +
  labs(y = "VE", x = "Time period") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
p

# ------------------------------------------------------------------------------

# summarise results for parameter value plot -----------------------------------
summary_wr_00 <- results_00$waning_rate %>%
  summarise_at(.vars = "slope", .funs = c("mean", "sd"), na.rm = TRUE) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         true_waning_rate = 0
  )

summary_wr_03 <- results_03$waning_rate %>%
  summarise_at(.vars = "slope", .funs = c("mean", "sd"), na.rm = TRUE) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         true_waning_rate = -0.03
  )

summary_wr_06 <- results_06$waning_rate %>%
  summarise_at(.vars = "slope", .funs = c("mean", "sd"), na.rm = TRUE) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd,
         true_waning_rate = -0.06
  )

wr_all <- bind_rows(summary_wr_00, summary_wr_03, summary_wr_06) %>%
  mutate(scenario = c("0%", "3%", "6%")) %>%
  select(-sd) %>%
  pivot_longer(cols = c(mean, true_waning_rate), names_to = "quantity",
               values_to = "value") %>%
  mutate(lower = ifelse(quantity == "true_waning_rate", NA, lower),
         upper = ifelse(quantity == "true_waning_rate", NA, upper),
         quantity = ifelse(quantity == "mean", "Mean", "True Waning Rate"))

# plot parameter values
p2 <- ggplot(data = wr_all, aes(x = scenario, y = value, color = quantity)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.05,
                              position=position_dodge(0.5)) +
  #geom_point(aes(x = waning_rate, y = true_waning_rate), color = "black") +
  labs(y = "Mean Waning Rate", x = "Scenario", color = "Quantity") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
p2

p2a <- plot_grid(p2, NULL, ncol = 2, rel_widths = c(3, 1))
p_sim <- plot_grid(p, p2a, ncol = 1, rel_heights = c(1.5, 1))
p_sim

ggsave(file = "inst/extdata/output/figure1.jpg", plot = p_sim)




