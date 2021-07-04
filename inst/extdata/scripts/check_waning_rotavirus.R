# --------------------------------------------------
# Check waning rate using rotavirus input files
# --------------------------------------------------

# load required packages
library(readr)
library(dplyr)
library(ggplot2)
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
    #print(i)
    # subset data for current simulation (i)
    dat <- outcomes_dat1 %>%
      filter(Sim == i)
    # fit ordinary Cox propotional hazards model
    flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data = dat)

    # test the proportional hazards assumption and compute the Schoenfeld residuals ($y)
    flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
    reject_h0_durham <- reject_h0_durham + ifelse(flu_zph$table[1,3] < 0.05, 1, 0)
    x <- flu_zph
    df = 4
    # calculate VE
    xx <- x$x
    yy <- x$y
    d <- nrow(yy)
    df <- max(df)
    nvar <- ncol(yy)
    pred.x <- seq(from = (params$NDJ/2), to = params$ND - (params$NDJ/2), length = params$NJ)
    temp <- c(pred.x, xx)
    lmat <- ns(temp, df = df, intercept = TRUE)
    pmat <- lmat[1:length(pred.x), ]
    xmat <- lmat[-(1:length(pred.x)), ]
    qmat <- qr(xmat)
    if (x$transform!="identity")
      stop("please re-fit the Cox model with the identity transform")
    if (qmat$rank < df)
      stop("Spline fit is singular, try a smaller degrees of freedom")
    # se
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- ((pmat %*% xtx) * pmat) %*% rep(1, df) # supposed to be multiplied by factor of d, but that causes MASSIVE errors

    # ve estimate
    yhat.beta <- pmat %*% qr.coef(qmat, yy)
    yhat <- 1-exp(yhat.beta) # VE estimate

    # se
    tmp <- 2 * sqrt(x$var[1,1] * seval)
    ylow.beta <- yhat.beta - tmp
    yup.beta <- yhat.beta + tmp
    ylow <- 1-exp(yup.beta)
    yup <- 1-exp(ylow.beta)
    yr <- range(yhat, yup, ylow)

    output <- data.frame(time = pred.x,
                         period = 1:params$NJ,
                         yhat = yhat,
                         ylow = ylow,
                         yup = yup,
                         sim = i)
    # store output in a list
    rtn[[i]] <- output
      # durham_ve(flu_zph, df = 4, n_days = params$ND, n_periods = params$NJ,
      #                     n_days_period = params$NDJ, var = "V") %>%
      # mutate(Sim = i, Method = "Durham")
  }
  print(reject_h0_durham/params$sim)
  if (j == 1){ rtn1 <- rtn
  } else if (j == 2) { rtn2 <- rtn
  } else { rtn3 <- rtn }
}
# bind the outputs of each simulation
results0 <- bind_rows(rtn1) %>%
  select(sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd)

results03 <- bind_rows(rtn2) %>%
  select(sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd)

results06 <- bind_rows(rtn3) %>%
  select(sim, period, V) %>%
  group_by(period) %>%
  summarise_at(.vars = "V", .funs = c("mean", "sd")) %>%
  mutate(lower = mean - 1.96 * sd,
         upper = mean + 1.96 * sd)

results_all <- bind_rows(results0,
                         results03,
                         results06,
                         .id = "waning_rate") %>%
  mutate(waning_rate = case_when(
    waning_rate == 1 ~ "no waning",
    waning_rate == 2 ~ "3% waning",
    waning_rate == 3 ~ "6% waning"
  ))

# fancier ggplot
p <- ggplot(data = results_all, aes(x = period, y = mean, color = waning_rate)) +
  geom_line() +
  geom_ribbon(aes(ymin=lower,ymax=upper, fill = waning_rate), alpha = 0.1, linetype = "dashed") +
  facet_wrap(~waning_rate, nrow = 3) +
  theme_bw()
p
