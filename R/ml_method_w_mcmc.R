# ----------------------------------------------------------
#' Estimate VE using a maximum likelihood (ML) method with MCMC
#'
#' This method is similar to that of Ainslie et al. 2017 (SIM). The ML method is based on calculating the
#' contribution to the likelihood function of each study participant. See XX for more details.
#' @param data data set
#' @param params input parameters
#' @param my_prior prior function for MCMC
#' @param file_name character string of output file name
#' @param par_tab data frame of parameter values
#' @param mcmc_pars vector of MCMC options
#' @return list with a tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period) and the maximu likelihood estimates of the parameters.
#' @keywords wave
#' @import dplyr
#' @import tidyr
#' @importFrom utils read.csv
#' @importFrom stats quantile
#' @export
ml_ve2 <- function(data, params, my_prior = NULL, file_name = "test", par_tab, mcmc_pars){

## The MCMC code uses the parameter table. Here, we should specify some random starting
## points in the "values" column.
startTab <- par_tab
## You could have done something like this:
startTab$values[1:3] <- runif(3, startTab$lower_bound[1:3], startTab$upper_bound[1:3])

output <- run_MCMC(parTab = startTab, data = data, mcmcPars = mcmc_pars,
                   filename = file_name,
                   CREATE_POSTERIOR_FUNC = my_creation_function, mvrPars = NULL,
                   PRIOR_FUNC = my_prior, OPT_TUNING = 0.2)

# plot results (exclude adaptive period)
chain <- read.csv(output$file)
# plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))

## Determine mle
chain1 <- chain[chain$sampno > mcmcPars["adaptive_period"],]
max_loglik <- which(chain1$lnlike == max(chain1$lnlike))
mles <- unique(chain1[max_loglik, c("alpha", "theta_0", "phi")])
mles[2:3, ] <- apply(chain1[,2:4], 2, function(x) quantile(x, probs = c(0.025,0.975)))
rownames(mles) <- c("mle", "2.5%", "97.5%")
mles$lambda <- mles$phi - 1

periods <- rep(1:params$NJ, each = params$NDJ)
ve_dat <- tibble(day = 1:params$ND, period = periods, ve = 1 - (mles$theta_0 + (mles$lambda * .data$day))) %>%
  select(-.data$day) %>%
  group_by(.data$period) %>%
  summarise_all(.funs = mean) %>%
  mutate(ve = ifelse(.data$ve > 1, 1, .data$ve))

# output
rtn <- list(param_est = mles, ve_dat = ve_dat)

return(rtn)
}
