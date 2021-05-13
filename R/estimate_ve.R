# ----------------------------------------------------------
#' Apply different VE estimation methods
#'
#' This function applies teh Durham method, Tian method, and ML method to data from a vaccine
#' efficacy/effectiveness study.
#' @param dat data set
#' @param params list of input parameters
#' @param write_to_file logical. If true, outputs are written to file
#' @param path path where files are to be written. Defaults to working directory
#' @param par_tab data frame of parameters for ML method
#' @param mcmc_pars vector of MCMC inputs
#' @return list of VE estimates from each method, maximum likelihood parameter estimates from the ML method for
#' each simulation, the proportion of times the null hypothesis was rejected for each method, mean VE estimate over simulations
#' for each time period, and the mean maximum likelihood parameters over all simulations.
#' @keywords wave
#' @import survival
#' @import dplyr
#' @import tidyr
#' @export
estimate_ve <- function(dat, params, write_to_file = TRUE, path = getwd(), par_tab, mcmc_pars){

# initialise count of number of simulations in which H0 is rejected --------------------------------------
   reject_h0_durham <- reject_h0_tian <- reject_h0_ml <- 0

# loop through simulations and apply each method ---------------------------------------------------------
   for (i in 1:params$sim){
     print(i)

    # subset data for each simulation---------------------------------------------------------------------
      dat1 <- dat %>% filter(.data$Sim == i)

# method from Durham et al. 1988 -------------------------------------------------------------------------

# fit ordinary Cox propotional hazards model -------------------------------------------------------------
  flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data=dat1)

# test the proportional hazards assumption and compute the Schoenfeld residuals ($y) ---------------------
  flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
  reject_h0_durham <- reject_h0_durham + ifelse(flu_zph$table[1,3] < 0.05, 1, 0)

# calculate VE -------------------------------------------------------------------------------------------
# the nsmo argument indicates the number of time points to calculate VE at
  temp <- durham_ve(flu_zph, n_days = params$ND, n_periods = params$NJ,
                    n_days_period = params$NDJ,var = "V") %>%
    mutate(Sim = i, Method = "Durham")

# bind results together ----------------------------------------------------------------------------------
  if (i > 1){
    ve_est <- bind_rows(ve_est,temp)
  } else {ve_est <- temp}


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

# outputs ------------------------------------------------------------------------------------------------
   prop_reject_h0 <- tibble(method = c('Durham','Tian', 'ML'),
                            count = c(reject_h0_durham, reject_h0_tian,reject_h0_ml),
                            proportion = count/params$sim)
   mean_ve <- ve_est %>%
     group_by(.data$Method, .data$period) %>%
     summarise_at("ve", c(mean, sd)) %>%
     rename(ve_mean = .data$fn1, ve_sd = .data$fn2)

   mean_mle_params <- mle_param_est %>%
     filter(value == "mle") %>%
     select(-.data$Sim, -.data$Method, -.data$value) %>%
     summarise_all(.funs = "mean")

   rtn <- list(ve_est = ve_est,
               mle_param_est = mle_param_est,
               prop_reject_h0 = prop_reject_h0,
               mean_ve = mean_ve,
               mean_mle_params = mean_mle_params)

   if(write_to_file){
     # write to one excel file and have each output as separate sheet
     list_output <- c("ve_est", "mle_param_est", "prop_reject_h0", "mean_ve", "mean_mle_params")

     # create workbook
     wb <- createWorkbook(paste0("sim_output_", params$title, ".xlsx"))

     # add sheets
     addWorksheet(wb, list_output[1])
     addWorksheet(wb, list_output[2])
     addWorksheet(wb, list_output[3])
     addWorksheet(wb, list_output[4])
     addWorksheet(wb, list_output[5])

     # add outputs to their sheet
     writeData(wb, 1, ve_est)
     writeData(wb, 2, mle_param_est)
     writeData(wb, 3, prop_reject_h0)
     writeData(wb, 4, mean_ve)
     writeData(wb, 5, mean_mle_params)

     # save workbook
     saveWorkbook(wb, file = paste0("sim_output_", params$title, ".xlsx"), overwrite = TRUE)

     # write.csv(prop_reject_h0,file = paste0(path,"reject_h0_prop_",params$title,".csv"))
     # write.csv(mle_param_est, file = paste0(path,"mle_parameter_estimates_",params$title,".csv"))
     # write.csv(mean_ve,file = paste0(path,"mean_ve_estimates_",params$title,".csv"))
   }

   return(rtn)
 }
# ----------------------------------------------------------
