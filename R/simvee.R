#' Simulate data from a vaccine efficacy trial
#'
#' This function simulates data from an ideal randomized placebo-controlled vaccine trial. Each study participant
#' is randomly assigned a binary vaccination status V, where V=1  indicates a vaccine recipient and V=0 indicates
#' a placebo recipient, i.e., an unvaccinated person. We assume that all study participants receive the vaccine or
#' the placebo on the same calendar day, just prior to the onset of the study. The vaccine is assumed ‘leaky’, i.e.,
#' the hazard of infection of a vaccinated person is a fraction of the hazard of an unvaccinated. The parameters in
#' the basic model are as follows:
#' * lambda_dv= the probability that a participant of vaccination status V=v  who was uninfected at the end of day d-1
#'   becomes infected on day d. Then the daily hazards of infection for a vaccinee and a non-vaccinees are lambda_d1 and
#'   lambda_d0, respectively.
#' * theta_d = lambda_d1/lambda_d0  is the ratio of the hazards of a vaccinee and a non-vaccinee on day d. Then the TVE on day d
#'   is 1- theta_d.
#'
#' The simulation program iterates over days. On each day, each susceptible study participant may become infected,
#' and the probability of this event equals to her/his hazard of infection on that day. The input parameters of the
#' simulation program are the daily hazards of infection for unvaccinated persons {lambda_d0} and the hazard ratios
#' {theta_d}. The output is an ‘outcomes file’ with one record for each study participant. This record includes the
#' binary vaccination status V, and a variable DINF which gives the day on which s/he became infected. For a
#' participant who did not become infected during the study, DINF=0.
#' @param params list of input parameters
#' @param simNum number of simulations
#' @return simulated data frame and output files specified in params
#' @importFrom stats runif sd
#' @importFrom utils write.csv
#' @keywords wave
#' @export
# main simulation function
simvee <- function(params, simNum) {
  ID <- seq(1, params$N)
  X <-  rep(0, params$N)
  V <-  rep(0, params$N)
  DINF <-  rep(0, params$N)

  subject <- data.frame(ID=ID, X=X, V=V, DINF=DINF)

  ## SubjectY goes from d=0 to d=ND
  subjectY <- matrix(NA, nrow=params$N, ncol=params$ND+1)

  NDINF = array(0, dim=c(2, 2, params$ND))

  ### Initialize
  subjectY[,1] <- rep(0, params$N)

  # Day
  d = 0
  d_period = 0
  period = 1

  params$beta_d11 = params$beta_d01 * params$theta_d
  params$beta_d00 = params$beta_d01 * params$phi
  params$beta_d10 = params$beta_d01 * params$theta_d * params$phi


  ## Set value of X for each subject
  subject$X = as.numeric(runif(params$N) < params$pai)

  ## Set vaccination status for each subject
  for (i in ID) {
    if (subject[i,"X"] == 0) {
      subject[i,"V"] = as.numeric(runif(1) < params$alpha_0)
    }
    if (subject[i,"X"] == 1) {
      subject[i,"V"] = as.numeric(runif(1) < params$alpha_1)
    }
  }

  ## This function returns the beta value that should be applied to a
  #  subject based on X and V values.
  getBetaForSubject = function (sub) {
    X = sub$X
    V = sub$V
    if (V == 0 && X == 0) return(params$beta_d00)
    if (V == 0 && X == 1) return(params$beta_d01)
    if (V == 1 && X == 0) return(params$beta_d10)
    if (V == 1 && X == 1) return(params$beta_d11)
  }

  ## Iterate over number of days.
  while (d < params$ND) {
    d = d + 1
    d_period = d_period + 1

    # print(paste("Day",d, "Period", period, "Day within period", d_period))

    for (i in ID) {
      if (subjectY[i, d] == 0) {
        subjectY[i, (d+1)] = as.numeric(runif(1) < getBetaForSubject(subject[i,])[period])
        if (subjectY[i, (d+1)] == 1) {
          subject[i, "DINF"] = d
          NDINF[(subject[i,"X"]+1), (subject[i,"V"]+1), d] =
            NDINF[(subject[i,"X"]+1), (subject[i,"V"]+1), d] + 1
        }
      }
      if (subjectY[i, d] == 1) {
        subjectY[i, (d+1)] = 2
      }
      if (subjectY[i, d] == 2) {
        subjectY[i, (d+1)] = 2
      }
    }
    if(d_period == params$NDJ){
      d_period = 0
      period = period + 1
    }
  }
  return(list(subject=subject,subjectY=subjectY,NDINF=NDINF))
}

# function to run simulation and generate incidence report
run_simvee <- function(params, path = getwd()){
  for (i in 1:params$sim) {
    print(i)
    results <- simvee(params, i)

    if (params$population_report_file == TRUE) {
      temp_population_report <- cbind("Sim"=i, results$subject)
      if (i > 1) {
        population_report <- rbind(population_report, temp_population_report)
      }
      else {
        population_report <- temp_population_report
      }
    }

    if (params$detailed_file == TRUE) {
      temp_detailed <- cbind(i, results$subject, results$subjectY)
      colnames(temp_detailed) <- c("sim", names(results$subject), paste0("D",seq(0,params$ND)))
      if (i > 1) {
        detailed = rbind(detailed, temp_detailed)
      }
      else {
        detailed = temp_detailed
      }
    }

    day_period <- rep(seq(1:params$NJ), each=params$NDJ)
    temp_daily <- cbind(i, seq(1:params$ND), day_period, t(apply(results$NDINF, 3L, c)))
    colnames(temp_daily) <- c("Sim", "Day", "Period", "X0V0", "X1V0", "X0V1", "X1V1")
    if (i > 1) {
      incidence = rbind(incidence, temp_daily)
    }
    else {
      incidence = temp_daily
    }
  }
  # Tibble for dplyr
  g_inc <- as_tibble(incidence)

  # write simulation results to different output files
  if (params$csv == TRUE) {
    if (params$population_report_file == TRUE) {
      write.csv(population_report, paste0(path,'/Outcomes_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$detailed_file == TRUE) {
      write.csv(detailed, paste0(path,'/Detailed_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$daily_each_sim_file == TRUE) {
      write.csv(incidence, paste0(path,'/Daily_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$daily_overall_file == TRUE) {
      inc_daily_overall <- g_inc %>% group_by(.data$Day,.data$Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(-.data$Sim)
      write.csv(inc_daily_overall, paste0(path,'/Daily_overall_',params$title,'.csv'), row.names = FALSE)
    }
    # if (params$period_each_sim_file == TRUE) {
    #   inc_period_sim <- g_inc %>% group_by(Sim,Period) %>% summarise_all(sum, na.rm = TRUE) %>%
    #                           select(-Day) %>% ungroup() %>% group_by(Period) %>%
    #                           summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
    #   write.csv(inc_period_sim, paste0(path,'Period_',params$title,'.csv'), row.names = FALSE)
    # }
    # if (params$period_overall_file == TRUE) {
    #   inc_period_overall <- g_inc %>% group_by(Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(c(-Day,-Sim))
    #   write.csv(inc_period_overall, paste0(path,'Period_overall_',params$title,'.csv'), row.names = FALSE)
    # }
    if (params$seasonal_each_sim_file == TRUE) {
      inc_seasonal_sim <- g_inc %>% group_by(.data$Sim) %>% summarise_all(sum, na.rm = TRUE) %>% select(c(-.data$Day,-.data$Period))
      write.csv(inc_seasonal_sim, paste0(path,'/Seasonal_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$seasonal_overall_file == TRUE) {
      inc_seasonal_overall <- g_inc %>% group_by(.data$Sim) %>% summarise_all(sum, na.rm = TRUE) %>%
                              select(c(-.data$Day,-.data$Period)) %>%  summarize_all(mean, na.rm=TRUE) %>%
                              select(c(-.data$Sim))
      write.csv(inc_seasonal_overall, paste0(path,'/Seasonal_overall_',params$title,'.csv'), row.names = FALSE)
    }
  }
  # return incidence tibble
  return(population_report)
}
#
# if (params$sas == TRUE) {
#   library(haven)
#   if (params$population_report_file == TRUE) {
#     write_sas(population_report, paste0('Outcomes_',params$title,'.sas7bdat'))
#   }
#   if (params$detailed_file == TRUE) {
#     write_sas(detailed, paste0('Detailed_',params$title,'.sas7bdat'))
#   }
#   if (params$daily_each_sim_file == TRUE) {
#     write_sas(incidence, paste0('Daily_',params$title,'.sas7bdat'))
#   }
#   if (params$daily_overall_file == TRUE) {
#     inc_daily_overall <- g_inc %>% group_by(Day,Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
#     write_sas(data.finc_daily_overall, paste0('Daily_overall_',params$title,'.sas7bdat'))
#   }
#   if (params$period_each_sim_file == TRUE) {
#     inc_period_sim <- g_inc %>% group_by(Sim,Period) %>% summarise_all(sum, na.rm = TRUE) %>%
#       select(-Day) %>% ungroup() %>% group_by(Period) %>%
#       summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
#     write_sas(inc_period_sim, paste0('Period_',params$title,'.sas7bdat'))
#   }
#   if (params$period_overall_file == TRUE) {
#     inc_period_overall <- g_inc %>% group_by(Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(c(-Day,-Sim))
#     write_sas(inc_period_overall, paste0('Period_overall_',params$title,'.sas7bdat'))
#   }
#   if (params$seasonal_each_sim_file == TRUE) {
#     inc_seasonal_sim <- g_inc %>% group_by(Sim) %>% summarise_all(sum, na.rm = TRUE) %>% select(c(-Day,-Period))
#     write_sas(inc_seasonal_sim, paste0('Seasonal_',params$title,'.sas7bdat'))
#   }
#   if (params$seasonal_overall_file == TRUE) {
#     inc_seasonal_overall <- g_inc %>% group_by(Sim) %>% summarise_all(sum, na.rm = TRUE) %>%
#       select(c(-Day,-Period)) %>%  summarize_all(mean, na.rm=TRUE) %>%
#       select(c(-Sim))
#     write_sas(inc_seasonal_overall, paste0('Seasonal_overall_',params$title,'.sas7bdat'))
#   }
# }



