#' Function to run simvee simulation function and generate incidence report
#' @param params list of input parameters
#' @param path directory where output files are to be saved. Defaults to working directory.
#' @return output files specified in params
#' @importFrom tibble as_tibble
#' @importFrom utils write.csv
#' @import tidyr
#' @import dplyr
#' @keywords wave
#' @export

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
      inc_daily_overall <- g_inc %>% group_by(.data$Day,.data$Period) %>% summarise_all(mean, na.rm = TRUE) %>% dplyr::select(-.data$Sim)
      write.csv(inc_daily_overall, paste0(path,'/Daily_overall_',params$title,'.csv'), row.names = FALSE)
    }
    # if (params$period_each_sim_file == TRUE) {
    #   inc_period_sim <- g_inc %>% group_by(Sim,Period) %>% summarise_all(sum, na.rm = TRUE) %>%
    #                           dplyr::select(-Day) %>% ungroup() %>% group_by(Period) %>%
    #                           summarise_all(mean, na.rm = TRUE) %>% dplyr::select(-Sim)
    #   write.csv(inc_period_sim, paste0(path,'Period_',params$title,'.csv'), row.names = FALSE)
    # }
    # if (params$period_overall_file == TRUE) {
    #   inc_period_overall <- g_inc %>% group_by(Period) %>% summarise_all(mean, na.rm = TRUE) %>% dplyr::select(c(-Day,-Sim))
    #   write.csv(inc_period_overall, paste0(path,'Period_overall_',params$title,'.csv'), row.names = FALSE)
    # }
    if (params$seasonal_each_sim_file == TRUE) {
      inc_seasonal_sim <- g_inc %>% group_by(.data$Sim) %>% summarise_all(sum, na.rm = TRUE) %>% dplyr::select(c(-.data$Day,-.data$Period))
      write.csv(inc_seasonal_sim, paste0(path,'/Seasonal_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$seasonal_overall_file == TRUE) {
      inc_seasonal_overall <- g_inc %>% group_by(.data$Sim) %>% summarise_all(sum, na.rm = TRUE) %>%
        dplyr::select(c(-.data$Day,-.data$Period)) %>%  summarize_all(mean, na.rm=TRUE) %>%
        dplyr::select(c(-.data$Sim))
      write.csv(inc_seasonal_overall, paste0(path,'/Seasonal_overall_',params$title,'.csv'), row.names = FALSE)
    }
  }
  # return incidence tibble
  return(population_report)
}

# for SAS output
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
#     inc_daily_overall <- g_inc %>% group_by(Day,Period) %>% summarise_all(mean, na.rm = TRUE) %>% dplyr::select(-Sim)
#     write_sas(data.finc_daily_overall, paste0('Daily_overall_',params$title,'.sas7bdat'))
#   }
#   if (params$period_each_sim_file == TRUE) {
#     inc_period_sim <- g_inc %>% group_by(Sim,Period) %>% summarise_all(sum, na.rm = TRUE) %>%
#       dplyr::select(-Day) %>% ungroup() %>% group_by(Period) %>%
#       summarise_all(mean, na.rm = TRUE) %>% dplyr::select(-Sim)
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

