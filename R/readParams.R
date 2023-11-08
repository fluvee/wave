#' Read in parameters
#'
#' Function reads in parameters for simulation from a csv file
#' @param filename path to csv file
#' @return list of parameters
#' @keywords wave
#' @export
readParams = function(filename = NULL) {
  if(is.null(filename))(filename = file.choose())
  para <- list()
  Y <- readLines(filename, warn = FALSE)
  separater <- substr(Y[[1]], nchar(Y[[1]])-1, nchar(Y[[1]])-1)

  para$title <- strsplit(Y[grep("Title",Y)],separater)[[1]][2]
  para$sim <- as.numeric(strsplit(Y[grep("Number of simulations",Y)],separater)[[1]][2])
  para$rseed <- as.numeric(strsplit(Y[grep("Seed",Y)],separater)[[1]][2])
  para$all_or_no_vaccine <- as.numeric(strsplit(Y[grep("Leaky or all-or-none vaccine",Y)],separater)[[1]][2])

  para$sas <- as.character(strsplit(Y[grep("Output files in sas format",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$csv <- as.character(strsplit(Y[grep("Output files in csv format",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$inputsfile <- ifelse(strsplit(Y[grep("Input and calculated parameters",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$detailed_file <- ifelse(strsplit(Y[grep("Detailed for each simulation",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$daily_each_sim_file <- ifelse(strsplit(Y[grep("Incidence-daily for each simulation",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
# para$period_each_sim_file <- ifelse(strsplit(Y[grep("Incidence-period for each simulation",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$seasonal_each_sim_file <- ifelse(strsplit(Y[grep("Incidence-total for each simulation",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$daily_overall_file <- ifelse(strsplit(Y[grep("Incidence-daily means over all simulations",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
# para$period_overall_file <- ifelse(strsplit(Y[grep("Incidence-period means over all simulations",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$seasonal_overall_file <- ifelse(strsplit(Y[grep("Incidence-total means over all simulations",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$population_report_file <- ifelse(strsplit(Y[grep("Outcomes for each simulation",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)
  para$timestamp <- ifelse(strsplit(Y[grep("Add timestamp to output file names",Y)],separater)[[1]][2]=="yes",TRUE,FALSE)

  para$NDJ <- as.numeric(strsplit(Y[grep("Number of days in each period",Y)],separater)[[1]][2])
  para$NJ <- as.numeric(strsplit(Y[grep("Number of periods in the study",Y)],separater)[[1]][2])
  para$ND <- as.numeric(strsplit(Y[grep("Number of days in the study",Y)],separater)[[1]][2])
  if (para$ND != para$NDJ * para$NJ) {
    stop("Number of days in the study isn't a multiple of number of periods")
  }
  para$N <- as.numeric(strsplit(Y[grep("Number of participants in study",Y)],separater)[[1]][2])
  para$pai <- as.numeric(strsplit(Y[grep("Probability of X=1",Y)],separater)[[1]][2])
  para$alpha_1 <- as.numeric(strsplit(Y[grep("Probability of vaccination for X=1",Y)],separater)[[1]][2])
  para$alpha_0 <- as.numeric(strsplit(Y[grep("Probability of vaccination for X=0",Y)],separater)[[1]][2])
  if (identical(Y[grep("Daily probability of infection for V=0 X=1",Y)], character(0))) {
    for (i in 1:para$NJ){
      para$beta_d01[i] <- as.numeric(strsplit(Y[grep("Daily probabilities of infection when V=0 X=1",Y)+i],separater)[[1]][3])
    }
  } else {
    para$beta_d01 <- as.numeric(strsplit(Y[grep("Daily probability of infection for V=0 X=1",Y)],separater)[[1]][2])
  }
  if (identical(Y[grep("Multiplier for daily probability of infection for V=1",Y)], character(0))) {
    for (i in 1:para$NJ){
      para$theta_d[i] <- as.numeric(strsplit(Y[grep("Daily multiplier for daily probability of infection when V=1 for each period",Y)+i],separater)[[1]][3])
    }
  } else {
    para$theta_d <- as.numeric(strsplit(Y[grep("Multiplier for daily probability of infection for V=1",Y)],separater)[[1]][2])
  }
  para$phi <- as.numeric(strsplit(Y[grep("Multiplier for daily probability of infection when X=0",Y)],separater)[[1]][2])
  para$eta <- as.numeric(strsplit(Y[grep("Waning rate",Y)],separater)[[1]][2])

  return(para)
}
