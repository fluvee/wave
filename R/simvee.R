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
simvee <- function(params) {
  ID <- seq(1, params$N)
  X <-  rep(0, params$N)
  V <-  rep(0, params$N)
  DINF <-  rep(0, params$N)

  subject <- data.frame(SIM = params$sim, ID=ID, X=X, V=V, DINF=DINF)

  ## SubjectY goes from d=0 to d=ND
  subjectY <- matrix(NA, nrow=params$N, ncol=params$ND+1)

  NDINF = array(0, dim=c(2, 2, params$ND))

  ### Initialize
  subjectY[,1] <- rep(0, params$N)

  # Day
  d = 0
  d_period = 0
  period = 1

  beta_d01 <- rep(params$beta_d01, each = params$NDJ)
  theta_d_wane <- params$theta_d + (params$eta * 1:params$ND)
  theta_d_wane <- ifelse(theta_d_wane > 1, 1, theta_d_wane)
  beta_d11 = beta_d01 * theta_d_wane
  beta_d00 = beta_d01 * params$phi
  beta_d10 = beta_d01 * theta_d_wane * params$phi

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
    if (V == 0 && X == 0) return(beta_d00)
    if (V == 0 && X == 1) return(beta_d01)
    if (V == 1 && X == 0) return(beta_d10)
    if (V == 1 && X == 1) return(beta_d11)
  }

  ## Iterate over number of days.
  while (d < params$ND) {
    d = d + 1
    d_period = d_period + 1

    # print(paste("Day",d, "Period", period, "Day within period", d_period))

    for (i in ID) {
      if (subjectY[i, d] == 0) {
        subjectY[i, (d+1)] = as.numeric(runif(1) < getBetaForSubject(subject[i,])[d])
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

