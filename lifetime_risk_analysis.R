###############################################################################
#  Code for estimating the lifetime risk of disease in the Swedish population
#
#  Purpose:
#
#    Estimate the cumulative lifetime risk (from index age 18 to age 95) of a given
#    diagnosis in the Swedish adult population, overall and by sex, while
#    accounting for
#
#      (i)  left truncation: 
#               individuals enter the risk set at age 18 or 
#               at their age at the start of follow-up, whichever is later and
#
#      (ii) the competing risk of death: 
#             death without the diagnosis is treated as a competing event.
#
#  Methods:
#
#    1. Cumulative incidence functions (CIFs) are estimated with the
#       Aalen-Johansen estimator (survival::survfit with a multi-state
#       outcome: censored / event / death), using age as the time scale.
#
#    2. Pseudo-observations of the CIF at age 95 are computed for every
#       individual (survival::pseudo).
#
#    3. The pseudo-observations are regressed on sex using a generalised
#       linear model with a complementary log-log link and robust
#       (sandwich) standard errors (survey::svyglm, one cluster per person).
#
#    4. Lifetime risks, the female-male difference, and 95% confidence
#       intervals are derived from the model using the delta method
#       (msm::deltamethod).
#
###############################################################################


###############################################################################
## Load Packages
###############################################################################

library(arrow)      
library(dplyr)      
library(lubridate)  
library(zoo)        
library(survival)   
library(survey)     
library(msm)        

set.seed(123)       


###############################################################################
## Configuration options
###############################################################################

# Paths in the TRE environment
paths <- list(
  disease_dir = "/safe/data/APNC/DATA/ANALYSIS DATA/",
  tpr         = "/safe/data/APNC/DATA/TPR_45.parquet",
  cdr         = "/safe/data/SCARP/DATA_2024/PARQUET/ut_r_dors_19328_2023.parquet",
  cdr_avi     = "/safe/data/SCARP/DATA_2024/PARQUET/ut_r_dors_avi_19328_2023.parquet",
  migration   = "/safe/data/SCARP/DATA_2024/PARQUET/ar_lev_migrationer.parquet",
  out_dir     = "/safe/data/APNC/SCARP-Lifecourse/Data/Derived/ltr/Diagnoses/"
)

END_OF_FOLLOWUP_YEAR <- 2024 # Last calendar year of follow-up



################################################################################
## Help functions
################################################################################

#' Bounded complementary log-log link
#' Pseudo-observations can fall outside (0, 1), where the standard cloglog
#' link function is undefined. This link keeps the usual cloglog inverse
#' link, mu = 1 - exp(-exp(eta)), but truncates mu to [edge, 1 - edge]
#' before applying the link function. glm() only uses the link function to
#' compute starting values, so the fitted estimates are unaffected.
#' 
#' @param edge Truncation bound. With the default (0.5), every observation
#'   receives the same starting value, eta = log(-log(0.5)).
#' @return A "link-glm" object for use in a family, e.g. gaussian(link = .)

bcloglog <- function(edge = 0.5) {
  link <- make.link("cloglog")
  link$linkfun <- function(mu) {
    x <- pmax(edge, pmin(mu, 1 - edge))
    log(-log(1 - x))
  }
  link$name <- "bcloglog"
  link
}

################################################################################
#' Read a Cause of Death Register file
#' 
#' @param path Path to the parquet file.
#' @return Data frame with lopnr, Death_Year and Death_dt (year-month).
read_cdr <- function(path) {
  open_dataset(path) %>%
    collect() %>%
    mutate(Death_Year = as.numeric(substr(DODSDAT, 1, 4)),
           Death_dt   = as.yearmon(DODSDAT, "%Y%m")) %>%
    select(lopnr = LopNr, Death_Year, Death_dt)
}

################################################################################
#' Format an estimate and its Wald 95% CI as percentages
#' 
#' @param est Point estimate on the probability scale.
#' @param se  Standard error on the probability scale.
#' @return Character string "est (lower - upper)".
fmt_ci <- function(est, se) {
  sprintf("%.2f (%.2f - %.2f)",
          est * 100, (est - 1.96 * se) * 100, (est + 1.96 * se) * 100)
}

################################################################################
#' Youngest age at which the cumulative incidence exceeds a threshold
#'
#' @param df        Data frame with columns time (age) and cif.
#' @param threshold Cumulative incidence threshold (e.g. 0.01 for 1%).
#' @return Age, or NA if the threshold is never reached.
first_age_above <- function(df, threshold) {
  ages <- df$time[df$cif > threshold]
  if (length(ages) == 0) NA_real_ else min(ages)
}


################################################################################
##  Main analysis function
################################################################################

#' Estimate the lifetime risk of a diagnosis, overall and by sex
#'
#' @param Disease    Name of the diagnosis (must match the input file name).
#' @param index_age  Age at study entry.
#' @param exit_age   Age at which lifetime risk is evaluated (95 years)
#' @param saving     If TRUE, write results and CIF curves to CSV.
#' @param start_year First calendar year of follow-up.
#' 
#' @return One-row data frame with lifetime risks (in %, with 95% CI), the
#'   female-male difference, its p-value, and the ages at which the
#'   cumulative incidence first exceeds 1% and 10%.

calculate_lifetime_risk <- function(Disease    = "MI",
                                    index_age  = 18,
                                    exit_age   = 95,
                                    saving     = TRUE,
                                    start_year = 2014) {

  message("Analysing ", Disease)

  ################################################################################
  ##  Read the data
  ################################################################################
  
  ## Diagnosis data: 
  
   file_path <- paste0(paths$disease_dir, "ICD_",
                      paste(Disease, "all_unique", "2", "1987", "2024",
                            sep = "_"),
                      ".parquet")
   
  if (!file.exists(file_path)) stop("Disease file not found: ", file_path)
  disease_data <- read_parquet(file_path)

  ## Total Population Register:
  
  message("Loading TPR data from ", paths$tpr)
  tpr <- read_parquet(paths$tpr) %>%
    filter(if_any(num_range("InSweden", start_year:END_OF_FOLLOWUP_YEAR),
                  ~ .x == 1))
  lopnr_list <- tpr$lopnr

  ## Cause of Death Register: 
  
  message("Loading CDR data")
  
  cdr_main <- read_cdr(paths$cdr)
  cdr_avi  <- read_cdr(paths$cdr_avi) %>%
    anti_join(cdr_main, by = "lopnr")
  cdr      <- bind_rows(cdr_main, cdr_avi)
  rm(cdr_main, cdr_avi)
  gc()

  ## Emigration data: 
  
  message("Loading migration data from ", paths$migration)
  
  emig <- open_dataset(paths$migration) %>%
    collect() %>%
    filter(posttyp == "Utv") %>%                     # emigrations only
    mutate(emig_date = as.Date(datum, "%Y%m%d"),
           emig_year = year(emig_date)) %>%
    select(lopnr, emig_date, emig_year) %>%
    left_join(tpr, by = "lopnr", relationship = "many-to-many") %>%
    mutate(age_start = start_year - as.numeric(Fodelsear),
           emig_age  = emig_year  - as.numeric(Fodelsear),
           emigrated = TRUE,
           
           # entry age = index_age, or age at start_year if older
           entry_age = if_else(age_start > index_age, age_start, index_age)) %>%
    filter(emig_age >= entry_age) %>%               # emigration after entry
    arrange(lopnr, emig_date) %>%
    group_by(lopnr) %>%
    slice(1) %>%                                    # first such emigration
    ungroup() %>%
    select(lopnr, emig_date, emig_year, emig_age, emigrated)

  ################################################################################
  ## Build the analysis data set
  ## status: 
  ##         0 = censored, 
  ##         1 = diagnosis
  ##         2 = death without diagnosis
  ################################################################################
  
  ## Cases: 
  
  cases <- disease_data %>%
    filter(lopnr %in% lopnr_list) %>%
    left_join(cdr, by = "lopnr") %>%
    mutate(
      # Age at death, or age at end of follow-up if alive (capped at exit_age)
      Death_Age = if_else(is.na(Death_Year),
                          pmin(END_OF_FOLLOWUP_YEAR - B_Year, exit_age),
                          pmin(Death_Year - B_Year, exit_age)),
      ev_Death  = if_else(is.na(Death_Year) | Death_Age >= exit_age, 0, 1),
      status    = 1,
      age_start = start_year - B_Year,
      entry_age = if_else(age_start > index_age, age_start, index_age)
    ) %>%
    select(lopnr, AR, B_Year, status, alder, Kon, ev_Death, entry_age) %>%
    distinct()

  ## Non-cases: 
  
  noncases <- tpr %>%
    anti_join(cases, by = "lopnr") %>%
    mutate(AR     = if_else(Doddatum == "",
                            END_OF_FOLLOWUP_YEAR,
                            as.numeric(substring(Doddatum, 1, 4))),
           B_Year = as.numeric(Fodelsear)) %>%
    mutate(alder    = AR - B_Year,
           ev_Death = if_else(Doddatum == "" | alder >= exit_age, 0, 1)) %>%
    mutate(status = 2 * ev_Death) %>%           # 2 = death, 0 = censored
    left_join(emig, by = "lopnr") %>%
    mutate(
      # Censor people who emigrated after entry at their emigration age
      alder     = if_else(!is.na(emigrated), emig_age, alder),
      status    = if_else(!is.na(emigrated), 0, status),
      age_start = start_year - B_Year,
      entry_age = if_else(age_start > index_age, age_start, index_age)
    ) %>%
    select(lopnr, AR, B_Year, status, alder, Kon, ev_Death, entry_age)

  ## Combine and apply the inclusion criteria:
  ##  - exit at or after index_age 
  ##  - exit in or after start_year (people diagnosed before start_year are
  ##    prevalent cases, are not at risk, and are excluded)
  ##  - non-missing sex
  
  analysis_data <- bind_rows(cases, noncases) %>%
    filter(alder >= index_age,
           AR    >= start_year,
           Kon   != "") %>%
    
  # Ages are in whole years. 
  # survfit() requires exit > entry, so people
  # who enter and exit at the same age are given one year of follow-up.

    mutate(alder = if_else(entry_age == alder, entry_age + 1, alder)) %>%
    mutate(status = factor(status, levels = c(0, 1, 2),
                           labels = c("Censored", "Event", "Death")),
           Kon    = relevel(factor(Kon, labels = c("Male", "Female")),
                            ref = "Male")) %>%
    distinct()
  rm(cases, noncases, disease_data, tpr, cdr, emig)
  gc()
  message("Combined dataset created: ", nrow(analysis_data), " individuals")

  ################################################################################
  ##  Cumulative incidence (Aalen-Johansen), age as time scale
  ################################################################################


  fit_overall <- survfit(Surv(entry_age, alder, status) ~ 1,
                         data = analysis_data, id = lopnr)
  message("Overall CIF calculated")

  fit_sex <- survfit(Surv(entry_age, alder, status) ~ Kon,
                     data = analysis_data, id = lopnr)
  message("CIF stratified by sex calculated")

  ################################################################################
  ## Pseudo-observation regression for the CIF at exit_age
  ################################################################################
  
  pseudo_obs <- pseudo(fit_overall, times = exit_age, data.frame = TRUE)
  message("Pseudo-observations created")

  pdata <- merge(analysis_data,
                 subset(pseudo_obs, state == "Event"),
                 by = "lopnr") %>%
    select(lopnr, status, entry_age, alder, time, state, pseudo, Kon)

  # One cluster per person gives robust (sandwich) standard errors
  mdesign <- svydesign(data = pdata, 
                       id = ~lopnr,
                       weights = NULL,
                       variables = ~ . - lopnr)

  message("Fitting models ...")
  fit_by_sex  <- svyglm(pseudo ~ Kon,
                        design = mdesign,
                        family = gaussian(link = bcloglog()))
  
  fit_overall_glm <- svyglm(pseudo ~ 1, 
                            design = mdesign,
                            family = gaussian(link = bcloglog()))

  ################################################################################
  ## Lifetime risks with delta-method confidence intervals
  ## With a cloglog link, risk = 1 - exp(-exp(linear predictor)).
  ################################################################################
  
  message("Calculating results")

  # By sex: 
  #     b0 = intercept (males),
  #     b1 = female vs male
  
  coefs <- coef(fit_by_sex)
  b0    <- unname(coefs[1])
  b1    <- unname(coefs[2])
  vbeta <- vcov(fit_by_sex)

  ltr_male   <- 1 - exp(-exp(b0))
  ltr_female <- 1 - exp(-exp(b0 + b1))
  ltr_diff   <- ltr_female - ltr_male

  se_male   <- deltamethod(~ 1 - exp(-exp(x1)),
                           c(b0, b1), vbeta)
  se_female <- deltamethod(~ 1 - exp(-exp(x1 + x2)),
                           c(b0, b1), vbeta)
  se_diff   <- deltamethod(~ (1 - exp(-exp(x1 + x2))) - (1 - exp(-exp(x1))),
                           c(b0, b1), vbeta)

  # Overall
  b0_over  <- unname(coef(fit_overall_glm)[1])
  ltr_over <- 1 - exp(-exp(b0_over))
  se_over  <- deltamethod(~ 1 - exp(-exp(x1)), b0_over, vcov(fit_overall_glm))

  ################################################################################
  ##  CIF curves for plotting
  ################################################################################

  event_col_overall <- which(fit_overall$states == "Event")
  event_col_sex     <- which(fit_sex$states == "Event")
  sex_labels <- rep(sub("^Kon=", "", names(fit_sex$strata)), fit_sex$strata)

  plot_data <- bind_rows(
    data.frame(time  = fit_overall$time,
               cif   = fit_overall$pstate[, event_col_overall],
               group = "Overall"),
    data.frame(time  = fit_sex$time,
               cif   = fit_sex$pstate[, event_col_sex],
               group = sex_labels)
  )

  df_over   <- subset(plot_data, group == "Overall")
  df_male   <- subset(plot_data, group == "Male")
  df_female <- subset(plot_data, group == "Female")

  ################################################################################
  ## Save results
  ################################################################################

  res <- data.frame(
    Disease       = Disease,
    ltr.overall   = fmt_ci(ltr_over,   se_over),
    ltr.Male      = fmt_ci(ltr_male,   se_male),
    ltr.Female    = fmt_ci(ltr_female, se_female),
    ltr.diff      = fmt_ci(ltr_diff,   se_diff),
    pval          = summary(fit_by_sex)$coefficients[2, 4],
    first.overall = first_age_above(df_over,   0.01),
    ten.overall   = first_age_above(df_over,   0.10),
    first.Male    = first_age_above(df_male,   0.01),
    ten.Male      = first_age_above(df_male,   0.10),
    first.Female  = first_age_above(df_female, 0.01),
    ten.Female    = first_age_above(df_female, 0.10)
  )

  if (saving) {
    res_file  <- paste0(paths$out_dir, Disease, "_Lifetime_risk_",
                        start_year, ".csv")
    plot_file <- paste0(paths$out_dir, "plot_data/", Disease,
                        "_Lifetime_risk_plot_", start_year, ".csv")
    write.csv(res,       res_file,  row.names = FALSE)
    write.csv(plot_data, plot_file, row.names = FALSE)
    message("Results saved to ", res_file)
    message("Plot data saved to ", plot_file)
  }

  return(res)
}


################################################################################
## Example of how the function is used
################################################################################


MI_ltr <- calculate_lifetime_risk (Disease    = "MI",
                                   index_age  = 18,
                                   exit_age   = 95,
                                   start_year = 1990)
