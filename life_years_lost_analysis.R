###############################################################################
#  Code for estimating excess life years lost after a diagnosis in the Swedish population
#
#  Purpose:
#
#    Estimate the life years lost(LYL) before age 95 among people diagnosed
#    with a given disease, and the excess LYL compared with the general
#    Swedish population of the same age, sex and calendar period. Results
#    are produced overall and by sex.
#
#  Methods:
#
#    1. Cases are people with a first recorded diagnosis from start_year
#       onwards, aged index_age or older at diagnosis. Follow-up runs from
#       diagnosis to death, emigration or end of follow-up.
#
#    2. LYL among cases between age at diagnosis and age 95 is estimated
#       with the method of Andersen (2013, 2017) as implemented in the
#       lillies package (lillies::lyl_range), with bootstrap confidence intervals
#       (lillies::lyl_ci). Estimates are averaged over the age at diagnosis
#       distribution of the cases.
#
#    3. LYL in the general population is derived from official Swedish
#       life tables (remaining life expectancy by sex, age and calendar
#       year), weighted by the calendar year and age distributions of the
#       cases' diagnoses.
#
#    4. Excess LYL = LYL among cases - LYL in the general population.
#
###############################################################################


################################################################################
## Load Packages
################################################################################

library(arrow)     
library(dplyr)     
library(lubridate)  
library(zoo)     
library(lillies)    

set.seed(123)     


################################################################################
## Configuration options
################################################################################
# paths in the TRE environment
paths <- list(
  disease_dir = "/safe/data/APNC/DATA/ANALYSIS DATA/",
  lifetables  = "/safe/data/APNC/SCARP-Lifecourse/Data/Derived/Lifetables_per_year.csv",
  cdr         = "/safe/data/SCARP/DATA_2024/PARQUET/ut_r_dors_19328_2023.parquet",
  cdr_avi     = "/safe/data/SCARP/DATA_2024/PARQUET/ut_r_dors_avi_19328_2023.parquet",
  migration   = "/safe/data/SCARP/DATA_2024/PARQUET/ar_lev_migrationer.parquet",
  out_dir     = "/safe/data/APNC/SCARP-Lifecourse/Data/Derived/ltr/45/Diagnoses/plot_data/"
)

# Last calendar year of follow-up
END_OF_FOLLOWUP_DATE <- as.Date("2024-01-01") 

# Labels used for each group in the life table file
LIFETABLE_GENDER <- c(Overall = "Overall", Male = "Man", Female = "Woman")


################################################################################
##  Helper functions
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

#' Life years lost before tau in the general population, by sex and age
#'
#' Life expectancy for each sex and age is averaged across calendar years,
#' weighted by the number of cases diagnosed in each year. Years are matched
#' by name, so diagnoses in years without a life table are ignored.
#'
#' @param lifetables Data frame with Gender, age and one column per year
#' @param diag_years Calendar years of diagnosis of the cases.
#' @param tau        Upper age limit (default 95).
#' 
#' @return Data frame with Gender, age, weighted_Life_exp and LYL.

reference_lyl <- function(lifetables, diag_years, tau = 95) {
  year_cols   <- grep("^X[0-9]{4}$", names(lifetables), value = TRUE)
  table_years <- as.integer(sub("^X", "", year_cols))
  year_w      <- as.numeric(table(factor(diag_years, levels = table_years)))

  e_matrix <- as.matrix(lifetables[, year_cols])
  e_weighted <- apply(e_matrix, 1, weighted.mean, w = year_w, na.rm = TRUE)

  lifetables %>%
    select(Gender, age) %>%
    mutate(weighted_Life_exp = e_weighted,
           LYL = tau - age - weighted_Life_exp)
}

################################################################################
#' General-population LYL averaged over the cases' ages at diagnosis
#'
#' @param ref    Output of reference_lyl().
#' @param gender Life table label ("Overall", "Man" or "Woman").
#' @param ages   Ages at diagnosis of the cases in this group.
#' 
#' @return Weighted mean LYL in the general population.

general_pop_lyl <- function(ref, gender, ages) {
  ref_g <- ref[ref$Gender == gender, ]
  n_at_age <- as.numeric(table(factor(ages, levels = ref_g$age)))
  weighted.mean(ref_g$LYL, w = n_at_age, na.rm = TRUE)
}

################################################################################
#' Excess life years lost for one group of cases
#'
#' @param data      Analysis data for the group (one row per person).
#' @param group     Group label: "Overall", "Male" or "Female".
#' @param ref       Output of reference_lyl().
#' @param index_age Age at diagnosis.
#' @param tau       Upper age limit.
#' @param n_iter    Number of bootstrap iterations.
#' 
#' @return One-row data frame with LYL among cases, general-population LYL,
#'   and excess LYL with 95% CI.

estimate_excess_lyl <- function(data, group, ref, index_age, tau, n_iter) {
  message("Calculating life years lost: ", group)

  # LYL before tau among cases, for each age at diagnosis in the range
  fit <- lyl_range(data = data,
                   t0 = alder, 
                   t = Death_age,
                   status = Died,
                   age_begin = index_age, 
                   age_end = tau - 1,
                   tau = tau)

  # Bootstrap CIs.
  fit_ci <- lyl_ci(fit, niter = n_iter)
  fit_ci$LYL_ci[is.na(fit_ci$LYL_ci)] <- 0

  # Averaged over the observed age-at-diagnosis distribution of the cases
  cases_lyl <- as.data.frame(summary(fit_ci, weights = data$alder)) %>%
    select(lyl_estimate.TotalLYL, 
           lyl_ci_left.TotalLYL,
           lyl_ci_right.TotalLYL)

  # General population with the same age distribution
  pop_lyl <- general_pop_lyl(ref, LIFETABLE_GENDER[[group]], data$alder)

  data.frame(
    Group            = group,
    N                = n_distinct(data$lopnr),
    LYL_cases        = cases_lyl$lyl_estimate.TotalLYL,
    LYL_general_pop  = pop_lyl,
    Excess_LYL       = cases_lyl$lyl_estimate.TotalLYL - pop_lyl,
    Excess_LYL_lower = cases_lyl$lyl_ci_left.TotalLYL  - pop_lyl,
    Excess_LYL_upper = cases_lyl$lyl_ci_right.TotalLYL - pop_lyl
  )
}


################################################################################
## Main analysis function
################################################################################

#' Estimate excess life years lost after a diagnosis, overall and by sex
#'
#' @param Disease    Name of the diagnosis; must match the input file name.
#' @param index_age  Age at diagnosis.
#' @param tau        Upper age limit for LYL.
#' @param start_year First calendar year of diagnosis.
#' @param groups     Any of "Overall", "Male", "Female".
#' @param n_iter     Bootstrap iterations for the CI.
#' @param saving     If TRUE, write the results to CSV.
#' 
#' @return Data frame with one row per group.
#' 

calculate_LYL <- function(Disease    = "MI",
                          index_age  = 18,
                          tau        = 95,
                          start_year = 1990,
                          groups     = c("Overall", "Male", "Female"),
                          n_iter     = 1000,
                          saving     = TRUE) {

  message("Analysing ", Disease)

  ##############################################################################
  ##  Read the data
  ##############################################################################
  
  ## Diagnosis data: 
  file_path <- paste0(paths$disease_dir, "ICD_",
                      paste(Disease, "all_unique", "2", "1987", "2024",
                            sep = "_"),
                      ".parquet")
  if (!file.exists(file_path)) stop("Disease file not found: ", file_path)

  disease_data <- read_parquet(file_path) %>%
    filter(AR >= start_year) %>%
    distinct()

  ## Life tables and general population LYL, weighted by year of diagnosis
  
  lifetables <- read.csv(paths$lifetables)[-1]      
  ref <- reference_lyl(lifetables, disease_data$AR, tau)
  message("General-population life years lost calculated")

  ## Cause of Death Register: 
  
  message("Loading CDR data")
  cdr_main <- read_cdr(paths$cdr)
  cdr_avi  <- read_cdr(paths$cdr_avi) %>% anti_join(cdr_main, by = "lopnr")
  cdr      <- bind_rows(cdr_main, cdr_avi)
  rm(cdr_main, cdr_avi); gc()

  ## Emigrations:
  
  message("Loading migration data")
  emig <- open_dataset(paths$migration) %>%
    collect() %>%
    filter(posttyp == "Utv") %>%
    transmute(lopnr, emig_date = as.Date(datum, "%Y%m%d"))

  ##############################################################################
  ##  Cases
  ##############################################################################
  
  cases <- disease_data %>%
    left_join(cdr, by = "lopnr") %>%
    filter(alder >= index_age) %>%
    mutate(Death_dt = as.Date(Death_dt)) %>%   # first day of month of death
    select(lopnr, AR, alder, INDATUM, Death_dt, B_Year, Kon) %>%
    distinct()

  # People with conflicting records (e.g. two different death dates)
  # cannot be resolved and are excluded.
  duplicated_lopnrs <- cases$lopnr[duplicated(cases$lopnr)]
  message("Excluded with conflicting records: ",
          length(unique(duplicated_lopnrs)))
  cases <- cases %>% 
    filter(!lopnr %in% duplicated_lopnrs)

  # First emigration on or after the date of diagnosis
  first_emig <- cases %>%
    select(lopnr, INDATUM) %>%
    inner_join(emig, by = "lopnr", relationship = "many-to-many") %>%
    filter(emig_date >= INDATUM) %>%
    group_by(lopnr) %>%
    summarise(emig_date = min(emig_date), .groups = "drop")

  ##############################################################################
  ## Analysis data set
  ##############################################################################
  
  ## Died      : 1 if a death is recorded, 0 otherwise
  ## Death_age : age at death, or age at censoring (emigration or end of
  ##             follow-up) for those without a recorded death

  analysis_set <- cases %>%
    left_join(first_emig, by = "lopnr") %>%
    mutate(
      censor_date = coalesce(emig_date, END_OF_FOLLOWUP_DATE),

      # Death dates are known only to the month and are set to the 1st.
      # A death in the same month as diagnosis is moved to the diagnosis
      # date so that early deaths are not excluded below.
      Death_dt = if_else(!is.na(Death_dt) & Death_dt < INDATUM &
                           format(Death_dt, "%Y%m") == format(INDATUM, "%Y%m"),
                         INDATUM, Death_dt),

      exit_date  = pmin(Death_dt, censor_date, na.rm = TRUE),
      Time_Death = as.numeric(exit_date - INDATUM + 1),

      Died      = if_else(is.na(Death_dt), 0, 1),
      Death_age = if_else(Died == 1,
                          year(Death_dt)    - B_Year,
                          year(censor_date) - B_Year),

      # Ages are in whole years. lillies requires exit > entry, so people
      # who exit at their age at diagnosis are given one year of follow-up.
      Death_age = if_else(Death_age == alder, Death_age + 1, Death_age)
    ) %>%
    filter(Time_Death > 0) %>%     # drop deaths recorded before diagnosis
    as.data.frame()

  message("Analysis data set created: ", nrow(analysis_set), " individuals")

  ##############################################################################
  ##  Excess life years lost by group
  ##############################################################################
  
  group_data <- list(
    Overall = analysis_set,
    Male    = subset(analysis_set, Kon == 1),
    Female  = subset(analysis_set, Kon == 2)
  )

  res <- bind_rows(lapply(groups, function(g) {
    estimate_excess_lyl(group_data[[g]], g, ref, index_age, tau, n_iter)
  }))
  res <- cbind(Disease = Disease, res)

  ##############################################################################
  ##  Save
  ##############################################################################
  if (saving) {
    out_file <- paste0(paths$out_dir, Disease, "_Life_years_lost.csv")
    write.csv(res, out_file, row.names = FALSE)
    message("Results saved to ", out_file)
  }

  return(res)
}


################################################################################
##  Example
################################################################################
MI_lyl <- calculate_LYL(Disease = "MI", 
                        index_age = 18,
                        tau = 95,
                        start_year = 1990,
                        n_iter = 1000)

