# CV Lifecourse

Analytical code for the lifetime risk and life years lost analyses in:

**The cardiovascular disease landscape across the life course: a nationwide population study**

The code estimates, for a given diagnosis in the Swedish adult population:

1. **Lifetime risk**: the cumulative risk of the diagnosis from index age (usually 18) to age 95, accounting for left truncation and the competing risk of death.
2. **Excess life years lost (LYL)**: the life years lost before age 95 after the diagnosis, compared with the general Swedish population of the same sex, age and calendar period.

Both analyses report results overall and by sex.

> **Data access:** This repository contains code only. The analyses use pseudonymised individual-level data from Swedish national registers, held in a secure trusted research environment (TRE). The data cannot be shared; access requires ethical approval and permission from the data holders.

---

## Repository structure

```
.
├── R code/
│   ├── lifetime_risk_analysis.R   # Lifetime risk (Aalen-Johansen + pseudo-observations)
│   └── LYL_analysis.R             # Excess life years lost (lillies)
└── README.md
```

Each script is self-contained: it loads its packages, defines file paths and helper functions, defines one main function, and ends with an example call for myocardial infarction (`"MI"`).

---

## Requirements

R with the following packages:

| Package     | Lifetime risk | LYL | Used for |
|-------------|:-------------:|:---:|----------|
| `arrow`     | ✓ | ✓ | Reading parquet register files |
| `dplyr`     | ✓ | ✓ | Data manipulation |
| `lubridate` | ✓ | ✓ | Date handling |
| `zoo`       | ✓ | ✓ | Year-month death dates (`as.yearmon`) |
| `survival`  | ✓ |   | Aalen-Johansen estimator (`survfit`) and pseudo-observations (`pseudo`) |
| `survey`    | ✓ |   | GLM with robust standard errors (`svyglm`) |
| `msm`       | ✓ |   | Delta-method standard errors (`deltamethod`) |
| `lillies`   |   | ✓ | Life years lost and bootstrap CIs (`lyl_range`, `lyl_ci`) |

```r
install.packages(c("arrow", "dplyr", "lubridate", "zoo",
                   "survival", "survey", "msm", "lillies"))
```

`survival::pseudo()` is only available in recent versions of `survival`. 


## Configuration

File paths and the end of follow-up are set at the top of each script and point to the TRE. Update them for your environment:

```r
paths <- list(
  disease_dir = "...",   # folder with diagnosis files
  tpr         = "...",   # Total Population Register (lifetime risk only)
  lifetables  = "...",   # life table CSV (LYL only)
  cdr         = "...",   # Cause of Death Register, main file
  cdr_avi     = "...",   # Cause of Death Register, supplementary file
  migration   = "...",   # migration register
  out_dir     = "..."    # output folder
)
```

| Script | Setting | Value |
|--------|---------|-------|
| `lifetime_risk_analysis.R` | `END_OF_FOLLOWUP_YEAR` | `2024` |
| `LYL_analysis.R` | `END_OF_FOLLOWUP_DATE` | `as.Date("2024-01-01")` |

---

## 1. Lifetime risk

### Usage

```r
source("R code/lifetime_risk_analysis.R")

MI_ltr <- calculate_lifetime_risk(
  Disease    = "MI",   # must match the diagnosis file name
  index_age  = 18,     # youngest age at study entry
  exit_age   = 95,     # age at which lifetime risk is evaluated
  start_year = 1990,   # first calendar year of follow-up (function default: 2014)
  saving     = TRUE    # write results and CIF curves to CSV
)
```


### Method

1. Cumulative incidence functions (CIFs) are estimated with the **Aalen-Johansen estimator** (`survival::survfit` with a multi-state outcome), overall and by sex.
2. **Pseudo-observations** of the CIF at age `exit_age` are computed for every individual (`survival::pseudo`).
3. The pseudo-observations are regressed on sex in a **generalised linear model with a complementary log-log link** and robust (sandwich) standard errors (`survey::svyglm`, one cluster per person). A bounded cloglog link (`bcloglog()`) is used because pseudo-observations can fall outside (0, 1); it only affects starting values, not the fitted estimates.
4. Lifetime risks, the female–male difference and 95% CIs are derived from the model with the **delta method** (`msm::deltamethod`).
---

## 2. Excess life years lost

### Usage

```r
source("R code/LYL_analysis.R")

MI_lyl <- calculate_LYL(
  Disease    = "MI",                          # must match the diagnosis file name
  index_age  = 18,                            # minimum age at diagnosis
  tau        = 95,                            # upper age limit for LYL
  start_year = 1990,                          # first calendar year of diagnosis
  groups     = c("Overall", "Male", "Female"),
  n_iter     = 1000,                          # bootstrap iterations
  saving     = TRUE                           # write results to CSV
)
```

### Method

1. **LYL among cases** between age at diagnosis and age `tau` is estimated with the method of Andersen (2013, 2017) using `lillies::lyl_range`, with bootstrap confidence intervals (`lillies::lyl_ci`). Estimates are averaged over the cases' distribution of age at diagnosis.
2. **LYL in the general population** is derived from Swedish life tables (remaining life expectancy by sex, age and calendar year), weighted by the cases' distributions of calendar year and age at diagnosis.
3. **Excess LYL** = LYL among cases − LYL in the general population.


---

---

## Citation

If you use this code, please cite the paper above. <!-- TODO: full citation -->

## Contact

Konstantinos Konstantinou, konstantinos.konstantinou@apnc.se

