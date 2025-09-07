# ---------------------------------------------------------------------------- #
# Clean Managing Anxiety Demographics Data and Create Table
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# The present script cleans demographics data, using information from (a)
# "Script0_Demographics.R" on the OSF project for the Managing Anxiety study main 
# outcomes paper (see https://doi.org/10.17605/OSF.IO/3B67V), (b) from "R34.ipynb"
# (likely written by Sonia Baee based on Claudia Calicho-Mamani's following script)
# and "R34_cleaning_script.R" (likely written by Claudia Calicho-Mamani) in the 
# study's Data Cleaning folder (https://bit.ly/3CLi5It) on GitHub.

# ---------------------------------------------------------------------------- #
# Store working directory, check correct R version, load packages ----
# ---------------------------------------------------------------------------- #

# Store working directory

wd_dir <- getwd()

# Load custom functions

source("./code/01a_define_functions.R")

# Check correct R version, load groundhog package, and specify groundhog_day

groundhog_day <- version_control()

# Load packages

pkgs <- c("flextable", "officer", "ftExtra")
groundhog.library(pkgs, groundhog_day)

# Set "flextable" package defaults and load "officer" package properties

source("./code/01b_set_flextable_defaults.R")
source("./code/01c_set_officer_properties.R")

# ---------------------------------------------------------------------------- #
# Import data ----
# ---------------------------------------------------------------------------- #

load("./data/ma/intermediate/dat_ma_dem.RData")

# Import raw demographics data

load("./data/ma/intermediate/raw_dat_son_a_dem.RData")
load("./data/ma/intermediate/raw_dat_son_b_dem.RData")

# ---------------------------------------------------------------------------- #
# Clean age ----
# ---------------------------------------------------------------------------- #

# "age" and corresponding "birthYear" are missing for 12 participants, but the
# "birthYear"s are present in raw data files. Although "R34.ipynb" (see above)
# presumably cleans "birthYear" (as does "R34_cleaning_script.R", see above), it
# seems that the "birthYears" for all these participants were recoded as NA.

sum(is.na(dat_ma_dem$age)) == 12
sum(is.na(dat_ma_dem$birthYear)) == 12

missing_birthYear_ids <- dat_ma_dem$participant_id[is.na(dat_ma_dem$birthYear)]

raw_dat_son_a_tmp <- raw_dat_son_a_dem[raw_dat_son_a_dem$participantRSA %in% missing_birthYear_ids, ]
raw_dat_son_b_tmp <- raw_dat_son_b_dem[raw_dat_son_b_dem$participantRSA %in% missing_birthYear_ids, ]

raw_dat_son_a_tmp$participantRSA <- as.integer(raw_dat_son_a_tmp$participantRSA)

# View(raw_dat_son_a_tmp[order(raw_dat_son_a_tmp$participantRSA), ])
# View(raw_dat_son_b_tmp[order(raw_dat_son_b_tmp$participantRSA), ])

# With one exception (i.e., "birthYear" of 2222 for participant 1344), the values
# for "birthYear" from "raw_dat_son_b" seem correct. Use those "birthYear"s and the
# "date"s of the responses so that "age"s can be computed below.

names(raw_dat_son_b_tmp)[names(raw_dat_son_b_tmp) == "participantRSA"] <- "participant_id"

dat_ma_dem$birthYear[match(missing_birthYear_ids, dat_ma_dem$participant_id)] <-
  raw_dat_son_b_tmp$birthYear[match(missing_birthYear_ids, raw_dat_son_b_tmp$participant_id)]

dat_ma_dem$date <- NA

dat_ma_dem$date[match(missing_birthYear_ids, dat_ma_dem$participant_id)] <-
  raw_dat_son_b_tmp$date[match(missing_birthYear_ids, raw_dat_son_b_tmp$participant_id)]

# Recode "2222" to NA, leaving only 2 values of NA (one was "2222"; another was "0")

dat_ma_dem$birthYear[dat_ma_dem$birthYear == "2222"] <- NA

sum(is.na(dat_ma_dem$birthYear)) == 2

# Compute "age" for participants with previously missing "birthYear"s as difference 
# between "birthYear" and year of "date" (following "R34.ipynb")

dat_ma_dem$date <- as.POSIXct(dat_ma_dem$date, tz = "EST", format = "%a, %d %b %Y %H:%M:%S %z")

dat_ma_dem$age[dat_ma_dem$participant_id %in% missing_birthYear_ids] <-
  as.numeric(format(dat_ma_dem$date[dat_ma_dem$participant_id %in% missing_birthYear_ids], "%Y")) -
  dat_ma_dem$birthYear[dat_ma_dem$participant_id %in% missing_birthYear_ids]

dat_ma_dem$date <- NULL

# Note: Range is reasonable

range(dat_ma_dem$age, na.rm = TRUE) == c(18, 91)

# ---------------------------------------------------------------------------- #
# Clean gender ----
# ---------------------------------------------------------------------------- #

pna  <- "Prefer not to answer"
pnam <- "Prefer not to answer or missing"

# Recode levels

dat_ma_dem$gender[dat_ma_dem$gender == ""]  <- NA

sum(is.na(dat_ma_dem$gender)) == 1
sum(dat_ma_dem$gender == pna, na.rm = TRUE) == 5

dat_ma_dem$gender[is.na(dat_ma_dem$gender) | dat_ma_dem$gender == pna] <- pnam

# Reorder levels

dat_ma_dem$gender <- 
  factor(dat_ma_dem$gender,
         levels = c("Female", "Male", "Transgender", "Other", pnam))

# ---------------------------------------------------------------------------- #
# Clean race ----
# ---------------------------------------------------------------------------- #

# Recode levels

dat_ma_dem$race[dat_ma_dem$race == ""] <- NA

sum(is.na(dat_ma_dem$race)) == 2
sum(dat_ma_dem$race == pna, na.rm = TRUE) == 30

dat_ma_dem$race[is.na(dat_ma_dem$race) | dat_ma_dem$race == pna] <- pnam

# Reorder levels

dat_ma_dem$race <-
  factor(dat_ma_dem$race,
         levels = c("American Indian/Alaska Native", "Black/African origin",
                    "East Asian", "Native Hawaiian/Pacific Islander", "South Asian",
                    "White/European origin", "Other or Unknown", pnam))

# ---------------------------------------------------------------------------- #
# Clean ethnicity ----
# ---------------------------------------------------------------------------- #

# Recode level

dat_ma_dem$ethnicity[dat_ma_dem$ethnicity == ""] <- NA

sum(is.na(dat_ma_dem$ethnicity)) == 3
sum(dat_ma_dem$ethnicity == pna, na.rm = TRUE) == 40

dat_ma_dem$ethnicity[is.na(dat_ma_dem$ethnicity) | dat_ma_dem$ethnicity == pna] <- pnam

# Reorder levels

dat_ma_dem$ethnicity <- 
  factor(dat_ma_dem$ethnicity,
         levels = c("Hispanic or Latino", "Not Hispanic or Latino", 
                    "Unknown", pnam))

# ---------------------------------------------------------------------------- #
# Clean education ----
# ---------------------------------------------------------------------------- #

# Note: "R34.ipynb" (see above) does not clean education (likely an oversight),
# whereas "R34_cleaning_script.R" (see above) does. Use information from the
# latter to clean education.

# Recode levels

dat_ma_dem$education[dat_ma_dem$education %in% c("", "????")] <- NA
dat_ma_dem$education[dat_ma_dem$education == "Diploma di scuola superiore"] <- "High School Graduate"
dat_ma_dem$education[dat_ma_dem$education == "Un lyc̩e"]                     <- "Some High School"

sum(is.na(dat_ma_dem$education)) == 2
sum(dat_ma_dem$education == pna, na.rm = TRUE) == 9

dat_ma_dem$education[is.na(dat_ma_dem$education) | dat_ma_dem$education == pna] <- pnam

# Reorder levels

dat_ma_dem$education <-
  factor(dat_ma_dem$education,
         levels = c("Elementary School", "Junior High", "Some High School", 
                    "High School Graduate", "Some College", "Associate's Degree", 
                    "Bachelor's Degree", "Some Graduate School", "Master's Degree", 
                    "M.B.A.", "J.D.", "M.D.", "Ph.D.", "Other Advanced Degree",
                    pnam))

# ---------------------------------------------------------------------------- #
# Clean employment status ----
# ---------------------------------------------------------------------------- #

# Recode levels

dat_ma_dem$employmentStatus[dat_ma_dem$employmentStatus == ""] <- NA
dat_ma_dem$employmentStatus[dat_ma_dem$employmentStatus == 
                              "Homemaker/keeping house or raising children full-time"] <- "Homemaker"

sum(is.na(dat_ma_dem$employmentStatus)) == 2
sum(dat_ma_dem$employmentStatus == pna, na.rm = TRUE) == 9

dat_ma_dem$employmentStatus[is.na(dat_ma_dem$employmentStatus) | dat_ma_dem$employmentStatus == pna] <- pnam

# Reorder levels

dat_ma_dem$employmentStatus <- 
  factor(dat_ma_dem$employmentStatus,
         levels = c("Student", "Homemaker", "Unemployed or laid off", "Looking for work",
                    "Working part-time", "Working full-time", "Retired", "Other", pnam))

# ---------------------------------------------------------------------------- #
# Clean income ----
# ---------------------------------------------------------------------------- #

# Recode levels

dat_ma_dem$income[dat_ma_dem$income == ""] <- NA
dat_ma_dem$income[dat_ma_dem$income == "Don't know"] <- "Unknown"

sum(is.na(dat_ma_dem$income)) == 1
sum(dat_ma_dem$income == pna, na.rm = TRUE) == 92

dat_ma_dem$income[is.na(dat_ma_dem$income) | dat_ma_dem$income == pna] <- pnam

# Reorder levels

dat_ma_dem$income <-
  factor(dat_ma_dem$income,
         levels = c("Less than $5,000", "$5,000 through $11,999", 
                    "$12,000 through $15,999", "$16,000 through $24,999", 
                    "$25,000 through $34,999", "$35,000 through $49,999",
                    "$50,000 through $74,999", "$75,000 through $99,999",
                    "$100,000 through $149,999", "$150,000 through $199,999",
                    "$200,000 through $249,999", "$250,000 or greater",
                    "Unknown", pnam))

# ---------------------------------------------------------------------------- #
# Clean marital status ----
# ---------------------------------------------------------------------------- #

# Recode levels

dat_ma_dem$maritalStatus[dat_ma_dem$maritalStatus == ""] <- NA
dat_ma_dem$maritalStatus[dat_ma_dem$maritalStatus == 
                           "Single, but casually dating"] <- 
  "Dating"
dat_ma_dem$maritalStatus[dat_ma_dem$maritalStatus == 
                           "Single, but currently engaged to be married"] <- 
  "Engaged"
dat_ma_dem$maritalStatus[dat_ma_dem$maritalStatus == 
                           "Single, but currently living with someone in a marriage-like relationship"] <- 
  "In marriage-like relationship"
dat_ma_dem$maritalStatus[dat_ma_dem$maritalStatus == 
                           "In a domestic or civil union"] <- 
  "In domestic or civil union"

sum(is.na(dat_ma_dem$maritalStatus)) == 3
sum(dat_ma_dem$maritalStatus == pna, na.rm = TRUE) == 13

dat_ma_dem$maritalStatus[is.na(dat_ma_dem$maritalStatus) | dat_ma_dem$maritalStatus == pna] <- pnam

# Reorder levels

dat_ma_dem$maritalStatus <-
  factor(dat_ma_dem$maritalStatus,
         levels = c("Single", "Dating", "Engaged", "In marriage-like relationship",
                    "Married", "In domestic or civil union", "Separated", "Divorced", 
                    "Widow/widower", "Other", pnam))

# ---------------------------------------------------------------------------- #
# Clean country ----
# ---------------------------------------------------------------------------- #

# Recode levels

dat_ma_dem$residenceCountry[dat_ma_dem$residenceCountry == ""] <- NA
dat_ma_dem$residenceCountry[dat_ma_dem$residenceCountry == "NoAnswer"] <- pna

sum(is.na(dat_ma_dem$residenceCountry)) == 1
sum(dat_ma_dem$residenceCountry == pna, na.rm = TRUE) == 2

dat_ma_dem$residenceCountry[is.na(dat_ma_dem$residenceCountry) | dat_ma_dem$residenceCountry == pna] <- pnam

# Define desired levels order (decreasing frequency ending with "Prefer not to answer or missing")

country_levels <- names(sort(table(dat_ma_dem$residenceCountry), decreasing = TRUE))
country_levels <- c(country_levels[country_levels != pnam], pnam)

# Reorder levels of "country"

dat_ma_dem$residenceCountry <- factor(dat_ma_dem$residenceCountry, levels = country_levels)

# Define "country_col", collapsing countries with fewer than 10 participants into "Other"

top_countries <- 
  names(table(dat_ma_dem$residenceCountry)[as.numeric(table(dat_ma_dem$residenceCountry)) > 10])

all(top_countries == c("United States", "Canada", "United Kingdom"))

dat_ma_dem$country_col <- NA

for (i in 1:nrow(dat_ma_dem)) {
  if (is.na(dat_ma_dem$residenceCountry)[i]) {
    dat_ma_dem$country_col[i] <- NA
  } else if (as.character(dat_ma_dem$residenceCountry)[i] %in% top_countries) {
    dat_ma_dem$country_col[i] <- as.character(dat_ma_dem$residenceCountry)[i]
  } else if (as.character(dat_ma_dem$residenceCountry)[i] == pnam) {
    dat_ma_dem$country_col[i] <- pnam
  } else {
    dat_ma_dem$country_col[i] <- "Other"
  }
}

# Reorder levels of "country_col"

dat_ma_dem$country_col <-
  factor(dat_ma_dem$country_col,
         levels = c(top_countries, "Other", pnam))

# ---------------------------------------------------------------------------- #
# Save cleaned data ----
# ---------------------------------------------------------------------------- #

dat_ma_dem2 <- dat_ma_dem

save(dat_ma_dem2, file = "./data/ma/intermediate/dat_ma_dem2.RData")

# ---------------------------------------------------------------------------- #
# Create demographics table ----
# ---------------------------------------------------------------------------- #

dem_tbl_ma <- dat_ma_dem2

# Define function to compute descriptives for Managing Anxiety

compute_desc_ma <- function(df) {
  # Compute sample size
  
  n <- data.frame(label = "n",
                  value = length(df$participant_id))
  
  # Compute mean and standard deviation for numeric variables
  
  num_res <- rbind(data.frame(label = "Age",
                              value = NA),
                   data.frame(label = "Years: M (SD)",
                              value = paste0(format(round(mean(df$age, na.rm = TRUE), 2),
                                                    nsmall = 2, trim = TRUE), 
                                             " (",
                                             format(round(sd(df$age, na.rm = TRUE), 2),
                                                    nsmall = 2, trim = TRUE),
                                             ")")))
  
  # Compute count and percentage "prefer not to answer or missing" for numeric variables
  
  num_res_pnam <- data.frame(label = "Prefer not to answer or missing: n (%)",
                             value = paste0(sum(is.na(df$age)),
                                            " (",
                                            format(round((sum(is.na(df$age))/length(df$age))*100, 1),
                                                   nsmall = 1, trim = TRUE),
                                            ")"))
  
  # Compute count and percentage for factor variables
  
  vars <- c("gender", "race", "ethnicity", "country_col", "education",
            "employmentStatus", "income", "maritalStatus")
  var_labels <- paste0(c("Gender", "Race", "Ethnicity", "Country", "Education",
                         "Employment Status", "Annual Income", "Marital Status"),
                       ": n (%)")
  
  fct_res <- data.frame()
  
  for (i in 1:length(vars)) {
    tbl <- table(df[, vars[i]])
    prop_tbl <- prop.table(tbl)*100
    
    tbl_res <- rbind(data.frame(label = var_labels[i],
                                value = NA),
                     data.frame(label = names(tbl),
                                value = paste0(as.numeric(tbl),
                                               " (", 
                                               format(round(as.numeric(prop_tbl), 1),
                                                      nsmall = 1, trim = TRUE),
                                               ")")))
    fct_res <- rbind(fct_res, tbl_res)
  }
  
  # Combine results
  
  res <- rbind(n, num_res, num_res_pnam, fct_res)
  
  return(res)
}

# Compute descriptives across conditions

res_across_cond_ma <- compute_desc_ma(dem_tbl_ma)

# Save table to CSV

dem_path <- "./results/demographics/"

dir.create(dem_path, recursive = TRUE)

write.csv(res_across_cond_ma,
          paste0(dem_path, "across_cond_ma.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Format demographics table ----
# ---------------------------------------------------------------------------- #

# "flextable" defaults are set in "set_flextable_defaults.R" above

# Section and text properties are sourced from "set_officer_properties.R" above

# Define table notes

gen_note <- as_paragraph_md("*Note.* Each characteristic has missing data for 1-3 participants due to a server error.")

footnotes <- list(country_other = "\\ Countries with fewer than 10 participants were collapsed into Other.")

# Run function

dem_tbl_ma_ft <-
  format_dem_tbl(res_across_cond_ma, footnotes, "Demographic Characteristics in Managing Anxiety Study", 29)

# Export flextable

save(dem_tbl_ma_ft, file = "./results/demographics/dem_tbl_ma_ft.RData")