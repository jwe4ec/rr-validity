# ---------------------------------------------------------------------------- #
# Clean Calm Thinking Demographics Data and Create Table
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Although the present script creates a demographics table for Calm Thinking, we 
# decided not to analyze Calm Thinking data as part of the dissertation's scope

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

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

load("./data/ct/intermediate_clean_further/dat_ct_rr.RData")

# ---------------------------------------------------------------------------- #
# Recode "Prefer not to answer" in categorical and ordinal variables ----
# ---------------------------------------------------------------------------- #

pna <- "Prefer not to answer"

# Recode 555 as "Prefer not to answer". 555 in "birth_year" (integer variable)
# is recoded after "age" is computed below.

target_vars <- c("education", "employment_stat", "ethnicity", "gender", "income",
                 "marital_stat")
dat_ct_rr$demographics[, target_vars][dat_ct_rr$demographics[, target_vars] == 555] <- pna

dat_ct_rr$demographics_race$race[dat_ct_rr$demographics_race$race == 555] <- pna

# Recode "NoAnswer" as "Prefer not to answer"

dat_ct_rr$demographics$country[dat_ct_rr$demographics$country == "NoAnswer"] <- pna

# ---------------------------------------------------------------------------- #
# Compute age ----
# ---------------------------------------------------------------------------- #

# Compute age (use NA where "birth_year" is 555 ["prefer not to answer"])

dat_ct_rr$demographics$age <- NA

for (i in 1:nrow(dat_ct_rr$demographics)) {
  if (dat_ct_rr$demographics$birth_year[i] != 555) {
    dat_ct_rr$demographics$age[i] <- 
      as.integer(format(dat_ct_rr$demographics$date_as_POSIXct[i], "%Y")) - 
      dat_ct_rr$demographics$birth_year[i]
  }
}

# Participants range from 18 to 75 years of age, which is reasonable

range(dat_ct_rr$demographics$age, na.rm = TRUE) == c(18, 75)

# All NAs in "age" are due to "birth_year" of 555

all(dat_ct_rr$demographics$birth_year[is.na(dat_ct_rr$demographics$age)] == 555)

# Recode 555 in "birth_year" to "Prefer not to answer"

dat_ct_rr$demographics$birth_year[dat_ct_rr$demographics$birth_year == 555] <- pna

# ---------------------------------------------------------------------------- #
# Clean gender ----
# ---------------------------------------------------------------------------- #

# Reorder levels and add "Transgender Female". Codebook indicates that on 8/5/2019, 
# "Transgender" option was replaced with "Transgender Male" and "Transgender Female".

dat_ct_rr$demographics$gender <- 
  factor(dat_ct_rr$demographics$gender,
         levels = c("Female", "Male", 
                    "Transgender", "Transgender Female", "Transgender Male",
                    "Other", "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean ethnicity ----
# ---------------------------------------------------------------------------- #

# Reorder levels

dat_ct_rr$demographics$ethnicity <- 
  factor(dat_ct_rr$demographics$ethnicity,
         levels = c("Hispanic or Latino", "Not Hispanic or Latino", 
                    "Unknown", "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean education ----
# ---------------------------------------------------------------------------- #

# Reorder levels

dat_ct_rr$demographics$education <-
  factor(dat_ct_rr$demographics$education,
         levels = c("Elementary School", "Junior High", "Some High School", 
                    "High School Graduate", "Some College", "Associate's Degree", 
                    "Bachelor's Degree", "Some Graduate School", "Master's Degree", 
                    "M.B.A.", "J.D.", "M.D.", "Ph.D.", "Other Advanced Degree",
                    "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean employment status ----
# ---------------------------------------------------------------------------- #

# Recode "Homemaker"

dat_ct_rr$demographics$employment_stat[dat_ct_rr$demographics$employment_stat ==
                                         "Homemaker/keeping house or raising children full-time"] <-
  "Homemaker"

# Reorder levels

dat_ct_rr$demographics$employment_stat <- 
  factor(dat_ct_rr$demographics$employment_stat,
         levels = c("Student", "Homemaker", "Unemployed or laid off", "Looking for work",
                    "Working part-time", "Working full-time", "Retired", "Other",
                    "Unknown", "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean income ----
# ---------------------------------------------------------------------------- #

# Recode "Don't know"

dat_ct_rr$demographics$income[dat_ct_rr$demographics$income == "Don't know"] <- "Unknown"

# Reorder levels

dat_ct_rr$demographics$income <-
  factor(dat_ct_rr$demographics$income,
         levels = c("Less than $5,000", "$5,000 through $11,999", 
                    "$12,000 through $15,999", "$16,000 through $24,999", 
                    "$25,000 through $34,999", "$35,000 through $49,999",
                    "$50,000 through $74,999", "$75,000 through $99,999",
                    "$100,000 through $149,999", "$150,000 through $199,999",
                    "$200,000 through $249,999", "$250,000 or greater",
                    "Unknown", "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean marital status ----
# ---------------------------------------------------------------------------- #

# Recode levels

dat_ct_rr$demographics$marital_stat[dat_ct_rr$demographics$marital_stat == "Single,dating"]       <- "Dating"
dat_ct_rr$demographics$marital_stat[dat_ct_rr$demographics$marital_stat == "Single,engaged"]      <- "Engaged"
dat_ct_rr$demographics$marital_stat[dat_ct_rr$demographics$marital_stat == "Single,marriagelike"] <- "In marriage-like relationship"
dat_ct_rr$demographics$marital_stat[dat_ct_rr$demographics$marital_stat == "civilunion"]          <- "In domestic or civil union"

# Reorder levels

dat_ct_rr$demographics$marital_stat <-
  factor(dat_ct_rr$demographics$marital_stat,
         levels = c("Single", "Dating", "Engaged", "In marriage-like relationship",
                    "Married", "In domestic or civil union", "Separated", "Divorced", 
                    "Widow/widower", "Other", "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean race ----
# ---------------------------------------------------------------------------- #

# No duplicated responses exist

nrow(dat_ct_rr$demographics_race[duplicated(dat_ct_rr$demographics_race[, c("participant_id", 
                                                                            "race")]), ]) == 0

# Compute number of responses per "participant_id"

race_ag <- aggregate(race ~ participant_id,
                     dat_ct_rr$demographics_race,
                     length)
names(race_ag)[names(race_ag) == "race"] <- "race_cnt"

dat_ct_rr$demographics_race <- merge(dat_ct_rr$demographics_race, race_ag,
                                     by = "participant_id", all.x = TRUE)

# Define "race_col", collapsing cases of more than once race

dat_ct_rr$demographics_race$race_col <- NA

for (i in 1:nrow(dat_ct_rr$demographics_race)) {
  if (dat_ct_rr$demographics_race$race_cnt[i] == 1) {
    dat_ct_rr$demographics_race$race_col[i] <- dat_ct_rr$demographics_race$race[i]
  } else if (dat_ct_rr$demographics_race$race_cnt[i] > 1) {
    dat_ct_rr$demographics_race$race_col[i] <- "More than one race"
  }
}

# TODO: Exclude participant 1992 from "demographics_race" table in centralized Calm 
# Thinking data cleaning. Then it will not need to be done here.

dat_ct_rr$demographics_race <- dat_ct_rr$demographics_race[dat_ct_rr$demographics_race$participant_id != 1992, ]





# Add "race_col" to "demographics" table

dat_ct_rr$demographics <- merge(dat_ct_rr$demographics,
                          unique(dat_ct_rr$demographics_race[, c("participant_id", 
                                                                 "race_col")]),
                          by = "participant_id",
                          all.x = TRUE)

# Reorder levels

dat_ct_rr$demographics$race_col <-
  factor(dat_ct_rr$demographics$race_col,
         levels = c("American Indian/Alaska Native", "Black/African origin",
                    "East Asian", "Native Hawaiian/Pacific Islander", "South Asian",
                    "White/European origin", "More than one race", "Other or Unknown", 
                    "Prefer not to answer"))

# ---------------------------------------------------------------------------- #
# Clean country ----
# ---------------------------------------------------------------------------- #

# Recode "Ã…land Islands"

dat_ct_rr$demographics$country[dat_ct_rr$demographics$country == "Ã…land Islands"] <- "Åland Islands"

# Define desired levels order (decreasing frequency ending with "prefer not to answer")

country_levels <- names(sort(table(dat_ct_rr$demographics$country), decreasing = TRUE))
country_levels <- c(country_levels[country_levels != pna], pna)

# Reorder levels of "country"

dat_ct_rr$demographics$country <- factor(dat_ct_rr$demographics$country, levels = country_levels)

# Define "country_col", collapsing countries with fewer than 10 participants into "Other"

top_countries <- 
  names(table(dat_ct_rr$demographics$country)[as.numeric(table(dat_ct_rr$demographics$country)) > 10])

all(top_countries == c("United States", "Australia", "Canada", "United Kingdom"))

dat_ct_rr$demographics$country_col <- NA

for (i in 1:nrow(dat_ct_rr$demographics)) {
  if (as.character(dat_ct_rr$demographics$country)[i] %in% top_countries) {
    dat_ct_rr$demographics$country_col[i] <- as.character(dat_ct_rr$demographics$country)[i]
  } else if (as.character(dat_ct_rr$demographics$country)[i] == pna) {
    dat_ct_rr$demographics$country_col[i] <- pna
  } else {
    dat_ct_rr$demographics$country_col[i] <- "Other"
  }
}

# Reorder levels of "country_col"

dat_ct_rr$demographics$country_col <-
  factor(dat_ct_rr$demographics$country_col,
         levels = c(top_countries, "Other", pna))

# ---------------------------------------------------------------------------- #
# Create demographics table ----
# ---------------------------------------------------------------------------- #

dem_tbl_ct <- dat_ct_rr$demographics

# Define function to compute descriptives for Calm Thinking

compute_desc_ct <- function(df) {
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
  
  # Compute count and percentage "prefer not to answer" for numeric variables
  
  num_res_pna <- data.frame(label = "Prefer not to answer: n (%)",
                            value = paste0(sum(is.na(df$age)),
                                           " (",
                                           format(round((sum(is.na(df$age))/length(df$age))*100, 1),
                                                  nsmall = 1, trim = TRUE),
                                           ")"))
  
  # Compute count and percentage for factor variables
  
  vars <- c("gender", "race_col", "ethnicity", "country_col", "education",
            "employment_stat", "income", "marital_stat")
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
  
  res <- rbind(n, num_res, num_res_pna, fct_res)
  
  return(res)
}

# Compute descriptives across conditions

res_across_cond_ct <- compute_desc_ct(dem_tbl_ct)

# Save table to CSV

write.csv(res_across_cond_ct, "./results/demographics/across_cond_ct.csv", row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Format demographics table ----
# ---------------------------------------------------------------------------- #

# "flextable" defaults are set in "set_flextable_defaults.R" above

# Section and text properties are sourced from "set_officer_properties.R" above

# Define table notes

gen_note <- NULL

footnotes <- list(country_other = "\\ Countries with fewer than 10 participants were collapsed into Other.")

# Run function

dem_tbl_ct_ft <-
  format_dem_tbl(res_across_cond_ct, footnotes, "Demographic Characteristics in Calm Thinking Study", 33)

# Export flextable

save(dem_tbl_ct_ft, file = "./results/demographics/dem_tbl_ct_ft.RData")