# ---------------------------------------------------------------------------- #
# Prepare Calm Thinking Data
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Although the present script prepares Calm Thinking data, we decided not to 
# analyze Calm Thinking data as part of the dissertation's scope

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# The present script imports selected intermediate clean data files (v1.0.1) from 
# the Public Component of the OSF project for the MindTrails Calm Thinking study 
# (https://osf.io/s8v3h/). The data have been centrally cleaned; for the cleaning 
# scripts, see https://doi.org/10.5281/zenodo.6192907.

# ---------------------------------------------------------------------------- #
# Store working directory, check correct R version, load packages ----
# ---------------------------------------------------------------------------- #

# Store working directory

wd_dir <- getwd()

# Load custom functions

source("./code/01a_define_functions.R")

# Check correct R version, load groundhog package, and specify groundhog_day

groundhog_day <- version_control()

# No packages loaded

# ---------------------------------------------------------------------------- #
# Import intermediate clean data ----
# ---------------------------------------------------------------------------- #

# Obtain file names of selected intermediate clean CSV data files and output a
# warning if they do not contain all those relevant to present paper

int_cln_data_dir <- paste0(wd_dir, "/data/ct/source/intermediate_clean")
filenames <- list.files(int_cln_data_dir, pattern = "\\.csv$", full.names = FALSE)

check_relevant_files(filenames)

# Import tables into list and name tables

dat <- lapply(paste0(int_cln_data_dir, "/", filenames), read.csv)
names(dat) <- sub(".csv", "", filenames)

# Convert system-generated timestamps to POSIXct data types

dat <- convert_POSIXct(dat)

# ---------------------------------------------------------------------------- #
# Import item names ----
# ---------------------------------------------------------------------------- #

load("./data/helper/dat_items.RData")
load("./data/helper/rr_item_map.RData")

# ---------------------------------------------------------------------------- #
# Note on filtering data ----
# ---------------------------------------------------------------------------- #

# Note: As noted in README for centralized data cleaning, end of data collection
# for Calm Thinking study was defined as 12/3/2020, but the last system-generated
# timestamp for a Calm Thinking participant was "2020-11-13 22:13:27 EST". We will
# follow the main outcomes paper and analyze data collected through 11/27/2020;
# thus, no filtering of data is needed.

# ---------------------------------------------------------------------------- #
# Compute participant flow ----
# ---------------------------------------------------------------------------- #

# Note: Numbers screened (3519), ineligible (for various reasons, 774 + 111 + 23), 
# eligible but not enrolled (2611), and enrolled (1748) were computed as part of
# centralized data cleaning. For ineligible participants with multiple screening
# attempts, their reason for ineligibility is based on their most recent attempt.

# Confirm number enrolled

enrolled_ids <- dat$participant$participant_id[!is.na(dat$participant$participant_id)]
length(enrolled_ids) == 1748

# Compute number Stage 1 randomized

stg1_ids <- dat$study$participant_id[dat$study$conditioning != "NONE"]
length(stg1_ids) == 1614

# Compute number enrolled but not Stage 1 randomized. These participants did not
# complete "demographic" table, whose data are required for Stage 1 randomization,
# which was stratified by gender and baseline anxiety symptom severity

length(enrolled_ids) - length(stg1_ids) == 134

nrow(dat$demographics[!(dat$demographics$participant_id %in% stg1_ids), ]) == 0

# Compute number who started S1 training

start_s1_train_ids <- dat$task_log[dat$task_log$session_only == "firstSession" &
                                     dat$task_log$task_name == "Affect" &
                                     dat$task_log$tag == "pre", "participant_id"]
length(start_s1_train_ids) == 1239

# Compute number who did not start S1 training

no_start_s1_ids <- setdiff(stg1_ids, start_s1_train_ids)
length(no_start_s1_ids) == 375

# ---------------------------------------------------------------------------- #
# Identify analysis exclusions ----
# ---------------------------------------------------------------------------- #

# Identify analysis exclusions due to having more than two unique sets of
# DASS-21-AS screening responses by condition (for more information, see
# https://github.com/TeachmanLab/MT-Data-CalmThinkingStudy#participant-flow-and-analysis-exclusions)

exclude_repeat_screen_ids <- 
  dat$participant$participant_id[dat$participant$exclude_analysis == 1]
length(exclude_repeat_screen_ids) == 6

#   Note: Participants 1453 and 1893 did not start Session 1 training so are 
#   already not ITT participants

setdiff(exclude_repeat_screen_ids, start_s1_train_ids) == c(1453, 1893)

# Identify analysis exclusions due to switching conditions, by condition 
# (for more information, see 
# https://github.com/TeachmanLab/MT-Data-CalmThinkingStudy#condition-switching)

exclude_switch_condition_ids <- 382

nrow(dat$study[dat$study$participant_id %in% exclude_switch_condition_ids &
                 dat$study$conditioning == "CONTROL", ]) == 1     # 382

#   Note: Participant 382 would otherwise be part of ITT sample if not excluded   

exclude_switch_condition_ids %in% start_s1_train_ids

# Collect analysis exclusions

exclude_ids <- c(exclude_repeat_screen_ids, exclude_switch_condition_ids)

# ---------------------------------------------------------------------------- #
# Define ITT analysis sample and restrict to ITT participants ----
# ---------------------------------------------------------------------------- #

# Identify intent-to-treat (ITT) analysis sample

itt_anlys_ids <- setdiff(start_s1_train_ids, exclude_ids)
length(itt_anlys_ids) == 1234

# Add indicator for ITT sample to "participant" table

dat$participant$itt_anlys <- 0
dat$participant[dat$participant$participant_id %in% itt_anlys_ids,
                "itt_anlys"] <- 1

# ---------------------------------------------------------------------------- #
# Identify scale items ----
# ---------------------------------------------------------------------------- #

# Identify items for RR, BBSIQ, OASIS, and DASS-21-AS. Note: BBSIQ item names in Calm 
# Thinking are the same as those in Managing Anxiety (stored in "dat_items" list)

rr_neg_thr_items_ct <- sub("_NS", "_ns", dat_items$rr_neg_thr_items_ma)
rr_neg_non_items_ct <- sub("_NF", "_nf", dat_items$rr_neg_non_items_ma)
rr_neg_items_ct <- c(rr_neg_thr_items_ct, rr_neg_non_items_ct)

rr_pos_thr_items_ct <- sub("_PS", "_ps", dat_items$rr_pos_thr_items_ma)
rr_pos_non_items_ct <- sub("_PF", "_pf", dat_items$rr_pos_non_items_ma)
rr_pos_items_ct <- c(rr_pos_thr_items_ct, rr_pos_non_items_ct)

rr_items_ct <- c(rr_neg_items_ct, rr_pos_items_ct)

length(rr_neg_thr_items_ct) == 9
length(rr_neg_non_items_ct) == 9
length(rr_pos_thr_items_ct) == 9
length(rr_pos_non_items_ct) == 9

oa_items_ct <- c("axf", "axs", "avo", "wrk", "soc")

length(oa_items_ct) == 5

dass21_as_items_ct <- c("bre", "dry", "hea", "pan", "sca", "tre", "wor")
dass21_as_mean_items_ct <- paste0(dass21_as_items_ct, "_mean")

length(dass21_as_items_ct) == 7

# Add items to "dat_items" helper list

dat_items$rr_items_ct             <- rr_items_ct
dat_items$rr_neg_items_ct         <- rr_neg_items_ct
dat_items$rr_neg_thr_items_ct     <- rr_neg_thr_items_ct
dat_items$rr_neg_non_items_ct     <- rr_neg_non_items_ct
dat_items$rr_pos_items_ct         <- rr_pos_items_ct
dat_items$rr_pos_thr_items_ct     <- rr_pos_thr_items_ct
dat_items$rr_pos_non_items_ct     <- rr_pos_non_items_ct

dat_items$oa_items_ct             <- oa_items_ct
dat_items$dass21_as_items_ct      <- dass21_as_items_ct
dat_items$dass21_as_mean_items_ct <- dass21_as_mean_items_ct

# ---------------------------------------------------------------------------- #
# Rename columns ----
# ---------------------------------------------------------------------------- #

# Define RR item names to reflect positive/negative, threat/nonthreat, scenario, 
# scenario number (i.e., numeral), and item number (i.e., letter), as was already
# done for Managing Anxiety in prior script

# First, add Calm Thinking items to "rr_item_map" and ensure that they correctly
# map onto Managing Anxiety items

rr_item_map <- rr_item_map[order(rr_item_map$items_ma), ]

rr_item_map$items_ct <- sort(rr_items_ct)

target_cols <- c("items_ma", "items_ct")
rr_item_map <- rr_item_map[, c(target_cols, names(rr_item_map)[!(names(rr_item_map) %in% target_cols)])]

all(tolower(rr_item_map$items_ma) == rr_item_map$items_ct)

rr_item_map <- rr_item_map[order(rr_item_map$valence,
                                 rr_item_map$threat_relevance,
                                 rr_item_map$scenario_number), ]

# Rename RR items

names(dat$rr)[match(rr_items_ct, names(dat$rr))] <- 
  rr_item_map$items_rename[match(rr_items_ct, rr_item_map$items_ct)]

# ---------------------------------------------------------------------------- #
# Handle repeated screenings ----
# ---------------------------------------------------------------------------- #

# Centralized Calm Thinking data cleaning revealed that "dass21_as" table has 
# repeated attempts at screening for some participants. Participants with more 
# than two unique sets of item responses will be excluded from analysis, and 
# "dass21_as_total_anal" (i.e., total score based on row mean of available column 
# means of available multiple item responses) will serve as basis for analysis.

# View(dat$dass21_as[!is.na(dat$dass21_as$n_eligibility_rows) &
#                             dat$dass21_as$n_eligibility_rows > 1, ])

# Remove attempt-specific columns except "X", "id", and date-related columns

dass21_as_exclude_cols <- c("session_and_eligibility_status", "time_on_page",
                            "over18", "dass21_as_total", "dass21_as_total_interp",
                            "dass21_as_eligible")

dat$dass21_as <- dat$dass21_as[, !(names(dat$dass21_as) %in% dass21_as_exclude_cols)]

# Temporarily remove original items

tmp_dass21_as <- dat$dass21_as[, !(names(dat$dass21_as) %in% dass21_as_items_ct)]

# For duplicated rows on "participant_id" and "session_only" when excluding the
# original items above, keep only last row

nrow(tmp_dass21_as) == nrow(dat$dass21_as)

  # TODO: Why does the following line of code leave only one row for participant_id of NA?





# View(dat$dass21_as[is.na(dat$dass21_as$participant_id), ])

dat$dass21_as <- 
  dat$dass21_as[!duplicated(tmp_dass21_as[, c("participant_id", "session_only")], 
                            fromLast = TRUE), ]

# View(dat$dass21_as[is.na(dat$dass21_as$participant_id), ])





# Recode original items as NA at "Eligibility" given that items with "_mean"
# appended are the ones that comprise "dass21_as_total_anal" at that time point

dat$dass21_as[dat$dass21_as$session_only == "Eligibility", dass21_as_items_ct] <- NA

# ---------------------------------------------------------------------------- #
# Handle unexpected multiple entries ----
# ---------------------------------------------------------------------------- #

# Centralized Calm Thinking data cleaning revealed that "bbsiq" and "oa" tables 
# have unexpected multiple entries but identical item responses (just different 
# "time_on_page", so use "time_on_page_mean" for analysis) and that "task_log" 
# table has unexpected multiple entries (some reflecting repeat screenings for 
# "dass21_as" table; so use "time_on_task_mean" for analysis). Note: Centralized
# cleaning did not check "angular_training" table for unexpected multiple entries.

# View(dat$bbsiq[dat$bbsiq$n_rows > 1, ])
# View(dat$oa[dat$oa$n_rows > 1, ])
# View(dat$task_log[dat$task_log$n_rows > 1, ])

# Remove attempt-specific columns except "X", "id", and date-related columns from
# "bbsiq" and "oa" tables

dat$bbsiq$time_on_page <- NULL
dat$oa$time_on_page <- NULL

# For duplicated rows on "participant_id" and "session_only" for "bbsiq" and "oa"
# tables, keep only last row. Leave "task_log" as it is.

dat$bbsiq <- dat$bbsiq[!duplicated(dat$bbsiq[, c("participant_id", "session_only")], 
                                   fromLast = TRUE), ]
dat$oa <-    dat$oa[!duplicated(dat$oa[, c("participant_id", "session_only")], 
                                fromLast = TRUE), ]

# ---------------------------------------------------------------------------- #
# Recode "prefer not to answer" values ----
# ---------------------------------------------------------------------------- #

# Check for outcomes. No relevant tables already have NA values except "dass21_as". 
# For "dass21_as" at "Eligibility", 555 ("prefer not to answer") was already recoded 
# to NA for "_mean" items, and original items were recoded as NA. For "dass21_as" at 
# other time points, "_mean" items are NA (as they were computed only at "Eligibility").

sum(is.na(dat$rr[, dat_items$rr_items])) == 0
sum(is.na(dat$bbsiq[, dat_items$bbsiq_items])) == 0
sum(is.na(dat$oa[, oa_items_ct])) == 0

sum(is.na(dat$dass21_as[dat$dass21_as$session_only != "Eligibility", dass21_as_items_ct]))
sum(is.na(dat$dass21_as[dat$dass21_as$session_only == "Eligibility", dass21_as_mean_items_ct]))

sum(dat$rr[, dat_items$rr_items] == 555, na.rm = TRUE)
sum(dat$bbsiq[, dat_items$bbsiq_items] == 555, na.rm = TRUE)
sum(dat$oa[, oa_items_ct] == 555, na.rm = TRUE)

sum(dat$dass21_as[, dass21_as_items_ct] == 555, na.rm = TRUE)
sum(dat$dass21_as[, dass21_as_mean_items_ct] == 555, na.rm = TRUE)

dat$rr[, dat_items$rr_items][dat$rr[, dat_items$rr_items] == 555] <- NA
dat$bbsiq[, dat_items$bbsiq_items][dat$bbsiq[, dat_items$bbsiq_items] == 555] <- NA
dat$oa[, oa_items_ct][dat$oa[, oa_items_ct] == 555] <- NA

dat$dass21_as[, dass21_as_items_ct][dat$dass21_as[, dass21_as_items_ct] == 555] <- NA

# ---------------------------------------------------------------------------- #
# Check response ranges ----
# ---------------------------------------------------------------------------- #

# Check for outcomes

all(sort(unique(as.vector(as.matrix(dat$rr[, dat_items$rr_items]))))       %in% 1:4)
all(sort(unique(as.vector(as.matrix(dat$bbsiq[, dat_items$bbsiq_items])))) %in% 0:4)
all(sort(unique(as.vector(as.matrix(dat$oa[, oa_items_ct]))))              %in% 0:4)

all(sort(unique(as.vector(as.matrix(dat$dass21_as[, dass21_as_items_ct])))) %in% 0:3)

dass21_as_mean_range <- 
  range(sort(unique(as.vector(as.matrix(dat$dass21_as[, dass21_as_mean_items_ct])))))
all(dass21_as_mean_range >= 0 & dass21_as_mean_range <= 3)

# ---------------------------------------------------------------------------- #
# Restrict to ITT sample and remove participants with no RR data at baseline ----
# ---------------------------------------------------------------------------- #

# Restrict to ITT sample

dat_itt <- lapply(dat, function(x) x[x$participant_id %in% itt_anlys_ids, ])

# 2 participants have no RR items at baseline (and 0 participants are missing BBSIQ,
# DASS-21-AS, and OASIS items), leaving 1232 participants with RR items

nrow(dat_itt$rr[dat_itt$rr$session_only == "preTest" &
                  rowSums(is.na(dat_itt$rr[, dat_items$rr_items])) == length(dat_items$rr_items), ]) == 2
nrow(dat_itt$bbsiq[dat_itt$bbsiq$session_only == "preTest" &
                  rowSums(is.na(dat_itt$bbsiq[, dat_items$bbsiq_items])) == length(dat_items$bbsiq_items), ]) == 0
nrow(dat_itt$dass21_as[dat_itt$dass21_as$session_only == "Eligibility" &
                  rowSums(is.na(dat_itt$dass21_as[, dat_items$dass21_as_mean_items_ct])) == length(dat_items$dass21_as_mean_items_ct), ]) == 0
nrow(dat_itt$oa[dat_itt$oa$session_only == "preTest" &
                  rowSums(is.na(dat_itt$oa[, dat_items$oa_items_ct])) == length(dat_items$oa_items_ct), ]) == 0

length(itt_anlys_ids) - 2 == 1232

ids_missing_bl_rr_items <- dat_itt$rr[dat_itt$rr$session_only == "preTest" &
                                    rowSums(is.na(dat_itt$rr[, dat_items$rr_items])) == length(dat_items$rr_items), "participant_id"]

# Remove participants with no RR item-level data at baseline

dat_ct_rr <- lapply(dat_itt, function(x) x[!(x$participant_id %in% ids_missing_bl_rr_items), ])

# ---------------------------------------------------------------------------- #
# Compute rate of item-level missingness ----
# ---------------------------------------------------------------------------- #

# Compute for RR

round(sum(is.na(dat_ct_rr$rr[, dat_items$rr_items])) / (length(dat_items$rr_items) * 1232), 4) == 0.0141

# ---------------------------------------------------------------------------- #
# Export data ----
# ---------------------------------------------------------------------------- #

save(dat_items, file = "./data/helper/dat_items.RData")
save(rr_item_map, file = "./data/helper/rr_item_map.RData")

dir.create("./data/ct/intermediate_clean_further")

save(dat_ct_rr, file = "./data/ct/intermediate_clean_further/dat_ct_rr.RData")