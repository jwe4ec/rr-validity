# ---------------------------------------------------------------------------- #
# Prepare Managing Anxiety Data
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# The present script imports a few potential sources of data for the Managing Anxiety
# study and determines which source(s) to use

# 1. Clean data files ("R34_Cronbach.csv" and "R34_FinalData_New_v02.csv") from the 
#    OSF project (https://osf.io/3b67v) for the Managing Anxiety Study main outcomes 
#    paper (https://doi.org/g62s). "R34_Cronbach.csv" has item-level data at baseline, 
#    whereas "R34_FinalData_New_v02.csv" has scale-level data at multiple time points.

#    "R34_FinalData_New_v02.csv" is generated from "FinalData-28Feb20_v02.csv" using 
#    "Script1_DataPrep.R" from the OSF project. "FinalData-28Feb20_v02.csv", which is 
#    not on the OSF project, is a newer version of the data file "FinalData-28Feb20.csv" 
#    outputted by the script "R34.ipynb" in the study's Data Cleaning folder 
#    (https://bit.ly/3CLi5It) on GitHub. Sonia Baee indicated that the files in the 
#    Data Cleaning folder and the Main Outcomes folder (https://bit.ly/3FHRz4G) on 
#    GitHub are out of date given file losses when switching laptops. Code to generate 
#    "FinalData-28Feb20_v02.csv" is not presently available as a result.

#    Code used to generate "R34_Cronbach.csv", likely written by Sonia Baee, is not
#    available. However, the file was likely created by adjusting the script "R34.ipynb" 
#    in the study's Data Cleaning folder (https://bit.ly/3CLi5It) on GitHub.

# 2. Raw data files obtained from Sonia Baee on 9/3/2020, who stated on that date that
#    they represent the latest version of the database on the R34 server and that she 
#    obtained them from Claudia Calicho-Mamani

# 3. Raw data files from the OSF project for the Managing Anxiety study, which Yuhan
#    Hou uploaded on 3/16/2023 after dumping the CSV files from the "teachmanlab" SQL
#    Data Server using Grafana (because the "angular_training" table was too large to
#    obtain via Grafana, Yuhan dumped it directly from the Data Server).

# Although the present script mentions Calm Thinking analyses, we decided not to
# analyze Calm Thinking data as part of the dissertation's scope

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
# Import clean data files from main outcomes paper ----
# ---------------------------------------------------------------------------- #

# Import baseline items and longitudinal scales from OSF project for Managing Anxiety 
# main outcomes paper (see Source 1 above for more info)

dat_main_bl_items  <- read.csv("./data/ma/source/clean_from_main_paper/R34_Cronbach.csv")
dat_main_lg_scales <- read.csv("./data/ma/source/clean_from_main_paper/R34_FinalData_New_v02.csv")

# ---------------------------------------------------------------------------- #
# Import raw data files from Sonia and from MA study OSF project ----
# ---------------------------------------------------------------------------- #

# Obtain file names of raw CSV data files from Sonia (see Source 2 above) and
# from MA study OSF project (see Source 3 above)

raw_data_dir_sb  <- paste0(wd_dir, "/data/ma/source/raw_from_sonia")
raw_data_dir_osf <- paste0(wd_dir, "/data/ma/source/raw_from_osf")

raw_filenames_sb  <- list.files(raw_data_dir_sb, 
                                pattern = "*.csv", full.names = FALSE)
raw_filenames_osf <- list.files(raw_data_dir_osf, 
                                pattern = "*.csv", full.names = FALSE)

# Import raw data files and store them in list

raw_dat_sb  <- lapply(paste0(raw_data_dir_sb,  "/", raw_filenames_sb),  read.csv)
raw_dat_osf <- lapply(paste0(raw_data_dir_osf, "/", raw_filenames_osf), read.csv)

# Name each raw data file in list

split_char <- ".csv"

names(raw_dat_sb)  <- unlist(lapply(raw_filenames_sb,
                                    function(f) {
                                      unlist(strsplit(f,
                                                      split = split_char,
                                                      fixed = FALSE))[1]
                                    }))
names(raw_dat_osf) <- unlist(lapply(raw_filenames_osf,
                                    function(f) {
                                      unlist(strsplit(f,
                                                      split = split_char,
                                                      fixed = FALSE))[1]
                                    }))

# ---------------------------------------------------------------------------- #
# Decide which data to use for measurement analyses ----
# ---------------------------------------------------------------------------- #

# For measurement analyses, use baseline items from OSF project for Managing Anxiety 
# main outcomes paper (Source 1 above), given that this presumably is the cleanest
# file with item-level data for intent-to-treat participants

dat <- dat_main_bl_items

# Similarly, use demographics data from longitudinal scales at Source 1

# ---------------------------------------------------------------------------- #
# Identify scale items ----
# ---------------------------------------------------------------------------- #

# Identify items for RR, BBSIQ, OASIS, and DASS-21-AS. Note: Names of items for
# RR, OASIS, and DASS-21-AS are different in the Managing Anxiety study versus
# the Calm Thinking study

rr_scenarios <- c("blood_test", "elevator", "job", "lunch", "meeting_friend", 
                  "noise", "scrape", "shopping", "wedding")

rr_neg_thr_items_ma <- paste0(rr_scenarios, "_NS")
rr_neg_non_items_ma <- paste0(rr_scenarios, "_NF")
rr_neg_items_ma <- c(rr_neg_thr_items_ma, rr_neg_non_items_ma)

rr_pos_thr_items_ma <- paste0(rr_scenarios, "_PS")
rr_pos_non_items_ma <- paste0(rr_scenarios, "_PF")
rr_pos_items_ma <- c(rr_pos_thr_items_ma, rr_pos_non_items_ma)

rr_items_ma <- c(rr_neg_items_ma, rr_pos_items_ma)

length(rr_neg_thr_items_ma) == 9
length(rr_neg_non_items_ma) == 9
length(rr_pos_thr_items_ma) == 9
length(rr_pos_non_items_ma) == 9

bbsiq_neg_int_items <- c("breath_suffocate", "vision_illness", "lightheaded_faint", "chest_heart", 
                         "heart_wrong", "confused_outofmind", "dizzy_ill")
bbsiq_neg_ext_items <- c("visitors_bored", "shop_irritating", "smoke_house", "friend_incompetent", 
                         "jolt_burglar", "party_boring", "urgent_died")
bbsiq_neg_items <- c(bbsiq_neg_int_items, bbsiq_neg_ext_items)

bbsiq_ben_int_items <- c("breath_flu", "breath_physically", "vision_glasses", "vision_strained",
                         "lightheaded_eat", "lightheaded_sleep", "chest_indigestion", "chest_sore",
                         "heart_active", "heart_excited", "confused_cold", "confused_work",
                         "dizzy_ate", "dizzy_overtired")
bbsiq_ben_ext_items <- c("visitors_engagement", "visitors_outstay", "shop_bored", "shop_concentrating",
                         "smoke_cig", "smoke_food", "friend_helpful", "friend_moreoften", "jolt_dream",
                         "jolt_wind", "party_hear", "party_preoccupied", "urgent_bill", "urgent_junk")
bbsiq_ben_items <- c(bbsiq_ben_int_items, bbsiq_ben_ext_items)

bbsiq_items <- c(bbsiq_neg_items, bbsiq_ben_items)

length(bbsiq_neg_int_items) == 7
length(bbsiq_neg_ext_items) == 7
length(bbsiq_neg_items)     == 14

length(bbsiq_ben_int_items) == 14
length(bbsiq_ben_ext_items) == 14
length(bbsiq_ben_items)     == 28 

oa_items_ma <- c("OA_anxious_freq", "OA_anxious_sev", "OA_avoid", "OA_interfere", "OA_interfere_social")

length(oa_items_ma) == 5

dass21_as_items_ma <- c("DASSA_breathing", "DASSA_dryness", "DASSA_heart", "DASSA_panic", "DASSA_scared", 
                     "DASSA_trembling", "DASSA_worry")

length(dass21_as_items_ma) == 7

# ---------------------------------------------------------------------------- #
# Rename columns ----
# ---------------------------------------------------------------------------- #

# Rename participant ID

names(dat)[names(dat) == "participantID"] <- "participant_id"

# Define RR item names to reflect positive/negative, threat/nonthreat, scenario, 
# scenario number (i.e., numeral), and item number (i.e., letter). These renamed
# items will be used in both Managing Anxiety and Calm Thinking analyses.

rr_items_rename <-
  c("neg_thr_blood_test_9c", "neg_thr_elevator_1a", "neg_thr_job_3c", "neg_thr_lunch_6a", "neg_thr_meeting_friend_5d", 
    "neg_thr_noise_4b", "neg_thr_scrape_7a", "neg_thr_shopping_8a", "neg_thr_wedding_2c",
    "neg_non_blood_test_9b", "neg_non_elevator_1c", "neg_non_job_3d", "neg_non_lunch_6b", "neg_non_meeting_friend_5b", 
    "neg_non_noise_4d", "neg_non_scrape_7b", "neg_non_shopping_8d", "neg_non_wedding_2d", 
    "pos_thr_blood_test_9d", "pos_thr_elevator_1d", "pos_thr_job_3a", "pos_thr_lunch_6c", "pos_thr_meeting_friend_5c", 
    "pos_thr_noise_4c", "pos_thr_scrape_7c", "pos_thr_shopping_8b", "pos_thr_wedding_2a", 
    "pos_non_blood_test_9a", "pos_non_elevator_1b", "pos_non_job_3b", "pos_non_lunch_6d", "pos_non_meeting_friend_5a", 
    "pos_non_noise_4a", "pos_non_scrape_7d", "pos_non_shopping_8c", "pos_non_wedding_2b")

valence <- c(rep("neg", 18),
             rep("pos", 18))

threat_relevance <- c(rep("thr", 9),
                      rep("non", 9),
                      rep("thr", 9),
                      rep("non", 9))

scenario_number <- rep(c(9, 1, 3, 6, 5, 
                         4, 7, 8, 2), 4)

rr_item_map <- data.frame(items_ma         = rr_items_ma,
                          items_rename     = rr_items_rename,
                          valence          = valence,
                          threat_relevance = threat_relevance,
                          scenario         = rr_scenarios,
                          scenario_number  = scenario_number)

rr_item_map <- rr_item_map[order(rr_item_map$valence,
                                 rr_item_map$threat_relevance,
                                 rr_item_map$scenario_number), ]

# Rename RR items

names(dat)[match(rr_items_ma, names(dat))] <- 
  rr_item_map$items_rename[match(rr_items_ma, rr_item_map$items_ma)]

rr_neg_thr_items <- rr_item_map$items_rename[rr_item_map$valence == "neg" &
                                               rr_item_map$threat_relevance == "thr"]
rr_neg_non_items <- rr_item_map$items_rename[rr_item_map$valence == "neg" &
                                               rr_item_map$threat_relevance == "non"]
rr_pos_thr_items <- rr_item_map$items_rename[rr_item_map$valence == "pos" &
                                               rr_item_map$threat_relevance == "thr"]
rr_pos_non_items <- rr_item_map$items_rename[rr_item_map$valence == "pos" &
                                               rr_item_map$threat_relevance == "non"]

rr_neg_items <- c(rr_neg_thr_items, rr_neg_non_items)
rr_pos_items <- c(rr_pos_thr_items, rr_pos_non_items)
rr_items     <- c(rr_neg_items, rr_pos_items)

# Store RR items in list

dat_items <- list(rr_items_ma         = rr_items_ma,
                  rr_neg_items_ma     = rr_neg_items_ma,
                  rr_neg_thr_items_ma = rr_neg_thr_items_ma,
                  rr_neg_non_items_ma = rr_neg_non_items_ma,
                  rr_pos_items_ma     = rr_pos_items_ma,
                  rr_pos_thr_items_ma = rr_pos_thr_items_ma,
                  rr_pos_non_items_ma = rr_pos_non_items_ma,
                  rr_items            = rr_items,
                  rr_neg_items        = rr_neg_items,
                  rr_neg_thr_items    = rr_neg_thr_items,
                  rr_neg_non_items    = rr_neg_non_items,
                  rr_pos_items        = rr_pos_items,
                  rr_pos_thr_items    = rr_pos_thr_items,
                  rr_pos_non_items    = rr_pos_non_items,
                  bbsiq_items         = bbsiq_items,
                  oa_items_ma         = oa_items_ma,
                  dass21_as_items_ma  = dass21_as_items_ma)

# ---------------------------------------------------------------------------- #
# Confirm ITT sample size and lack of unexpected multiple entries ----
# ---------------------------------------------------------------------------- #

length(dat$participant_id) == 807
length(unique(dat$participant_id)) == 807

# ---------------------------------------------------------------------------- #
# Recode "prefer not to answer" values ----
# ---------------------------------------------------------------------------- #

# Recode "prefer not to answer" (coded as 555) as NA

target_items <- c(dat_items$rr_items,
                  dat_items$bbsiq_items,
                  dat_items$oa_items_ma,
                  dat_items$dass21_as_items_ma)

dat[, target_items][dat[, target_items] == 555] <- NA

# ---------------------------------------------------------------------------- #
# Check response ranges ----
# ---------------------------------------------------------------------------- #

# Note: In present Managing Anxiety study, RR response options were 0:3, whereas 
# in Calm Thinking study, they were 1:4

all(sort(unique(as.vector(as.matrix(dat[, rr_items]))))           %in% 0:3)
all(sort(unique(as.vector(as.matrix(dat[, bbsiq_items]))))        %in% 0:4)
all(sort(unique(as.vector(as.matrix(dat[, oa_items_ma]))))        %in% 0:4)
all(sort(unique(as.vector(as.matrix(dat[, dass21_as_items_ma])))) %in% 0:3)

# ---------------------------------------------------------------------------- #
# Remove participants with no RR data at baseline ----
# ---------------------------------------------------------------------------- #

# 58 participants have no RR items at baseline (and 8, 9, and 2 are missing BBSIQ,
# DASS-21-AS, and OASIS items, respectively), leaving 749 participants with RR items

nrow(dat[rowSums(is.na(dat[, dat_items$rr_items]))           == length(dat_items$rr_items), ])           == 58
nrow(dat[rowSums(is.na(dat[, dat_items$bbsiq_items]))        == length(dat_items$bbsiq_items), ])        == 8
nrow(dat[rowSums(is.na(dat[, dat_items$dass21_as_items_ma])) == length(dat_items$dass21_as_items_ma), ]) == 9
nrow(dat[rowSums(is.na(dat[, dat_items$oa_items_ma]))        == length(dat_items$oa_items_ma), ])        == 2

807 - 58 == 749

ids_missing_bl_rr_items <- dat[rowSums(is.na(dat[, dat_items$rr_items])) == 
                                 length(dat_items$rr_items), "participant_id"]

# The other clean data file with longitudinal scales from OSF project for Managing
# Anxiety main outcomes paper (see Source 1 above) similarly has RR scores for only
# 748 participants at baseline

sum(!is.na(dat_main_lg_scales$RR_negative_nf_score_PRE)) == 748
sum(!is.na(dat_main_lg_scales$RR_negative_ns_score_PRE)) == 748
sum(!is.na(dat_main_lg_scales$RR_positive_pf_score_PRE)) == 748
sum(!is.na(dat_main_lg_scales$RR_positive_ps_score_PRE)) == 748

sum(!is.na(dat_main_lg_scales$bbsiq_physical_score_PRE)) == 797
sum(!is.na(dat_main_lg_scales$bbsiq_threat_score_PRE))   == 797
sum(!is.na(dat_main_lg_scales$negativeBBSIQ_PRE))        == 797

# Moreover, additional baseline RR and BBSIQ data are not in raw data from Sonia or 
# in raw data from OSF project for Managing Anxiety study (see Sources 2 and 3 above)

# Remove participants with no RR item-level data at baseline

dat_ma_rr <- dat[!(dat$participant_id %in% ids_missing_bl_rr_items), ]

# ---------------------------------------------------------------------------- #
# Compute rate of item-level missingness ----
# ---------------------------------------------------------------------------- #

# Compute for RR

round(sum(is.na(dat_ma_rr[, rr_items])) / (length(rr_items) * 749), 4) == .0026

# ---------------------------------------------------------------------------- #
# Extract demographics data  ----
# ---------------------------------------------------------------------------- #

# Extract demographics data from other clean data file with longitudinal scales
# (see Source 1 above)

  # Rename participant ID

names(dat_main_lg_scales)[names(dat_main_lg_scales) == "participantID"] <- "participant_id"

  # Remove participants with no RR item-level data at baseline

dat_ma_dem <- dat_main_lg_scales[!(dat_main_lg_scales$participant_id %in% ids_missing_bl_rr_items), ]

  # Confirm expected number of participants

length(unique(dat_ma_dem$participant_id)) == 749
all(dat_ma_dem$participant_id %in% dat_ma_rr$participant_id)

  # Restrict to demographic columns

dem_cols <- c("demographic_age", "demographic_birthYear", "demographic_residenceCountry",
              "demographic_gender", "demographic_race", "demographic_ethnicity", 
              "demographic_maritalStatus", "demographic_education", "demographic_employmentStatus",
              "demographic_income")

dat_ma_dem <- dat_ma_dem[, c("participant_id", dem_cols)]

  # Rename demographic columns

names(dat_ma_dem) <- sub("demographic_", "", names(dat_ma_dem))

# Extract raw demographics data to investigate "birthYear" in "clean_ma_demog_data_create_tbl.R"

raw_dat_sb_dem  <- raw_dat_sb$Demographic_recovered_Feb_02_2019
raw_dat_osf_dem <- raw_dat_osf$Demographics_02_02_2019

# ---------------------------------------------------------------------------- #
# Export data ----
# ---------------------------------------------------------------------------- #

dir.create("./data/helper")

save(dat_items,   file = "./data/helper/dat_items.RData")
save(rr_item_map, file = "./data/helper/rr_item_map.RData")

dir.create("./data/ma/intermediate")

save(dat_ma_rr,  file = "./data/ma/intermediate/dat_ma_rr.RData")
save(dat_ma_dem, file = "./data/ma/intermediate/dat_ma_dem.RData")

# Export raw demographics tables

save(raw_dat_sb_dem,  file = "./data/ma/intermediate/raw_dat_sb_dem.RData")
save(raw_dat_osf_dem, file = "./data/ma/intermediate/raw_dat_osf_dem.RData")