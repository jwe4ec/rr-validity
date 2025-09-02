# ---------------------------------------------------------------------------- #
# Create Tables and Figures
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

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

groundhog.library("lavaan", groundhog_day)

# ---------------------------------------------------------------------------- #
# Import results ----
# ---------------------------------------------------------------------------- #

load("./results/cfa/initial/res_ma_cfa_initial.RData")
load("./results/cfa/revised/res_ma_cfa_28.RData")
load("./results/cfa/revised/res_ma_cfa_thr.RData")

load("./results/efa/res_ma_efa.RData")

# ---------------------------------------------------------------------------- #
# Create CFA fit statistics tables ----
# ---------------------------------------------------------------------------- #

# Extract fit statistics for key models (use NA for improper solutions)

tar_stats <- c("chisq.scaled", "df.scaled", "pvalue.scaled", 
               "srmr", "rmsea.robust", "cfi.robust", "tli.robust")

ma_cfa_cf_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_initial$ma_cfa_cf_wlsmv_fit, tar_stats)
ma_cfa_bf_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_initial$ma_cfa_bf_wlsmv_fit, tar_stats)
ma_cfa_ho_wlsmv_fit_stats <- NA

tbl_fit_ma_cfa_initial <- as.data.frame(rbind(ma_cfa_cf_wlsmv_fit_stats, ma_cfa_bf_wlsmv_fit_stats, 
                                              ma_cfa_ho_wlsmv_fit_stats))

ma_cfa_28_3cf_efa_wlsmv_fit_stats    <- fitMeasures(res_ma_cfa_28$ma_cfa_28_3cf_efa_wlsmv_fit,    tar_stats)
ma_cfa_28_3cf_wlsmv_fit_stats        <- fitMeasures(res_ma_cfa_28$ma_cfa_28_3cf_wlsmv_fit,        tar_stats)
ma_cfa_28_2cf_wlsmv_fit_stats        <- fitMeasures(res_ma_cfa_28$ma_cfa_28_2cf_wlsmv_fit,        tar_stats)
ma_cfa_28_3cf_ce_thr_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_28$ma_cfa_28_3cf_ce_thr_wlsmv_fit, tar_stats)
ma_cfa_28_3cf_ce_all_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_28$ma_cfa_28_3cf_ce_all_wlsmv_fit, tar_stats)

tbl_fit_ma_cfa_28 <- as.data.frame(rbind(ma_cfa_28_3cf_efa_wlsmv_fit_stats, ma_cfa_28_3cf_wlsmv_fit_stats,
                                         ma_cfa_28_2cf_wlsmv_fit_stats, ma_cfa_28_3cf_ce_thr_wlsmv_fit_stats,
                                         ma_cfa_28_3cf_ce_all_wlsmv_fit_stats))

ma_cfa_thr_2cf_efa_wlsmv_fit_stats            <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_2cf_efa_wlsmv_fit,            tar_stats)
ma_cfa_thr_2cf_wlsmv_fit_stats                <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_2cf_wlsmv_fit,                tar_stats)
ma_cfa_thr_2cf_ce_wlsmv_fit_stats             <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_2cf_ce_wlsmv_fit,             tar_stats)
ma_cfa_thr_1f_ce_scen_wlsmv_fit_stats         <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_ce_scen_wlsmv_fit,         tar_stats)
ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit_stats    <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit,    tar_stats)
ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit, tar_stats)
ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit, tar_stats)
ma_cfa_thr_bf_ce_scen_wlsmv_fit_stats         <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_bf_ce_scen_wlsmv_fit,         tar_stats)
ma_cfa_thr_ho_ce_wlsmv_fit_stats              <- NA

tbl_fit_ma_cfa_thr <- as.data.frame(rbind(ma_cfa_thr_2cf_efa_wlsmv_fit_stats, ma_cfa_thr_2cf_wlsmv_fit_stats,
                                          ma_cfa_thr_2cf_ce_wlsmv_fit_stats, ma_cfa_thr_1f_ce_scen_wlsmv_fit_stats,
                                          ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit_stats, ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit_stats,
                                          ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit_stats, ma_cfa_thr_bf_ce_scen_wlsmv_fit_stats,
                                          ma_cfa_thr_ho_ce_wlsmv_fit_stats))

# Extract fit statistics for supplemental models (use NA for improper solutions)

ma_cfa_28_2cf_4cmf_wlsmv_fit_stats                        <- NA
ma_cfa_28_2cf_2cmf_thr_wlsmv_fit_stats                    <- NA
ma_cfa_28_2cf_2mf_thr_wlsmv_fit_stats                     <- NA
ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit_stats     <- NA
ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_28$ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit, tar_stats)

tbl_fit_ma_cfa_28_suppl <- as.data.frame(rbind(ma_cfa_28_2cf_4cmf_wlsmv_fit_stats, ma_cfa_28_2cf_2cmf_thr_wlsmv_fit_stats, 
                                               ma_cfa_28_2cf_2mf_thr_wlsmv_fit_stats, ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit_stats,
                                               ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit_stats))

ma_cfa_thr_1f_ce_val_scen_wlsmv_fit_stats <- NA
ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit, tar_stats)
ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit_stats <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit, tar_stats)
ma_cfa_thr_1f_ce_val_wlsmv_fit_stats      <- NA
ma_cfa_thr_1f_2cmf_wlsmv_fit_stats        <- fitMeasures(res_ma_cfa_thr$ma_cfa_thr_1f_2cmf_wlsmv_fit,        tar_stats)

tbl_fit_ma_cfa_thr_suppl <- as.data.frame(rbind(ma_cfa_thr_1f_ce_val_scen_wlsmv_fit_stats, ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit_stats, 
                                                ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit_stats, ma_cfa_thr_1f_ce_val_wlsmv_fit_stats,
                                                ma_cfa_thr_1f_2cmf_wlsmv_fit_stats))

# Define function to format tables

format_fit_tbl <- function(df) {
  df$model <- row.names(df)
  df$model <- sub("_fit_stats", "", df$model)
  
  row.names(df) <- 1:nrow(df)
  
  df <- df[c("model", names(df)[names(df) != "model"])]
  
  df$chisq.scaled <- round(df$chisq.scaled, 2)
  
  tar_cols <- c("srmr", "rmsea.robust", "cfi.robust", "tli.robust")
  df[tar_cols] <- round(df[tar_cols], 3)
  
  return(df)
}

# Run function

tbl_fit_ma_cfa_initial   <- format_fit_tbl(tbl_fit_ma_cfa_initial)
tbl_fit_ma_cfa_28        <- format_fit_tbl(tbl_fit_ma_cfa_28)
tbl_fit_ma_cfa_thr       <- format_fit_tbl(tbl_fit_ma_cfa_thr)

tbl_fit_ma_cfa_28_suppl  <- format_fit_tbl(tbl_fit_ma_cfa_28_suppl)
tbl_fit_ma_cfa_thr_suppl <- format_fit_tbl(tbl_fit_ma_cfa_thr_suppl)

# Export tables

cfa_tbls_path <- "./results/cfa/tables/"

dir.create(cfa_tbls_path, recursive = TRUE)

write.csv(tbl_fit_ma_cfa_initial,   paste0(cfa_tbls_path, "tbl_fit_ma_cfa_initial.csv"),   row.names = FALSE)
write.csv(tbl_fit_ma_cfa_28,        paste0(cfa_tbls_path, "tbl_fit_ma_cfa_28.csv"),        row.names = FALSE)
write.csv(tbl_fit_ma_cfa_thr,       paste0(cfa_tbls_path, "tbl_fit_ma_cfa_thr.csv"),       row.names = FALSE)

write.csv(tbl_fit_ma_cfa_28_suppl,  paste0(cfa_tbls_path, "tbl_fit_ma_cfa_28_suppl.csv"),  row.names = FALSE)
write.csv(tbl_fit_ma_cfa_thr_suppl, paste0(cfa_tbls_path, "tbl_fit_ma_cfa_thr_suppl.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Create EFA fit statistics table ----
# ---------------------------------------------------------------------------- #

# Extract fit statistics for key models

tar_stats <- c("chisq.scaled", "df.scaled", "pvalue.scaled", 
               "rmsea.robust", "cfi.robust")

ma_efa_oblimin_wlsmv_fit_stats     <- t(fitMeasures(res_ma_efa$ma_efa_oblimin_wlsmv_fit,     tar_stats))
ma_efa_35_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_35_oblimin_wlsmv_fit,  tar_stats))
ma_efa_34_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_34_oblimin_wlsmv_fit,  tar_stats))
ma_efa_33_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_33_oblimin_wlsmv_fit,  tar_stats))
ma_efa_32_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_32_oblimin_wlsmv_fit,  tar_stats))
ma_efa_31_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_31_oblimin_wlsmv_fit,  tar_stats))
ma_efa_30_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_30_oblimin_wlsmv_fit,  tar_stats))
ma_efa_29_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_29_oblimin_wlsmv_fit,  tar_stats))
ma_efa_28_oblimin_wlsmv_fit_stats  <- t(fitMeasures(res_ma_efa$ma_efa_28_oblimin_wlsmv_fit,  tar_stats))
ma_efa_thr_oblimin_wlsmv_fit_stats <- t(fitMeasures(res_ma_efa$ma_efa_thr_oblimin_wlsmv_fit, tar_stats))

tbl_fit_ma_efa <- as.data.frame(rbind(ma_efa_oblimin_wlsmv_fit_stats, ma_efa_35_oblimin_wlsmv_fit_stats, 
                                      ma_efa_34_oblimin_wlsmv_fit_stats, ma_efa_33_oblimin_wlsmv_fit_stats, 
                                      ma_efa_32_oblimin_wlsmv_fit_stats, ma_efa_31_oblimin_wlsmv_fit_stats, 
                                      ma_efa_30_oblimin_wlsmv_fit_stats, ma_efa_29_oblimin_wlsmv_fit_stats, 
                                      ma_efa_28_oblimin_wlsmv_fit_stats, ma_efa_thr_oblimin_wlsmv_fit_stats))

tbl_fit_ma_efa$model <- c(rep("ma_efa_oblimin_wlsmv_fit_stats",     nrow(ma_efa_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_35_oblimin_wlsmv_fit_stats",  nrow(ma_efa_35_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_34_oblimin_wlsmv_fit_stats",  nrow(ma_efa_34_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_33_oblimin_wlsmv_fit_stats",  nrow(ma_efa_33_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_32_oblimin_wlsmv_fit_stats",  nrow(ma_efa_32_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_31_oblimin_wlsmv_fit_stats",  nrow(ma_efa_31_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_30_oblimin_wlsmv_fit_stats",  nrow(ma_efa_30_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_29_oblimin_wlsmv_fit_stats",  nrow(ma_efa_29_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_28_oblimin_wlsmv_fit_stats",  nrow(ma_efa_28_oblimin_wlsmv_fit_stats)),
                          rep("ma_efa_thr_oblimin_wlsmv_fit_stats", nrow(ma_efa_thr_oblimin_wlsmv_fit_stats)))

# Format table

tbl_fit_ma_efa$nfactors <- row.names(tbl_fit_ma_efa)
tbl_fit_ma_efa$nfactors <- sub("nfactors...", "", tbl_fit_ma_efa$nfactors)
tbl_fit_ma_efa$nfactors <- gsub("\\..*$", "", tbl_fit_ma_efa$nfactors)

idx_vars <- c("model", "nfactors")

tbl_fit_ma_efa <- tbl_fit_ma_efa[c(idx_vars, names(tbl_fit_ma_efa)[!(names(tbl_fit_ma_efa) %in% idx_vars)])]

row.names(tbl_fit_ma_efa) <- 1:nrow(tbl_fit_ma_efa)

tbl_fit_ma_efa$chisq.scaled <- round(tbl_fit_ma_efa$chisq.scaled, 2)

tar_cols <- c("rmsea.robust", "cfi.robust")

tbl_fit_ma_efa[tar_cols] <- round(tbl_fit_ma_efa[tar_cols], 3)

# Export table

efa_tbl_path <- "./results/efa/table/"

dir.create(efa_tbl_path, recursive = TRUE)

write.csv(tbl_fit_ma_efa, paste0(efa_tbl_path, "tbl_fit_ma_efa.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# TODO: Write CFA figures to MS Word ----
# ---------------------------------------------------------------------------- #

# "officer" package does not seem to be able to add PDF files to Word and had an 
# issue adding "qgraph"s to Word. Also, "qgraph" PNGs seem to differ from PDFs. So
# it's currently unclear how to write "qgraph"s to Word accurately.





