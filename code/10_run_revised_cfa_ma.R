# ---------------------------------------------------------------------------- #
# Run Revised CFA in Managing Anxiety Study
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# Note: Random seed must be set for each CFA analysis for reproducible results

# ---------------------------------------------------------------------------- #
# Store working directory, check correct R version, load packages ----
# ---------------------------------------------------------------------------- #

# Store working directory

wd_dir <- getwd()

# Load custom functions

source("./code/01a_define_functions.R")

# Check correct R version, load groundhog package, and specify groundhog_day

groundhog_day <- version_control()

# Load packages (note: "dynamic" requires "devtools" and "GenOrd")

pkgs <- c("lavaan", "lavaanExtra", "semPlot", "qgraph", "corrplot", "devtools", "GenOrd")
groundhog.library(pkgs, groundhog_day)

  # TODO (install package from CRAN once updated): Categorical data functions ("catOne",
  # "catHB") for "dynamic" package are not yet on CRAN; install from latest GitHub commit 
  # "bf9a430" ("likertHB2 minor update") as of 1/24/24 (for documentation on functions, 
  # see McNeish, 2023, preprint, https://doi.org/10.31234/osf.io/tp35s; and vignette 
  # https://rpubs.com/dmcneish/1025400)

devtools::install_github("melissagwolf/dynamic", ref = "bf9a430") # Skip prompt to update packages

# Allow printing of more lines

getOption("max.print")
options(max.print=10000)

# ---------------------------------------------------------------------------- #
# Import data and create CFA revised results path ----
# ---------------------------------------------------------------------------- #

load("./data/ma/final/dat_ma_rr_only.RData")

load("./data/helper/dat_items.RData")
load("./data/helper/rr_item_map.RData")

cfa_revised_path <- "./results/cfa/revised/"
dir.create(cfa_revised_path)

# ---------------------------------------------------------------------------- #
# Define revised scale items ----
# ---------------------------------------------------------------------------- #

# Based on EFA results in "run_efa_ma.R", define revised scales based on 28 items

rr_28_exclude_items <- c("pos_thr_wedding_2a", "neg_thr_scrape_7a", "pos_thr_elevator_1d",
                         "pos_thr_job_3a", "pos_thr_meeting_friend_5c", "pos_non_scrape_7d",
                         "pos_non_lunch_6d", "pos_non_blood_test_9a")

length(rr_28_exclude_items) == 8

dat_items$rr_28_neg_thr_items <- setdiff(dat_items$rr_neg_thr_items, rr_28_exclude_items)
dat_items$rr_28_neg_non_items <- setdiff(dat_items$rr_neg_non_items, rr_28_exclude_items)
dat_items$rr_28_pos_thr_items <- setdiff(dat_items$rr_pos_thr_items, rr_28_exclude_items)
dat_items$rr_28_pos_non_items <- setdiff(dat_items$rr_pos_non_items, rr_28_exclude_items)

length(dat_items$rr_28_neg_thr_items) == 8
length(dat_items$rr_28_neg_non_items) == 9
length(dat_items$rr_28_pos_thr_items) == 5
length(dat_items$rr_28_pos_non_items) == 6

dat_items$rr_28_neg_items <- c(dat_items$rr_28_neg_thr_items, dat_items$rr_28_neg_non_items)
dat_items$rr_28_pos_items <- c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_pos_non_items)
dat_items$rr_28_items     <- c(dat_items$rr_28_neg_items, dat_items$rr_28_pos_items)

length(dat_items$rr_28_neg_items) == 17
length(dat_items$rr_28_pos_items) == 11
length(dat_items$rr_28_items)     == 28

# Save items

save(dat_items, file = "./data/helper/dat_items.RData")

# ---------------------------------------------------------------------------- #
# Run CFA for revised 3 correlated-factors model based on 28 items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_28_pos_thr_items,
            neg_thr = dat_items$rr_28_neg_thr_items,
            non     = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items))

ma_cfa_28_3cf <- write_lavaan(latent = lat)
cat(ma_cfa_28_3cf)

# ---------------------------------------------------------------------------- #
# Run using PML-AC ----
# ---------------------------------------------------------------------------- #

# Note: Negative thresholds are OK

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_3cf_pml_ac.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_3cf_pml_ac_fit <- cfa(ma_cfa_28_3cf, data = dat_ma_rr_only, std.lv = TRUE,
                                ordered = TRUE, estimator = "PML", missing = "available.cases")
(ma_cfa_28_3cf_pml_ac_summ <- summary(ma_cfa_28_3cf_pml_ac_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_28_3cf_pml_ac_std_all <- standardizedsolution(ma_cfa_28_3cf_pml_ac_fit))
sink(type = "message")
sink()

# TODO: No Heywood cases, but note that some fit statistics may be off





# Inspect residual correlations (focus on at least small ones; i.e., > .1; Flora, 2020, p. 490)

ma_cfa_28_3cf_pml_ac_resid_cor <- round(residuals(ma_cfa_28_3cf_pml_ac_fit, type = "cor")$cov, 2)

ma_cfa_28_3cf_pml_ac_resid_cor_small <- ma_cfa_28_3cf_pml_ac_resid_cor
ma_cfa_28_3cf_pml_ac_resid_cor_small[abs(ma_cfa_28_3cf_pml_ac_resid_cor_small) < .1 ] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_pml_ac_resid_cor_small_plot.pdf"))
ma_cfa_28_3cf_pml_ac_resid_cor_small_p <- 
  corrplot(ma_cfa_28_3cf_pml_ac_resid_cor_small, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey",
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Plot completely standardized solution

ma_cfa_28_3cf_pml_ac_p <- 
  semPaths(ma_cfa_28_3cf_pml_ac_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "neg_thr", "pos_thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 30, 1, 12),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_3cf_pml_ac_plot")

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_3cf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_3cf_wlsmv_fit <- cfa(ma_cfa_28_3cf, data = dat_ma_rr_only, std.lv = TRUE,
                               ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_3cf_wlsmv_summ <- summary(ma_cfa_28_3cf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_3cf_wlsmv_std_all <- standardizedsolution(ma_cfa_28_3cf_wlsmv_fit))
cat("\n", "DYNAMIC FIT INDEX CUTOFFS:", "\n\n")
(ma_cfa_28_3cf_wlsmv_dfi <- dynamic::catHB(ma_cfa_28_3cf_wlsmv_fit, plot = TRUE, 
                                           estimator = "WLSMV", reps = 250))
sink(type = "message")
sink()

# Save dynamic fit index cutoff plots

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_wlsmv_dfi_plots.pdf"))
ma_cfa_28_3cf_wlsmv_dfi$plots
dev.off()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_28_3cf_wlsmv_resid_std <- lavResiduals(ma_cfa_28_3cf_wlsmv_fit, type = "raw")$cov.z

ma_cfa_28_3cf_wlsmv_resid_std_sig <- ma_cfa_28_3cf_wlsmv_resid_std
ma_cfa_28_3cf_wlsmv_resid_std_sig[abs(ma_cfa_28_3cf_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_28_3cf_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_28_3cf_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_28_3cf_wlsmv_mi <- modificationIndices(ma_cfa_28_3cf_wlsmv_fit)

ma_cfa_28_3cf_wlsmv_mi_sig <- ma_cfa_28_3cf_wlsmv_mi[ma_cfa_28_3cf_wlsmv_mi$mi > 3.84, ]
ma_cfa_28_3cf_wlsmv_mi_sig <- ma_cfa_28_3cf_wlsmv_mi_sig[order(ma_cfa_28_3cf_wlsmv_mi_sig$mi,
                                                               decreasing = TRUE), ]

write.csv(ma_cfa_28_3cf_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_28_3cf_wlsmv_mi_sig.csv"))

# Inspect residual correlations (focus on at least small ones; i.e., > .1; Flora, 2020, p. 490)

ma_cfa_28_3cf_wlsmv_resid_cor <- round(residuals(ma_cfa_28_3cf_wlsmv_fit, type = "cor")$cov, 2)

ma_cfa_28_3cf_wlsmv_resid_cor_small <- ma_cfa_28_3cf_wlsmv_resid_cor
ma_cfa_28_3cf_wlsmv_resid_cor_small[abs(ma_cfa_28_3cf_wlsmv_resid_cor_small) < .1] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_wlsmv_resid_cor_small_plot.pdf"))
ma_cfa_28_3cf_wlsmv_resid_cor_small_p <- 
  corrplot(ma_cfa_28_3cf_wlsmv_resid_cor_small, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey",
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Plot completely standardized solution

ma_cfa_28_3cf_wlsmv_p <- 
  semPaths(ma_cfa_28_3cf_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "neg_thr", "pos_thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 30, 1, 12),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_3cf_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for revised 28-item model with 2 correlated negative threat/not negative threat factors ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(neg_thr     = dat_items$rr_28_neg_thr_items,
            not_neg_thr = c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items))

# TODO: Rename model to "ma_cfa_28_2cf_neg_thr_not"





ma_cfa_28_2cf <- write_lavaan(latent = lat)
cat(ma_cfa_28_2cf)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_2cf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_2cf_wlsmv_fit <- cfa(ma_cfa_28_2cf, data = dat_ma_rr_only, std.lv = TRUE,
                               ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_2cf_wlsmv_summ <- summary(ma_cfa_28_2cf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_2cf_wlsmv_std_all <- standardizedsolution(ma_cfa_28_2cf_wlsmv_fit))
cat("\n", "DYNAMIC FIT INDEX CUTOFFS:", "\n\n")
(ma_cfa_28_2cf_wlsmv_dfi <- dynamic::catHB(ma_cfa_28_2cf_wlsmv_fit, plot = TRUE, 
                                           estimator = "WLSMV", reps = 250))
sink(type = "message")
sink()

# Save dynamic fit index cutoff plots

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_2cf_wlsmv_dfi_plots.pdf"))
ma_cfa_28_2cf_wlsmv_dfi$plots
dev.off()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_28_2cf_wlsmv_resid_std <- lavResiduals(ma_cfa_28_2cf_wlsmv_fit, type = "raw")$cov.z

ma_cfa_28_2cf_wlsmv_resid_std_sig <- ma_cfa_28_2cf_wlsmv_resid_std
ma_cfa_28_2cf_wlsmv_resid_std_sig[abs(ma_cfa_28_2cf_wlsmv_resid_std_sig) < 1.96] <- NA

# TODO: Revise lines





pdf(file = paste0(cfa_revised_path, "ma_cfa_28_2cf_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_28_2cf_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_28_2cf_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("neg_thr_elevator_1a", "pos_thr_noise_4c", "neg_non_blood_test_9b"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_28_2cf_wlsmv_mi <- modificationIndices(ma_cfa_28_2cf_wlsmv_fit)

ma_cfa_28_2cf_wlsmv_mi_sig <- ma_cfa_28_2cf_wlsmv_mi[ma_cfa_28_2cf_wlsmv_mi$mi > 3.84, ]
ma_cfa_28_2cf_wlsmv_mi_sig <- ma_cfa_28_2cf_wlsmv_mi_sig[order(ma_cfa_28_2cf_wlsmv_mi_sig$mi,
                                                               decreasing = TRUE), ]

write.csv(ma_cfa_28_2cf_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_28_2cf_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_28_2cf_wlsmv_p <- 
  semPaths(ma_cfa_28_2cf_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2.5, curvature = 4.5, cardinal = TRUE, reorder = FALSE,
           latents = c("not_neg_thr", "neg_thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items)),
           sizeLat = 9, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 10, 1, 32),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_2cf_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run EFA in CFA for revised 3 correlated-factors 28-item model ----
# ---------------------------------------------------------------------------- #

# Define model (see Brown, 2015, p. 168). Based on high loadings and low cross-
# loadings relative to other items across rotations in EFA based on 28 items, use 
# "pos_thr_blood_test_9d", "neg_thr_meeting_friend_5d", and "neg_non_meeting_friend_5b" 
# as anchor items for "pos_thr", "neg_thr", and "non" factors, respectively.

rr_28_items_ordered <- c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_neg_thr_items,
                         dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)

lat <- list(pos_thr = setdiff(rr_28_items_ordered,
                              c("neg_thr_meeting_friend_5d", "neg_non_meeting_friend_5b")),
            neg_thr = setdiff(rr_28_items_ordered,
                              c("pos_thr_blood_test_9d", "neg_non_meeting_friend_5b")),
            non     = setdiff(rr_28_items_ordered,
                              c("pos_thr_blood_test_9d", "neg_thr_meeting_friend_5d")))

ma_cfa_28_3cf_efa <- write_lavaan(latent = lat)
cat(ma_cfa_28_3cf_efa)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_3cf_efa_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_3cf_efa_wlsmv_fit <- cfa(ma_cfa_28_3cf_efa, data = dat_ma_rr_only, std.lv = TRUE,
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_3cf_efa_wlsmv_summ <- summary(ma_cfa_28_3cf_efa_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_3cf_efa_wlsmv_std_all <- standardizedsolution(ma_cfa_28_3cf_efa_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_28_3cf_efa_wlsmv_resid_std <- lavResiduals(ma_cfa_28_3cf_efa_wlsmv_fit, type = "raw")$cov.z

ma_cfa_28_3cf_efa_wlsmv_resid_std_sig <- ma_cfa_28_3cf_efa_wlsmv_resid_std
ma_cfa_28_3cf_efa_wlsmv_resid_std_sig[abs(ma_cfa_28_3cf_efa_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_efa_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_28_3cf_efa_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_28_3cf_efa_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_28_3cf_efa_wlsmv_mi <- modificationIndices(ma_cfa_28_3cf_efa_wlsmv_fit)

ma_cfa_28_3cf_efa_wlsmv_mi_sig <- ma_cfa_28_3cf_efa_wlsmv_mi[ma_cfa_28_3cf_efa_wlsmv_mi$mi > 3.84, ]
ma_cfa_28_3cf_efa_wlsmv_mi_sig <- ma_cfa_28_3cf_efa_wlsmv_mi_sig[order(ma_cfa_28_3cf_efa_wlsmv_mi_sig$mi,
                                                                       decreasing = TRUE), ]

write.csv(ma_cfa_28_3cf_efa_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_28_3cf_efa_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_28_3cf_efa_wlsmv_p <- 
  semPaths(ma_cfa_28_3cf_efa_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "neg_thr", "pos_thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 10, 1, 5),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_3cf_efa_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for revised 3 correlated-factors 28-item model with correlated errors per scenario for threat items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_28_pos_thr_items,
            neg_thr = dat_items$rr_28_neg_thr_items,
            non     = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items))

cov <- list(pos_thr_noise_4c      = "neg_thr_noise_4b",
            pos_thr_lunch_6c      = "neg_thr_lunch_6a",
            pos_thr_shopping_8b   = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d = "neg_thr_blood_test_9c")

ma_cfa_28_3cf_ce_thr <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_3cf_ce_thr)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_thr_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_3cf_ce_thr_wlsmv_fit <- cfa(ma_cfa_28_3cf_ce_thr, data = dat_ma_rr_only, std.lv = TRUE,
                                      ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_3cf_ce_thr_wlsmv_summ <- summary(ma_cfa_28_3cf_ce_thr_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_3cf_ce_thr_wlsmv_std_all <- standardizedsolution(ma_cfa_28_3cf_ce_thr_wlsmv_fit))
cat("\n", "DYNAMIC FIT INDEX CUTOFFS:", "\n\n")
(ma_cfa_28_3cf_ce_thr_wlsmv_dfi <- dynamic::catHB(ma_cfa_28_3cf_ce_thr_wlsmv_fit, plot = TRUE, 
                                                  estimator = "WLSMV", reps = 250))
sink(type = "message")
sink()

# Save dynamic fit index cutoff plots

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_thr_wlsmv_dfi_plots.pdf"))
ma_cfa_28_3cf_ce_thr_wlsmv_dfi$plots
dev.off()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_28_3cf_ce_thr_wlsmv_resid_std <- lavResiduals(ma_cfa_28_3cf_ce_thr_wlsmv_fit, type = "raw")$cov.z

ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig <- ma_cfa_28_3cf_ce_thr_wlsmv_resid_std
ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig[abs(ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_28_3cf_ce_thr_wlsmv_mi <- modificationIndices(ma_cfa_28_3cf_ce_thr_wlsmv_fit)

ma_cfa_28_3cf_ce_thr_wlsmv_mi_sig <- ma_cfa_28_3cf_ce_thr_wlsmv_mi[ma_cfa_28_3cf_ce_thr_wlsmv_mi$mi > 3.84, ]
ma_cfa_28_3cf_ce_thr_wlsmv_mi_sig <- ma_cfa_28_3cf_ce_thr_wlsmv_mi_sig[order(ma_cfa_28_3cf_ce_thr_wlsmv_mi_sig$mi,
                                                                             decreasing = TRUE), ]

write.csv(ma_cfa_28_3cf_ce_thr_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_thr_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_28_3cf_ce_thr_wlsmv_p <- 
  semPaths(ma_cfa_28_3cf_ce_thr_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "neg_thr", "pos_thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 20, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_3cf_ce_thr_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for revised 3 correlated-factors 28-item model with correlated errors per scenario for all items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_28_pos_thr_items,
            neg_thr = dat_items$rr_28_neg_thr_items,
            non     = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items))

cov <- list(neg_thr_elevator_1a       = c("pos_non_elevator_1b", "neg_non_elevator_1c"),
            pos_non_elevator_1b       = "neg_non_elevator_1c",
            neg_thr_wedding_2c        = c("pos_non_wedding_2b", "neg_non_wedding_2d"),
            pos_non_wedding_2b        = "neg_non_wedding_2d",
            neg_thr_job_3c            = c("pos_non_job_3b", "neg_non_job_3d"),
            pos_non_job_3b            = "neg_non_job_3d",
            pos_thr_noise_4c          = c("neg_thr_noise_4b", "pos_non_noise_4a", "neg_non_noise_4d"),
            neg_thr_noise_4b          = c("pos_non_noise_4a", "neg_non_noise_4d"),
            pos_non_noise_4a          = "neg_non_noise_4d",
            neg_thr_meeting_friend_5d = c("pos_non_meeting_friend_5a", "neg_non_meeting_friend_5b"),
            pos_non_meeting_friend_5a = "neg_non_meeting_friend_5b",
            pos_thr_lunch_6c          = c("neg_thr_lunch_6a", "neg_non_lunch_6b"),
            neg_thr_lunch_6a          = "neg_non_lunch_6b",
            pos_thr_scrape_7c         = "neg_non_scrape_7b",
            pos_thr_shopping_8b       = c("neg_thr_shopping_8a", "pos_non_shopping_8c", "neg_non_shopping_8d"),
            neg_thr_shopping_8a       = c("pos_non_shopping_8c", "neg_non_shopping_8d"),
            pos_non_shopping_8c       = "neg_non_shopping_8d",
            pos_thr_blood_test_9d     = c("neg_thr_blood_test_9c", "neg_non_blood_test_9b"),
            neg_thr_blood_test_9c     = "neg_non_blood_test_9b")

ma_cfa_28_3cf_ce_all <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_3cf_ce_all)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_all_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_3cf_ce_all_wlsmv_fit <- cfa(ma_cfa_28_3cf_ce_all, data = dat_ma_rr_only, std.lv = TRUE,
                                      ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_3cf_ce_all_wlsmv_summ <- summary(ma_cfa_28_3cf_ce_all_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_3cf_ce_all_wlsmv_std_all <- standardizedsolution(ma_cfa_28_3cf_ce_all_wlsmv_fit))
cat("\n", "DYNAMIC FIT INDEX CUTOFFS:", "\n\n")
(ma_cfa_28_3cf_ce_all_wlsmv_dfi <- dynamic::catHB(ma_cfa_28_3cf_ce_all_wlsmv_fit, plot = TRUE, 
                                                  estimator = "WLSMV", reps = 250))
sink(type = "message")
sink()

# Save dynamic fit index cutoff plots

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_all_wlsmv_dfi_plots.pdf"))
ma_cfa_28_3cf_ce_all_wlsmv_dfi$plots
dev.off()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_28_3cf_ce_all_wlsmv_resid_std <- lavResiduals(ma_cfa_28_3cf_ce_all_wlsmv_fit, type = "raw")$cov.z

ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig <- ma_cfa_28_3cf_ce_all_wlsmv_resid_std
ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig[abs(ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_28_3cf_ce_all_wlsmv_mi <- modificationIndices(ma_cfa_28_3cf_ce_all_wlsmv_fit)

ma_cfa_28_3cf_ce_all_wlsmv_mi_sig <- ma_cfa_28_3cf_ce_all_wlsmv_mi[ma_cfa_28_3cf_ce_all_wlsmv_mi$mi > 3.84, ]
ma_cfa_28_3cf_ce_all_wlsmv_mi_sig <- ma_cfa_28_3cf_ce_all_wlsmv_mi_sig[order(ma_cfa_28_3cf_ce_all_wlsmv_mi_sig$mi,
                                                                             decreasing = TRUE), ]

write.csv(ma_cfa_28_3cf_ce_all_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_28_3cf_ce_all_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_28_3cf_ce_all_wlsmv_p <- 
  semPaths(ma_cfa_28_3cf_ce_all_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2.5, curvature = 4.5, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "neg_thr", "pos_thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 10, 1, 32),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_3cf_ce_all_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for revised 28-item model with 2 correlated threat/nonthreat factors and 4 correlated positive/negative method factors ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr     = c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_neg_thr_items),
            non     = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items),
            pos_thr = dat_items$rr_28_pos_thr_items,
            neg_thr = dat_items$rr_28_neg_thr_items,
            pos_non = dat_items$rr_28_pos_non_items,
            neg_non = dat_items$rr_28_neg_non_items)

cov <- list(thr     = "non",
            pos_thr = c("neg_thr", "pos_non", "neg_non"),
            neg_thr = c("pos_non", "neg_non"),
            pos_non = "neg_non")

# TODO: Rename model to "ma_cfa_28_2cf_thr_non_4cmf"





ma_cfa_28_2cf_4cmf <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_2cf_4cmf)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_2cf_4cmf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_2cf_4cmf_wlsmv_fit <- cfa(ma_cfa_28_2cf_4cmf, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                    ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_2cf_4cmf_wlsmv_summ <- summary(ma_cfa_28_2cf_4cmf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_2cf_4cmf_wlsmv_std_all <- standardizedsolution(ma_cfa_28_2cf_4cmf_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve issues below (but note that Brown, 2015, book, p. 199, states that 
# such a correlated methods model is usually empirically underidentified)

# Heywood case: Standardized factor correlation between "thr" and "non" > 1

# Warning message:
#   In lav_object_post_check(object) :
#   lavaan WARNING: covariance matrix of latent variables
# is not positive definite;
# use lavInspect(fit, "cov.lv") to investigate.

raw_cov                               <- cov(dat_ma_rr_only[dat_items$rr_28_items], use = "pairwise.complete.obs")
ma_cfa_28_2cf_4cmf_wlsmv_fit_samp_cov <- lavInspect(ma_cfa_28_2cf_4cmf_wlsmv_fit, "sampstat")$cov
ma_cfa_28_2cf_4cmf_wlsmv_fit_cov_ov   <- lavInspect(ma_cfa_28_2cf_4cmf_wlsmv_fit, "cov.ov")
ma_cfa_28_2cf_4cmf_wlsmv_fit_cov_lv   <- lavInspect(ma_cfa_28_2cf_4cmf_wlsmv_fit, "cov.lv")

eigen(raw_cov)$values
eigen(ma_cfa_28_2cf_4cmf_wlsmv_fit_samp_cov)$values
eigen(ma_cfa_28_2cf_4cmf_wlsmv_fit_cov_ov)$values
eigen(ma_cfa_28_2cf_4cmf_wlsmv_fit_cov_lv)$values # TODO: Negative eigenvalue

ma_cfa_28_2cf_4cmf_wlsmv_fit_information <- lavInspect(ma_cfa_28_2cf_4cmf_wlsmv_fit, "information")
eigen(ma_cfa_28_2cf_4cmf_wlsmv_fit_information)$values
det(ma_cfa_28_2cf_4cmf_wlsmv_fit_information)   # Determinant is not 0 (i.e., matrix not singular)
# solve(ma_cfa_28_2cf_4cmf_wlsmv_fit_information)

# Create plot without parameter estimates just to show structure

ma_cfa_28_2cf_4cmf_wlsmv_p <- 
  semPaths(ma_cfa_28_2cf_4cmf_wlsmv_fit, what = "path", whatLabels = "no", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = c("non", "thr"),
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr", "non", "thr"),
           manifests = c(rev(dat_items$rr_28_neg_items), rev(dat_items$rr_28_pos_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 10),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_2cf_4cmf_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_28_2cf_4cmf_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:28, 2]  <- 1.5*pi
ma_cfa_28_2cf_4cmf_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[29:56, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_28_2cf_4cmf_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_28_2cf_4cmf_wlsmv_p$Arguments$edge.color == 
                                                      "transparent"] <- "transparent"

# Make edges black

ma_cfa_28_2cf_4cmf_wlsmv_p$graphAttributes$Edges$color[ma_cfa_28_2cf_4cmf_wlsmv_p$graphAttributes$Edges$color == 
                                                "#808080FF"] <- "black"

plot(ma_cfa_28_2cf_4cmf_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for revised 28-item model with 2 correlated threat/nonthreat factors and 2 correlated positive/negative method factors for threat items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr     = c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_neg_thr_items),
            non     = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items),
            pos_thr = dat_items$rr_28_pos_thr_items,
            neg_thr = dat_items$rr_28_neg_thr_items)

cov <- list(thr     = "non",
            pos_thr = "neg_thr")

# TODO: Rename model to "ma_cfa_28_2cf_thr_non_2cmf_thr"





ma_cfa_28_2cf_2cmf_thr <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_2cf_2cmf_thr)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_2cf_2cmf_thr_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_2cf_2cmf_thr_wlsmv_fit <- cfa(ma_cfa_28_2cf_2cmf_thr, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                        ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_2cf_2cmf_thr_wlsmv_summ <- summary(ma_cfa_28_2cf_2cmf_thr_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_2cf_2cmf_thr_wlsmv_std_all <- standardizedsolution(ma_cfa_28_2cf_2cmf_thr_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below

# Warning message:
# In lavaan::lavaan(model = ma_cfa_28_2cf_2cmf_thr, data = dat_ma_rr_only,  :
#   lavaan WARNING:
#     the optimizer (NLMINB) claimed the model converged, but not all
#     elements of the gradient are (near) zero; the optimizer may not
#     have found a local solution use check.gradient = FALSE to skip
#     this check.
# SUMMARY: 
# 
# lavaan 0.6.17 did NOT end normally after 1086 iterations
# ** WARNING ** Estimates below are most likely unreliable

# Warning message:
# In lav_object_summary(object = object, header = header, fit.measures = fit.measures,  :
#   lavaan WARNING: fit measures not available if model did not converge

# Warning message:
# In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = object@SampleStats,  :
#   lavaan WARNING:
#     Could not compute standard errors! The information matrix could
#     not be inverted. This may be a symptom that the model is not
#     identified.

# Create plot without parameter estimates just to show structure

ma_cfa_28_2cf_2cmf_thr_wlsmv_p <- 
  semPaths(ma_cfa_28_2cf_2cmf_thr_wlsmv_fit, what = "path", whatLabels = "no", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = c("non", "thr"),
           latents = c("neg_thr", "pos_thr", "non", "thr"),
           manifests = c(rev(dat_items$rr_28_neg_items), rev(dat_items$rr_28_pos_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 10),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_2cf_2cmf_thr_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_28_2cf_2cmf_thr_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:28, 2]  <- 1.5*pi
ma_cfa_28_2cf_2cmf_thr_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[29:41, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_28_2cf_2cmf_thr_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_28_2cf_2cmf_thr_wlsmv_p$Arguments$edge.color == 
                                                               "transparent"] <- "transparent"

# Make edges black

ma_cfa_28_2cf_2cmf_thr_wlsmv_p$graphAttributes$Edges$color[ma_cfa_28_2cf_2cmf_thr_wlsmv_p$graphAttributes$Edges$color == 
                                                         "#808080FF"] <- "black"

plot(ma_cfa_28_2cf_2cmf_thr_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for revised 28-item model with 2 correlated threat/nonthreat factors and 2 uncorrelated positive/negative method factors for threat items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr     = c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_neg_thr_items),
            non     = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items),
            pos_thr = dat_items$rr_28_pos_thr_items,
            neg_thr = dat_items$rr_28_neg_thr_items)

cov <- list(thr     = "non")

# TODO: Rename model to "ma_cfa_28_2cf_thr_non_2mf_thr"





ma_cfa_28_2cf_2mf_thr <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_2cf_2mf_thr)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_2cf_2mf_thr_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_2cf_2mf_thr_wlsmv_fit <- cfa(ma_cfa_28_2cf_2mf_thr, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                       ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_2cf_2mf_thr_wlsmv_summ <- summary(ma_cfa_28_2cf_2mf_thr_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_2cf_2mf_thr_wlsmv_std_all <- standardizedsolution(ma_cfa_28_2cf_2mf_thr_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below

# Warning message:
# In lavaan::lavaan(model = ma_cfa_28_2cf_2mf_thr, data = dat_ma_rr_only,  :
#   lavaan WARNING:
#     the optimizer (NLMINB) claimed the model converged, but not all
#     elements of the gradient are (near) zero; the optimizer may not
#     have found a local solution use check.gradient = FALSE to skip
#     this check.
# SUMMARY: 
# 
# lavaan 0.6.17 did NOT end normally after 1136 iterations
# ** WARNING ** Estimates below are most likely unreliable

# Warning message:
# In lav_object_summary(object = object, header = header, fit.measures = fit.measures,  :
#   lavaan WARNING: fit measures not available if model did not converge

# Warning message:
# In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = object@SampleStats,  :
#   lavaan WARNING:
#     Could not compute standard errors! The information matrix could
#     not be inverted. This may be a symptom that the model is not
#     identified.

# Create plot without parameter estimates just to show structure

ma_cfa_28_2cf_2mf_thr_wlsmv_p <- 
  semPaths(ma_cfa_28_2cf_2mf_thr_wlsmv_fit, what = "path", whatLabels = "no", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = c("non", "thr"),
           latents = c("neg_thr", "pos_thr", "non", "thr"),
           manifests = c(rev(dat_items$rr_28_neg_items), rev(dat_items$rr_28_pos_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 10),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_2cf_2mf_thr_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_28_2cf_2mf_thr_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:28, 2]  <- 1.5*pi
ma_cfa_28_2cf_2mf_thr_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[29:41, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_28_2cf_2mf_thr_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_28_2cf_2mf_thr_wlsmv_p$Arguments$edge.color == 
                                                                   "transparent"] <- "transparent"

# Make edges black

ma_cfa_28_2cf_2mf_thr_wlsmv_p$graphAttributes$Edges$color[ma_cfa_28_2cf_2mf_thr_wlsmv_p$graphAttributes$Edges$color == 
                                                             "#808080FF"] <- "black"

plot(ma_cfa_28_2cf_2mf_thr_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for revised 28-item model with 2 correlated threat/nonthreat factors and correlated errors per scenario and valence for threat items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr = c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_neg_thr_items),
            non = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items))

cov <- list(pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c",
            
            pos_thr_noise_4c          = c("pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_lunch_6c          = c("pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_scrape_7c         = c("pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_shopping_8b       = "pos_thr_blood_test_9d",
            
            neg_thr_elevator_1a       = c("neg_thr_wedding_2c", "neg_thr_job_3c", 
                                          "neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_wedding_2c        = c("neg_thr_job_3c", "neg_thr_noise_4b", 
                                          "neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_job_3c            = c("neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_noise_4b          = c("neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_meeting_friend_5d = c("neg_thr_lunch_6a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_lunch_6a          = c("neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_shopping_8a       = "neg_thr_blood_test_9c")

ma_cfa_28_2cf_thr_non_ce_val_scen_thr <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_2cf_thr_non_ce_val_scen_thr)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit <- 
  cfa(ma_cfa_28_2cf_thr_non_ce_val_scen_thr, data = dat_ma_rr_only, std.lv = TRUE,
      ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_summ <- 
  summary(ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_std_all <- 
  standardizedsolution(ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below (also, e.g., residual variances of "thr" items are 1
# and standardized factor correlation between "thr" and "non" is > 1). Note: Including 
# such correlated errors for both valences (vs. one valence) tends to yield improper 
# solutions (see Marsh, 1996; Brown, 2003)

# Warning message:
# In lavaan::lavaan(model = ma_cfa_28_2cf_thr_non_ce_val_scen_thr,  :
#   lavaan WARNING:
#     the optimizer warns that a solution has NOT been found!
# SUMMARY: 
# 
# lavaan 0.6.17 did NOT end normally after 1050 iterations
# ** WARNING ** Estimates below are most likely unreliable

# Warning message:
# In lav_object_summary(object = object, header = header, fit.measures = fit.measures,  :
#   lavaan WARNING: fit measures not available if model did not converge

# Warning message:
# In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = object@SampleStats,  :
#   lavaan WARNING:
#     Could not compute standard errors! The information matrix could
#     not be inverted. This may be a symptom that the model is not
#     identified.

# Create plot without parameter estimates just to show structure

ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_p <- 
  semPaths(ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit, what = "path", whatLabels = "no", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 12, 1, 32),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for revised 28-item model with 2 correlated threat/nonthreat factors, correlated errors per scenario for threat items, and correlated errors among positive threat items ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr = c(dat_items$rr_28_pos_thr_items, dat_items$rr_28_neg_thr_items),
            non = c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items))

cov <- list(pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c",
            
            pos_thr_noise_4c          = c("pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_lunch_6c          = c("pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_scrape_7c         = c("pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_shopping_8b       = "pos_thr_blood_test_9d")

ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit <- 
  cfa(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr, data = dat_ma_rr_only, std.lv = TRUE,
      ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_summ <- 
    summary(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_std_all <- 
    standardizedsolution(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std <- 
  lavResiduals(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit, type = "raw")$cov.z

ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig <- ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std
ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig[abs(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_noise_4c", "neg_thr_elevator_1a", 
                               "pos_non_elevator_1b", "neg_non_blood_test_9b"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi <- modificationIndices(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit)

ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi_sig <- 
  ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi[ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi$mi > 3.84, ]
ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi_sig <- 
  ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi_sig[order(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi_sig$mi,
                                                                             decreasing = TRUE), ]

write.csv(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_p <- 
  semPaths(ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = c("non", "thr"),
           manifests = c(rev(c(dat_items$rr_28_pos_non_items, dat_items$rr_28_neg_non_items)), 
                         rev(dat_items$rr_28_neg_thr_items), rev(dat_items$rr_28_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 12, 1, 32),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 2 correlated-factors all threat-items model ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items)

ma_cfa_thr_2cf <- write_lavaan(latent = lat)
cat(ma_cfa_thr_2cf)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_2cf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_2cf_wlsmv_fit <- cfa(ma_cfa_thr_2cf, data = dat_ma_rr_only, std.lv = TRUE,
                                ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_2cf_wlsmv_summ <- summary(ma_cfa_thr_2cf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_2cf_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_2cf_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_2cf_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_2cf_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_2cf_wlsmv_resid_std_sig <- ma_cfa_thr_2cf_wlsmv_resid_std
ma_cfa_thr_2cf_wlsmv_resid_std_sig[abs(ma_cfa_thr_2cf_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_2cf_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_2cf_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_2cf_wlsmv_mi <- modificationIndices(ma_cfa_thr_2cf_wlsmv_fit)

ma_cfa_thr_2cf_wlsmv_mi_sig <- ma_cfa_thr_2cf_wlsmv_mi[ma_cfa_thr_2cf_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_2cf_wlsmv_mi_sig <- ma_cfa_thr_2cf_wlsmv_mi_sig[order(ma_cfa_thr_2cf_wlsmv_mi_sig$mi,
                                                                 decreasing = TRUE), ]

write.csv(ma_cfa_thr_2cf_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_2cf_wlsmv_p <- 
  semPaths(ma_cfa_thr_2cf_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_thr", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 20, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_2cf_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run EFA in CFA for 2 correlated-factors all threat-items model ----
# ---------------------------------------------------------------------------- #

# Define model (see Brown, 2015, p. 168). Based on high loadings and low cross-
# loadings relative to other items across rotations in EFA based on threat items, 
# use "pos_thr_blood_test_9d" and "neg_thr_meeting_friend_5d" as anchor items for 
# "pos_thr" and "neg_thr" factors, respectively.

rr_thr_items_ordered <- c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items)

lat <- list(pos_thr = setdiff(rr_thr_items_ordered, "neg_thr_meeting_friend_5d"),
            neg_thr = setdiff(rr_thr_items_ordered, "pos_thr_blood_test_9d"))

ma_cfa_thr_2cf_efa <- write_lavaan(latent = lat)
cat(ma_cfa_thr_2cf_efa)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_2cf_efa_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_2cf_efa_wlsmv_fit <- cfa(ma_cfa_thr_2cf_efa, data = dat_ma_rr_only, std.lv = TRUE,
                                    ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_2cf_efa_wlsmv_summ <- summary(ma_cfa_thr_2cf_efa_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_2cf_efa_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_2cf_efa_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value of 
# standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_2cf_efa_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_2cf_efa_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig <- ma_cfa_thr_2cf_efa_wlsmv_resid_std
ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig[abs(ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_meeting_friend_5d"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_2cf_efa_wlsmv_mi <- modificationIndices(ma_cfa_thr_2cf_efa_wlsmv_fit)

ma_cfa_thr_2cf_efa_wlsmv_mi_sig <- ma_cfa_thr_2cf_efa_wlsmv_mi[ma_cfa_thr_2cf_efa_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_2cf_efa_wlsmv_mi_sig <- ma_cfa_thr_2cf_efa_wlsmv_mi_sig[order(ma_cfa_thr_2cf_efa_wlsmv_mi_sig$mi,
                                                                         decreasing = TRUE), ]

write.csv(ma_cfa_thr_2cf_efa_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_efa_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_2cf_efa_wlsmv_p <- 
  semPaths(ma_cfa_thr_2cf_efa_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_thr", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 10, 1, 5),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_2cf_efa_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 2 correlated-factors all threat-items model with correlated errors per scenario ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items)

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_2cf_ce <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_2cf_ce)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_2cf_ce_wlsmv_fit <- cfa(ma_cfa_thr_2cf_ce, data = dat_ma_rr_only, std.lv = TRUE,
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_2cf_ce_wlsmv_summ <- summary(ma_cfa_thr_2cf_ce_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_2cf_ce_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_2cf_ce_wlsmv_fit))
cat("\n", "DYNAMIC FIT INDEX CUTOFFS:", "\n\n")
(ma_cfa_thr_2cf_ce_wlsmv_dfi <- dynamic::catHB(ma_cfa_thr_2cf_ce_wlsmv_fit, plot = TRUE,
                                               estimator = "WLSMV", reps = 250))
sink(type = "message")
sink()

# Save dynamic fit index cutoff plots

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_wlsmv_dfi_plots.pdf"))
ma_cfa_thr_2cf_ce_wlsmv_dfi$plots
dev.off()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_2cf_ce_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_2cf_ce_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig <- ma_cfa_thr_2cf_ce_wlsmv_resid_std
ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig[abs(ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_2cf_ce_wlsmv_mi <- modificationIndices(ma_cfa_thr_2cf_ce_wlsmv_fit)

ma_cfa_thr_2cf_ce_wlsmv_mi_sig <- ma_cfa_thr_2cf_ce_wlsmv_mi[ma_cfa_thr_2cf_ce_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_2cf_ce_wlsmv_mi_sig <- ma_cfa_thr_2cf_ce_wlsmv_mi_sig[order(ma_cfa_thr_2cf_ce_wlsmv_mi_sig$mi,
                                                                       decreasing = TRUE), ]

write.csv(ma_cfa_thr_2cf_ce_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_2cf_ce_wlsmv_p <- 
  semPaths(ma_cfa_thr_2cf_ce_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_thr", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 20, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_2cf_ce_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items))

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_ce_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_ce_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_ce_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_ce_scen, data = dat_ma_rr_only, std.lv = TRUE,
                                       ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_ce_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_ce_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_ce_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_ce_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_ce_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_ce_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig <- ma_cfa_thr_1f_ce_scen_wlsmv_resid_std
ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_ce_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_ce_scen_wlsmv_fit)

ma_cfa_thr_1f_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_scen_wlsmv_mi[ma_cfa_thr_1f_ce_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_scen_wlsmv_mi_sig[order(ma_cfa_thr_1f_ce_scen_wlsmv_mi_sig$mi,
                                                                               decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_ce_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_ce_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_ce_scen_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = "thr",
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_ce_scen_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario and per valence ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items))

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c",
            
            pos_thr_elevator_1d       = c("pos_thr_wedding_2a", "pos_thr_job_3a", 
                                          "pos_thr_noise_4c", "pos_thr_meeting_friend_5c", 
                                          "pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_wedding_2a        = c("pos_thr_job_3a", "pos_thr_noise_4c", 
                                          "pos_thr_meeting_friend_5c", "pos_thr_lunch_6c", 
                                          "pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_job_3a            = c("pos_thr_noise_4c", "pos_thr_meeting_friend_5c", 
                                          "pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_noise_4c          = c("pos_thr_meeting_friend_5c", "pos_thr_lunch_6c", 
                                          "pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_meeting_friend_5c = c("pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_lunch_6c          = c("pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_scrape_7c         = c("pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_shopping_8b       = "pos_thr_blood_test_9d",
            
            neg_thr_elevator_1a       = c("neg_thr_wedding_2c", "neg_thr_job_3c", 
                                          "neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_wedding_2c        = c("neg_thr_job_3c", "neg_thr_noise_4b", 
                                          "neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_job_3c            = c("neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_noise_4b          = c("neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_meeting_friend_5d = c("neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_lunch_6a          = c("neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_scrape_7a         = c("neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_shopping_8a       = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_ce_val_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_ce_val_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_val_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_ce_val_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_ce_val_scen, data = dat_ma_rr_only, std.lv = TRUE,
                                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_ce_val_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_ce_val_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_ce_val_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_ce_val_scen_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below

# Warning messages:
# 1: In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = lavsamplestats,  :
#   lavaan WARNING:
#     Could not compute standard errors! The information matrix could
#     not be inverted. This may be a symptom that the model is not
#     identified.
# 2: In lav_test_satorra_bentler(lavobject = NULL, lavsamplestats = lavsamplestats,  :
#   lavaan WARNING: could not invert information matrix needed for robust test statistic

# Create plot without parameter estimates just to show structure

ma_cfa_thr_1f_ce_val_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_ce_val_scen_wlsmv_fit, what = "path", whatLabels = "no", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = "thr",
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_ce_val_scen_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario and for positive valence ----
# ---------------------------------------------------------------------------- #

# Define model (e.g., see Brown, 2003, "one factor with method effects", p. 1417)

lat <- list(thr = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items))

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c",
            
            pos_thr_elevator_1d       = c("pos_thr_wedding_2a", "pos_thr_job_3a", 
                                          "pos_thr_noise_4c", "pos_thr_meeting_friend_5c", 
                                          "pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_wedding_2a        = c("pos_thr_job_3a", "pos_thr_noise_4c", 
                                          "pos_thr_meeting_friend_5c", "pos_thr_lunch_6c", 
                                          "pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_job_3a            = c("pos_thr_noise_4c", "pos_thr_meeting_friend_5c", 
                                          "pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_noise_4c          = c("pos_thr_meeting_friend_5c", "pos_thr_lunch_6c", 
                                          "pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_meeting_friend_5c = c("pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_lunch_6c          = c("pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_scrape_7c         = c("pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_shopping_8b       = "pos_thr_blood_test_9d")

ma_cfa_thr_1f_ce_pos_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_ce_pos_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_pos_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_ce_pos_scen, data = dat_ma_rr_only, std.lv = TRUE,
                                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_ce_pos_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_ce_pos_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig <- ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std
ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit)

ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi[ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi_sig[order(ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi_sig$mi,
                                                                                       decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_ce_pos_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = "thr",
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_ce_pos_scen_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario and for negative valence ----
# ---------------------------------------------------------------------------- #

# Define model (e.g., see Brown, 2003, "one factor with method effects", p. 1417)

lat <- list(thr = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items))

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c",
            
            neg_thr_elevator_1a       = c("neg_thr_wedding_2c", "neg_thr_job_3c", 
                                          "neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_wedding_2c        = c("neg_thr_job_3c", "neg_thr_noise_4b", 
                                          "neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_job_3c            = c("neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_noise_4b          = c("neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_meeting_friend_5d = c("neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_lunch_6a          = c("neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_scrape_7a         = c("neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_shopping_8a       = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_ce_neg_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_ce_neg_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_neg_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_ce_neg_scen, data = dat_ma_rr_only, std.lv = TRUE,
                                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_ce_neg_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_ce_neg_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig <- ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std
ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit)

ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi[ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi_sig[order(ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi_sig$mi,
                                                                                       decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_ce_neg_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = "thr",
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_ce_neg_scen_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per valence ----
# ---------------------------------------------------------------------------- #

# Define model (e.g., see Bachman & O'Malley, 1986, p. 40; although using correlated
# errors for both valences typically yield improper solutions; Brown, 2003, p. 1421;
# Marsh, 1996, pp. 813 and 815)

lat <- list(thr = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items))

cov <- list(pos_thr_elevator_1d       = c("pos_thr_wedding_2a", "pos_thr_job_3a", 
                                          "pos_thr_noise_4c", "pos_thr_meeting_friend_5c", 
                                          "pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_wedding_2a        = c("pos_thr_job_3a", "pos_thr_noise_4c", 
                                          "pos_thr_meeting_friend_5c", "pos_thr_lunch_6c", 
                                          "pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_job_3a            = c("pos_thr_noise_4c", "pos_thr_meeting_friend_5c", 
                                          "pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_noise_4c          = c("pos_thr_meeting_friend_5c", "pos_thr_lunch_6c", 
                                          "pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_meeting_friend_5c = c("pos_thr_lunch_6c", "pos_thr_scrape_7c", 
                                          "pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_lunch_6c          = c("pos_thr_scrape_7c", "pos_thr_shopping_8b", 
                                          "pos_thr_blood_test_9d"),
            pos_thr_scrape_7c         = c("pos_thr_shopping_8b", "pos_thr_blood_test_9d"),
            pos_thr_shopping_8b       = "pos_thr_blood_test_9d",
            
            neg_thr_elevator_1a       = c("neg_thr_wedding_2c", "neg_thr_job_3c", 
                                          "neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_wedding_2c        = c("neg_thr_job_3c", "neg_thr_noise_4b", 
                                          "neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_job_3c            = c("neg_thr_noise_4b", "neg_thr_meeting_friend_5d", 
                                          "neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_noise_4b          = c("neg_thr_meeting_friend_5d", "neg_thr_lunch_6a", 
                                          "neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_meeting_friend_5d = c("neg_thr_lunch_6a", "neg_thr_scrape_7a", 
                                          "neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_lunch_6a          = c("neg_thr_scrape_7a", "neg_thr_shopping_8a", 
                                          "neg_thr_blood_test_9c"),
            neg_thr_scrape_7a         = c("neg_thr_shopping_8a", "neg_thr_blood_test_9c"),
            neg_thr_shopping_8a       = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_ce_val <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_ce_val)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_val_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_ce_val_wlsmv_fit <- cfa(ma_cfa_thr_1f_ce_val, data = dat_ma_rr_only, std.lv = TRUE,
                                      ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_ce_val_wlsmv_summ <- summary(ma_cfa_thr_1f_ce_val_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_ce_val_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_ce_val_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below

# Warning messages:
# 1: In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = lavsamplestats,  :
#   lavaan WARNING:
#     Could not compute standard errors! The information matrix could
#     not be inverted. This may be a symptom that the model is not
#     identified.
# 2: In lav_test_satorra_bentler(lavobject = NULL, lavsamplestats = lavsamplestats,  :
#   lavaan WARNING: could not invert information matrix needed for robust test statistic

# Create plot without parameter estimates just to show structure

ma_cfa_thr_1f_ce_val_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_ce_val_wlsmv_fit, what = "path", whatLabels = "no", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = "thr",
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_ce_val_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated positive/negative method factors ----
# ---------------------------------------------------------------------------- #

# Define model (single-trait multimethod model; see Morin et al., 2020; e.g., see
# Rodebaugh et al., 2006, Rodebaugh et al., 2007, p. 196, and Hazlett-Stevens et al., 
# 2004--although it's unclear they correlated the method factors, Morin et al. says 
# they should be allowed to correlate)

lat <- list(thr     = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items),
            pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items)

cov <- list(pos_thr = "neg_thr")

ma_cfa_thr_1f_2cmf <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_2cmf)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_2cmf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_2cmf_wlsmv_fit <- cfa(ma_cfa_thr_1f_2cmf, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                    ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_2cmf_wlsmv_summ <- summary(ma_cfa_thr_1f_2cmf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_2cmf_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_2cmf_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_2cmf_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_2cmf_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig <- ma_cfa_thr_1f_2cmf_wlsmv_resid_std
ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_2cmf_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_2cmf_wlsmv_fit)

ma_cfa_thr_1f_2cmf_wlsmv_mi_sig <- ma_cfa_thr_1f_2cmf_wlsmv_mi[ma_cfa_thr_1f_2cmf_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_2cmf_wlsmv_mi_sig <- ma_cfa_thr_1f_2cmf_wlsmv_mi_sig[order(ma_cfa_thr_1f_2cmf_wlsmv_mi_sig$mi,
                                                                   decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_2cmf_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_2cmf_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_2cmf_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_2cmf_wlsmv_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = "thr",
           latents = c("neg_thr", "pos_thr", "thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 3, 1, 3),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_2cmf_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_thr_1f_2cmf_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:18, 2]  <- 1.5*pi
ma_cfa_thr_1f_2cmf_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[19:36, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_thr_1f_2cmf_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_thr_1f_2cmf_wlsmv_p$Arguments$edge.color == 
                                                               "transparent"] <- "transparent"

# Make edges black

ma_cfa_thr_1f_2cmf_wlsmv_p$graphAttributes$Edges$color[ma_cfa_thr_1f_2cmf_wlsmv_p$graphAttributes$Edges$color == 
                                                         "#808080FF"] <- "black"

plot(ma_cfa_thr_1f_2cmf_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated positive/negative method factors and correlated errors per scenario ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr     = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items),
            pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items)

cov <- list(pos_thr = "neg_thr",
            
            pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_2cmf_ce_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_2cmf_ce_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_2cmf_ce_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_2cmf_ce_scen, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                            ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig <- ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std
ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit)

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi[ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi_sig[order(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi_sig$mi,
                                                                                   decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 4, curvature = 2.5, reorder = FALSE,
           bifactor = "thr",
           latents = c("neg_thr", "pos_thr", "thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 5, 1, 5),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:18, 2]  <- 1.5*pi
ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[19:36, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p$Arguments$edge.color == 
                                                                       "transparent"] <- "transparent"

# Make edges black

ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p$graphAttributes$Edges$color[ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p$graphAttributes$Edges$color == 
                                                                 "#808080FF"] <- "black"

plot(ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario and positive method factor ----
# ---------------------------------------------------------------------------- #

# Define model (equivalent to single-trait correlated methods - 1 model, see Morin et al., 
# 2020, p. 1056--but also a bifactor -1 model in this case; see also Marsh, 1996, Models 4-5, 
# p. 814; Rodebaugh et al., 2006; Rodebaugh et al., 2007, p. 196; and Hazlett-Stevens et al., 
# 2004, "negative method factor model", p. 364)

lat <- list(thr     = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items),
            pos_thr = dat_items$rr_pos_thr_items)

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_1mf_pos_ce_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_1mf_pos_ce_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_1mf_pos_ce_scen, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                               ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig <- ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std
ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit)

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi[ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi_sig[order(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi_sig$mi,
                                                                                               decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 4, curvature = 2.5, reorder = FALSE,
           bifactor = "thr",
           latents = c("pos_thr", "thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 10, 1, 5),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:18, 2]  <- 1.5*pi
ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[19:36, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p$Arguments$edge.color == 
                                                                        "transparent"] <- "transparent"

# Make edges black

ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p$graphAttributes$Edges$color[ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p$graphAttributes$Edges$color == 
                                                                  "#808080FF"] <- "black"

plot(ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario and negative method factor ----
# ---------------------------------------------------------------------------- #

# Define model (see above for references)

lat <- list(thr     = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items),
            neg_thr = dat_items$rr_neg_thr_items)

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_1f_1mf_neg_ce_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_1mf_neg_ce_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit <- cfa(ma_cfa_thr_1f_1mf_neg_ce_scen, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                               ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_summ <- summary(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig <- ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std
ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit)

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi[ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi_sig[order(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi_sig$mi,
                                                                                               decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 4, curvature = 2.5, reorder = FALSE,
           bifactor = "thr",
           latents = c("neg_thr", "thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 10, 1, 5),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:18, 2]  <- 1.5*pi
ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[19:36, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p$Arguments$edge.color == 
                                                                          "transparent"] <- "transparent"

# Make edges black

ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p$graphAttributes$Edges$color[ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p$graphAttributes$Edges$color == 
                                                                    "#808080FF"] <- "black"

plot(ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for bifactor all threat-items model with correlated errors per scenario (i.e., positive/negative factors uncorrelated) ----
# ---------------------------------------------------------------------------- #

# Define model (e.g., see Morin et al., 2020)

lat <- list(thr     = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items),
            pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items)

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_bf_ce_scen <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_bf_ce_scen)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_bf_ce_scen_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_bf_ce_scen_wlsmv_fit <- cfa(ma_cfa_thr_bf_ce_scen, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                                       ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_bf_ce_scen_wlsmv_summ <- summary(ma_cfa_thr_bf_ce_scen_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_bf_ce_scen_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_bf_ce_scen_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_bf_ce_scen_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_bf_ce_scen_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig <- ma_cfa_thr_bf_ce_scen_wlsmv_resid_std
ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig[abs(ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a",
                               "neg_thr_blood_test_9c"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_bf_ce_scen_wlsmv_mi <- modificationIndices(ma_cfa_thr_bf_ce_scen_wlsmv_fit)

ma_cfa_thr_bf_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_bf_ce_scen_wlsmv_mi[ma_cfa_thr_bf_ce_scen_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_bf_ce_scen_wlsmv_mi_sig <- ma_cfa_thr_bf_ce_scen_wlsmv_mi_sig[order(ma_cfa_thr_bf_ce_scen_wlsmv_mi_sig$mi,
                                                                                         decreasing = TRUE), ]

write.csv(ma_cfa_thr_bf_ce_scen_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_bf_ce_scen_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_bf_ce_scen_wlsmv_p <- 
  semPaths(ma_cfa_thr_bf_ce_scen_wlsmv_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 4, curvature = 2.5, reorder = FALSE,
           bifactor = "thr",
           latents = c("neg_thr", "pos_thr", "thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 5, 1, 5),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_bf_ce_scen_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_thr_bf_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:18, 2]  <- 1.5*pi
ma_cfa_thr_bf_ce_scen_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[19:36, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_thr_bf_ce_scen_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_thr_bf_ce_scen_wlsmv_p$Arguments$edge.color == 
                                                                "transparent"] <- "transparent"

# Make edges black

ma_cfa_thr_bf_ce_scen_wlsmv_p$graphAttributes$Edges$color[ma_cfa_thr_bf_ce_scen_wlsmv_p$graphAttributes$Edges$color == 
                                                                 "#808080FF"] <- "black"

plot(ma_cfa_thr_bf_ce_scen_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run CFA for higher-order all threat-items model with correlated errors per scenario -----
# ---------------------------------------------------------------------------- #

# Define model (for identification given only two first-order factors, constrain
# higher-order factor loadings to equality, which will yield a just-identified
# solution--see Brown, 2015, p. 292)

lat <- list(pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items,
            thr = "ho*pos_thr + ho*neg_thr")

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c",
            pos_thr_job_3a            = "neg_thr_job_3c",
            pos_thr_noise_4c          = "neg_thr_noise_4b",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c")

ma_cfa_thr_ho_ce <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_ho_ce)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_ho_ce_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_ho_ce_wlsmv_fit <- cfa(ma_cfa_thr_ho_ce, data = dat_ma_rr_only, std.lv = TRUE,
                                  ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_ho_ce_wlsmv_summ <- summary(ma_cfa_thr_ho_ce_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_ho_ce_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_ho_ce_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warning below (and higher-order factor loadings are 0)

# Warning message:
# In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = lavsamplestats,  :
#   lavaan WARNING:
#     The variance-covariance matrix of the estimated parameters (vcov)
#     does not appear to be positive definite! The smallest eigenvalue
#     (= 1.409122e-18) is close to zero. This may be a symptom that the
#     model is not identified.

raw_cov <- cov(dat_ma_rr_only[c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items)], 
               use = "pairwise.complete.obs")
ma_cfa_thr_ho_ce_wlsmv_fit_samp_cov <- lavInspect(ma_cfa_thr_ho_ce_wlsmv_fit, "sampstat")$cov
ma_cfa_thr_ho_ce_wlsmv_fit_vcov     <- lavInspect(ma_cfa_thr_ho_ce_wlsmv_fit, "vcov")
ma_cfa_thr_ho_ce_wlsmv_fit_cov_ov   <- lavInspect(ma_cfa_thr_ho_ce_wlsmv_fit, "cov.ov")
ma_cfa_thr_ho_ce_wlsmv_fit_cov_lv   <- lavInspect(ma_cfa_thr_ho_ce_wlsmv_fit, "cov.lv")

eigen(raw_cov)$values
eigen(ma_cfa_thr_ho_ce_wlsmv_fit_samp_cov)$values
eigen(ma_cfa_thr_ho_ce_wlsmv_fit_vcov)$values
all(eigen(ma_cfa_thr_ho_ce_wlsmv_fit_vcov)$values > 0) # TODO: In addition to warning, one is < 0
eigen(ma_cfa_thr_ho_ce_wlsmv_fit_cov_ov)$values
eigen(ma_cfa_thr_ho_ce_wlsmv_fit_cov_lv)$values

# Create plot without parameter estimates just to show structure

ma_cfa_thr_ho_ce_wlsmv_p <- 
  semPaths(ma_cfa_thr_ho_ce_wlsmv_fit, what = "path", whatLabels = "no", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_thr", "pos_thr", "thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 3, 1, 10),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_ho_ce_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Reverse score negative threat items (for testing) ----
# ---------------------------------------------------------------------------- #

dat_test <- dat_ma_rr_only[c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items)]

dat_test[dat_items$rr_neg_thr_items] <- 3 - dat_ma_rr_only[dat_items$rr_neg_thr_items]

dat_items$rr_neg_thr_items_rev <- paste0(dat_items$rr_neg_thr_items, "_rev")

names(dat_test)[names(dat_test) %in% dat_items$rr_neg_thr_items] <- dat_items$rr_neg_thr_items_rev

# ---------------------------------------------------------------------------- #
# Run CFA for 2 correlated-factors all threat-items model with correlated errors per scenario (reverse scored) ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_pos_thr_items,
            neg_thr = dat_items$rr_neg_thr_items_rev)

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a_rev",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c_rev",
            pos_thr_job_3a            = "neg_thr_job_3c_rev",
            pos_thr_noise_4c          = "neg_thr_noise_4b_rev",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d_rev",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a_rev",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a_rev",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a_rev",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c_rev")

ma_cfa_thr_2cf_ce_rev <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_2cf_ce_rev)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_rev_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_2cf_ce_rev_wlsmv_fit <- cfa(ma_cfa_thr_2cf_ce_rev, data = dat_test, std.lv = TRUE,
                                       ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_2cf_ce_rev_wlsmv_summ <- summary(ma_cfa_thr_2cf_ce_rev_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_2cf_ce_rev_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_2cf_ce_rev_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_2cf_ce_rev_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig <- ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std
ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig[abs(ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a_rev",
                               "neg_thr_blood_test_9c_rev"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_2cf_ce_rev_wlsmv_mi <- modificationIndices(ma_cfa_thr_2cf_ce_rev_wlsmv_fit)

ma_cfa_thr_2cf_ce_rev_wlsmv_mi_sig <- ma_cfa_thr_2cf_ce_rev_wlsmv_mi[ma_cfa_thr_2cf_ce_rev_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_2cf_ce_rev_wlsmv_mi_sig <- ma_cfa_thr_2cf_ce_rev_wlsmv_mi_sig[order(ma_cfa_thr_2cf_ce_rev_wlsmv_mi_sig$mi,
                                                                         decreasing = TRUE), ]

write.csv(ma_cfa_thr_2cf_ce_rev_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_2cf_ce_rev_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_2cf_ce_rev_wlsmv_p <- 
  semPaths(ma_cfa_thr_2cf_ce_rev_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_thr", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_thr_items_rev), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 20, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_2cf_ce_rev_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for 1 factor all threat-items model with correlated errors per scenario (reverse scored) ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(thr = c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items_rev))

cov <- list(pos_thr_elevator_1d       = "neg_thr_elevator_1a_rev",
            pos_thr_wedding_2a        = "neg_thr_wedding_2c_rev",
            pos_thr_job_3a            = "neg_thr_job_3c_rev",
            pos_thr_noise_4c          = "neg_thr_noise_4b_rev",
            pos_thr_meeting_friend_5c = "neg_thr_meeting_friend_5d_rev",
            pos_thr_lunch_6c          = "neg_thr_lunch_6a_rev",
            pos_thr_scrape_7c         = "neg_thr_scrape_7a_rev",
            pos_thr_shopping_8b       = "neg_thr_shopping_8a_rev",
            pos_thr_blood_test_9d     = "neg_thr_blood_test_9c_rev")

ma_cfa_thr_1f_ce_scen_rev <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_thr_1f_ce_scen_rev)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_scen_rev_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit <- cfa(ma_cfa_thr_1f_ce_scen_rev, data = dat_test, std.lv = TRUE,
                                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
cat("SUMMARY:", "\n\n")
(ma_cfa_thr_1f_ce_scen_rev_wlsmv_summ <- summary(ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
cat("COMPLETELY STANDARDIZED SOLUTION:", "\n\n")
(ma_cfa_thr_1f_ce_scen_rev_wlsmv_std_all <- standardizedsolution(ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit))
sink(type = "message")
sink()

# Inspect standardized residuals (focus on sig. ones; i.e., absolute value
# of standardized residual [z score] > 1.96 [sig. at .05], or in large samples use
# > 2.58 [sig. at .01]; see Brown, 2015, p. 99)

ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std <- lavResiduals(ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit, type = "raw")$cov.z

ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig <- ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std
ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig[abs(ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig) < 1.96] <- NA

pdf(file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig_plot.pdf"))
ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig_p <- 
  corrplot(ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig, 
           type = "lower", method = "shade", cl.pos = "n", addgrid.col = "grey", is.corr = FALSE,
           tl.col = "black", tl.srt = 45, tl.cex = .45,
           number.cex = .35, addCoef.col = "black", na.label.col = "transparent") |>
  corrRect(lwd = 1.5, name = c("pos_thr_elevator_1d", "neg_thr_elevator_1a_rev",
                               "neg_thr_blood_test_9c_rev"))
dev.off()

# Inspect modification indices (focus on sig. ones--i.e., modification index
# [chi-square difference with 1 df] > 3.84--that correspond to completely standardized 
# expected parameter change values ["sepc.all"] of nontrivial size; Brown, 2015, p. 102)

ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi <- modificationIndices(ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit)

ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi[ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi$mi > 3.84, ]
ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi_sig <- ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi_sig[order(ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi_sig$mi,
                                                                                 decreasing = TRUE), ]

write.csv(ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi_sig, row.names = FALSE,
          file = paste0(cfa_revised_path, "ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi_sig.csv"))

# Plot completely standardized solution

ma_cfa_thr_1f_ce_scen_rev_wlsmv_p <- 
  semPaths(ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 4, curvature = 2, cardinal = TRUE, reorder = FALSE,
           latents = "thr",
           manifests = c(rev(dat_items$rr_neg_thr_items_rev), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 4, 1, 22),
           filetype = "pdf", filename = "./results/cfa/revised/ma_cfa_thr_1f_ce_scen_rev_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Export results ----
# ---------------------------------------------------------------------------- #

res_ma_cfa_28 <- 
  list(ma_cfa_28_3cf                                                   = ma_cfa_28_3cf,
       ma_cfa_28_3cf_pml_ac_fit                                        = ma_cfa_28_3cf_pml_ac_fit,
       ma_cfa_28_3cf_pml_ac_summ                                       = ma_cfa_28_3cf_pml_ac_summ,
       ma_cfa_28_3cf_pml_ac_std_all                                    = ma_cfa_28_3cf_pml_ac_std_all,
       ma_cfa_28_3cf_pml_ac_resid_cor                                  = ma_cfa_28_3cf_pml_ac_resid_cor,
       ma_cfa_28_3cf_pml_ac_resid_cor_small_p                          = ma_cfa_28_3cf_pml_ac_resid_cor_small_p,
       ma_cfa_28_3cf_pml_ac_p                                          = ma_cfa_28_3cf_pml_ac_p,
       ma_cfa_28_3cf_wlsmv_fit                                         = ma_cfa_28_3cf_wlsmv_fit,
       ma_cfa_28_3cf_wlsmv_summ                                        = ma_cfa_28_3cf_wlsmv_summ,
       ma_cfa_28_3cf_wlsmv_std_all                                     = ma_cfa_28_3cf_wlsmv_std_all,
       ma_cfa_28_3cf_wlsmv_dfi                                         = ma_cfa_28_3cf_wlsmv_dfi,
       ma_cfa_28_3cf_wlsmv_resid_std                                   = ma_cfa_28_3cf_wlsmv_resid_std,
       ma_cfa_28_3cf_wlsmv_resid_std_sig_p                             = ma_cfa_28_3cf_wlsmv_resid_std_sig_p,
       ma_cfa_28_3cf_wlsmv_mi                                          = ma_cfa_28_3cf_wlsmv_mi,
       ma_cfa_28_3cf_wlsmv_resid_cor                                   = ma_cfa_28_3cf_wlsmv_resid_cor,
       ma_cfa_28_3cf_wlsmv_resid_cor_small_p                           = ma_cfa_28_3cf_wlsmv_resid_cor_small_p,
       ma_cfa_28_3cf_wlsmv_p                                           = ma_cfa_28_3cf_wlsmv_p,
       ma_cfa_28_2cf                                                   = ma_cfa_28_2cf,
       ma_cfa_28_2cf_wlsmv_fit                                         = ma_cfa_28_2cf_wlsmv_fit,
       ma_cfa_28_2cf_wlsmv_summ                                        = ma_cfa_28_2cf_wlsmv_summ,
       ma_cfa_28_2cf_wlsmv_std_all                                     = ma_cfa_28_2cf_wlsmv_std_all,
       ma_cfa_28_2cf_wlsmv_dfi                                         = ma_cfa_28_2cf_wlsmv_dfi,
       ma_cfa_28_2cf_wlsmv_resid_std                                   = ma_cfa_28_2cf_wlsmv_resid_std,
       ma_cfa_28_2cf_wlsmv_resid_std_sig_p                             = ma_cfa_28_2cf_wlsmv_resid_std_sig_p,
       ma_cfa_28_2cf_wlsmv_mi                                          = ma_cfa_28_2cf_wlsmv_mi,
       ma_cfa_28_2cf_wlsmv_p                                           = ma_cfa_28_2cf_wlsmv_p,
       ma_cfa_28_3cf_efa                                               = ma_cfa_28_3cf_efa,
       ma_cfa_28_3cf_efa_wlsmv_fit                                     = ma_cfa_28_3cf_efa_wlsmv_fit,
       ma_cfa_28_3cf_efa_wlsmv_summ                                    = ma_cfa_28_3cf_efa_wlsmv_summ,
       ma_cfa_28_3cf_efa_wlsmv_std_all                                 = ma_cfa_28_3cf_efa_wlsmv_std_all,
       ma_cfa_28_3cf_efa_wlsmv_resid_std                               = ma_cfa_28_3cf_efa_wlsmv_resid_std,
       ma_cfa_28_3cf_efa_wlsmv_resid_std_sig_p                         = ma_cfa_28_3cf_efa_wlsmv_resid_std_sig_p,
       ma_cfa_28_3cf_efa_wlsmv_mi                                      = ma_cfa_28_3cf_efa_wlsmv_mi,
       ma_cfa_28_3cf_efa_wlsmv_p                                       = ma_cfa_28_3cf_efa_wlsmv_p,
       ma_cfa_28_3cf_ce_thr                                            = ma_cfa_28_3cf_ce_thr,
       ma_cfa_28_3cf_ce_thr_wlsmv_fit                                  = ma_cfa_28_3cf_ce_thr_wlsmv_fit,
       ma_cfa_28_3cf_ce_thr_wlsmv_summ                                 = ma_cfa_28_3cf_ce_thr_wlsmv_summ,
       ma_cfa_28_3cf_ce_thr_wlsmv_std_all                              = ma_cfa_28_3cf_ce_thr_wlsmv_std_all,
       ma_cfa_28_3cf_ce_thr_wlsmv_dfi                                  = ma_cfa_28_3cf_ce_thr_wlsmv_dfi,
       ma_cfa_28_3cf_ce_thr_wlsmv_resid_std                            = ma_cfa_28_3cf_ce_thr_wlsmv_resid_std,
       ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig_p                      = ma_cfa_28_3cf_ce_thr_wlsmv_resid_std_sig_p,
       ma_cfa_28_3cf_ce_thr_wlsmv_mi                                   = ma_cfa_28_3cf_ce_thr_wlsmv_mi,
       ma_cfa_28_3cf_ce_thr_wlsmv_p                                    = ma_cfa_28_3cf_ce_thr_wlsmv_p,
       ma_cfa_28_3cf_ce_all                                            = ma_cfa_28_3cf_ce_all,
       ma_cfa_28_3cf_ce_all_wlsmv_fit                                  = ma_cfa_28_3cf_ce_all_wlsmv_fit,
       ma_cfa_28_3cf_ce_all_wlsmv_summ                                 = ma_cfa_28_3cf_ce_all_wlsmv_summ,
       ma_cfa_28_3cf_ce_all_wlsmv_std_all                              = ma_cfa_28_3cf_ce_all_wlsmv_std_all,
       ma_cfa_28_3cf_ce_all_wlsmv_dfi                                  = ma_cfa_28_3cf_ce_all_wlsmv_dfi,
       ma_cfa_28_3cf_ce_all_wlsmv_resid_std                            = ma_cfa_28_3cf_ce_all_wlsmv_resid_std,
       ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig_p                      = ma_cfa_28_3cf_ce_all_wlsmv_resid_std_sig_p,
       ma_cfa_28_3cf_ce_all_wlsmv_mi                                   = ma_cfa_28_3cf_ce_all_wlsmv_mi,
       ma_cfa_28_3cf_ce_all_wlsmv_p                                    = ma_cfa_28_3cf_ce_all_wlsmv_p,
       ma_cfa_28_2cf_4cmf                                              = ma_cfa_28_2cf_4cmf,
       ma_cfa_28_2cf_4cmf_wlsmv_fit                                    = ma_cfa_28_2cf_4cmf_wlsmv_fit,
       ma_cfa_28_2cf_4cmf_wlsmv_summ                                   = ma_cfa_28_2cf_4cmf_wlsmv_summ,
       ma_cfa_28_2cf_4cmf_wlsmv_std_all                                = ma_cfa_28_2cf_4cmf_wlsmv_std_all,
       ma_cfa_28_2cf_4cmf_wlsmv_p                                      = ma_cfa_28_2cf_4cmf_wlsmv_p,
       ma_cfa_28_2cf_2cmf_thr                                          = ma_cfa_28_2cf_2cmf_thr,
       ma_cfa_28_2cf_2cmf_thr_wlsmv_fit                                = ma_cfa_28_2cf_2cmf_thr_wlsmv_fit,
       ma_cfa_28_2cf_2cmf_thr_wlsmv_summ                               = ma_cfa_28_2cf_2cmf_thr_wlsmv_summ,
       ma_cfa_28_2cf_2cmf_thr_wlsmv_std_all                            = ma_cfa_28_2cf_2cmf_thr_wlsmv_std_all,
       ma_cfa_28_2cf_2cmf_thr_wlsmv_p                                  = ma_cfa_28_2cf_2cmf_thr_wlsmv_p,
       ma_cfa_28_2cf_2mf_thr                                           = ma_cfa_28_2cf_2mf_thr,
       ma_cfa_28_2cf_2mf_thr_wlsmv_fit                                 = ma_cfa_28_2cf_2mf_thr_wlsmv_fit,
       ma_cfa_28_2cf_2mf_thr_wlsmv_summ                                = ma_cfa_28_2cf_2mf_thr_wlsmv_summ,
       ma_cfa_28_2cf_2mf_thr_wlsmv_std_all                             = ma_cfa_28_2cf_2mf_thr_wlsmv_std_all,
       ma_cfa_28_2cf_2mf_thr_wlsmv_p                                   = ma_cfa_28_2cf_2mf_thr_wlsmv_p,
       ma_cfa_28_2cf_thr_non_ce_val_scen_thr                           = ma_cfa_28_2cf_thr_non_ce_val_scen_thr,
       ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit                 = ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_fit,
       ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_summ                = ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_summ,
       ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_std_all             = ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_std_all,
       ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_p                   = ma_cfa_28_2cf_thr_non_ce_val_scen_thr_wlsmv_p,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr                       = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit             = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_fit,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_summ            = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_summ,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_std_all         = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_std_all,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std       = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig_p = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_resid_std_sig_p,
       ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_p               = ma_cfa_28_2cf_thr_non_ce_pos_val_scen_thr_wlsmv_p)

res_ma_cfa_thr <-
  list(ma_cfa_thr_2cf                                      = ma_cfa_thr_2cf,
       ma_cfa_thr_2cf_wlsmv_fit                            = ma_cfa_thr_2cf_wlsmv_fit,
       ma_cfa_thr_2cf_wlsmv_summ                           = ma_cfa_thr_2cf_wlsmv_summ,
       ma_cfa_thr_2cf_wlsmv_std_all                        = ma_cfa_thr_2cf_wlsmv_std_all,
       ma_cfa_thr_2cf_wlsmv_resid_std                      = ma_cfa_thr_2cf_wlsmv_resid_std ,
       ma_cfa_thr_2cf_wlsmv_resid_std_sig_p                = ma_cfa_thr_2cf_wlsmv_resid_std_sig_p,
       ma_cfa_thr_2cf_wlsmv_mi                             = ma_cfa_thr_2cf_wlsmv_mi,
       ma_cfa_thr_2cf_wlsmv_p                              = ma_cfa_thr_2cf_wlsmv_p,
       ma_cfa_thr_2cf_efa                                  = ma_cfa_thr_2cf_efa,
       ma_cfa_thr_2cf_efa_wlsmv_fit                        = ma_cfa_thr_2cf_efa_wlsmv_fit,
       ma_cfa_thr_2cf_efa_wlsmv_summ                       = ma_cfa_thr_2cf_efa_wlsmv_summ,
       ma_cfa_thr_2cf_efa_wlsmv_std_all                    = ma_cfa_thr_2cf_efa_wlsmv_std_all,
       ma_cfa_thr_2cf_efa_wlsmv_resid_std                  = ma_cfa_thr_2cf_efa_wlsmv_resid_std,
       ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig_p            = ma_cfa_thr_2cf_efa_wlsmv_resid_std_sig_p,
       ma_cfa_thr_2cf_efa_wlsmv_mi                         = ma_cfa_thr_2cf_efa_wlsmv_mi,
       ma_cfa_thr_2cf_efa_wlsmv_p                          = ma_cfa_thr_2cf_efa_wlsmv_p,
       ma_cfa_thr_2cf_ce                                   = ma_cfa_thr_2cf_ce,
       ma_cfa_thr_2cf_ce_wlsmv_fit                         = ma_cfa_thr_2cf_ce_wlsmv_fit,
       ma_cfa_thr_2cf_ce_wlsmv_summ                        = ma_cfa_thr_2cf_ce_wlsmv_summ,
       ma_cfa_thr_2cf_ce_wlsmv_std_all                     = ma_cfa_thr_2cf_ce_wlsmv_std_all,
       ma_cfa_thr_2cf_ce_wlsmv_dfi                         = ma_cfa_thr_2cf_ce_wlsmv_dfi,
       ma_cfa_thr_2cf_ce_wlsmv_resid_std                   = ma_cfa_thr_2cf_ce_wlsmv_resid_std,
       ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig_p             = ma_cfa_thr_2cf_ce_wlsmv_resid_std_sig_p,
       ma_cfa_thr_2cf_ce_wlsmv_mi                          = ma_cfa_thr_2cf_ce_wlsmv_mi,
       ma_cfa_thr_2cf_ce_wlsmv_p                           = ma_cfa_thr_2cf_ce_wlsmv_p,
       ma_cfa_thr_1f_ce_scen                               = ma_cfa_thr_1f_ce_scen,
       ma_cfa_thr_1f_ce_scen_wlsmv_fit                     = ma_cfa_thr_1f_ce_scen_wlsmv_fit,
       ma_cfa_thr_1f_ce_scen_wlsmv_summ                    = ma_cfa_thr_1f_ce_scen_wlsmv_summ,
       ma_cfa_thr_1f_ce_scen_wlsmv_std_all                 = ma_cfa_thr_1f_ce_scen_wlsmv_std_all,
       ma_cfa_thr_1f_ce_scen_wlsmv_resid_std               = ma_cfa_thr_1f_ce_scen_wlsmv_resid_std,
       ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig_p         = ma_cfa_thr_1f_ce_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_ce_scen_wlsmv_mi                      = ma_cfa_thr_1f_ce_scen_wlsmv_mi,
       ma_cfa_thr_1f_ce_scen_wlsmv_p                       = ma_cfa_thr_1f_ce_scen_wlsmv_p,
       ma_cfa_thr_1f_ce_val_scen                           = ma_cfa_thr_1f_ce_val_scen,
       ma_cfa_thr_1f_ce_val_scen_wlsmv_fit                 = ma_cfa_thr_1f_ce_val_scen_wlsmv_fit,
       ma_cfa_thr_1f_ce_val_scen_wlsmv_summ                = ma_cfa_thr_1f_ce_val_scen_wlsmv_summ,
       ma_cfa_thr_1f_ce_val_scen_wlsmv_std_all             = ma_cfa_thr_1f_ce_val_scen_wlsmv_std_all,
       ma_cfa_thr_1f_ce_val_scen_wlsmv_p                   = ma_cfa_thr_1f_ce_val_scen_wlsmv_p,
       ma_cfa_thr_1f_ce_pos_scen                           = ma_cfa_thr_1f_ce_pos_scen,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit                 = ma_cfa_thr_1f_ce_pos_scen_wlsmv_fit,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_summ                = ma_cfa_thr_1f_ce_pos_scen_wlsmv_summ,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_std_all             = ma_cfa_thr_1f_ce_pos_scen_wlsmv_std_all,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std           = ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig_p     = ma_cfa_thr_1f_ce_pos_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi                  = ma_cfa_thr_1f_ce_pos_scen_wlsmv_mi,
       ma_cfa_thr_1f_ce_pos_scen_wlsmv_p                   = ma_cfa_thr_1f_ce_pos_scen_wlsmv_p,
       ma_cfa_thr_1f_ce_neg_scen                           = ma_cfa_thr_1f_ce_neg_scen,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit                 = ma_cfa_thr_1f_ce_neg_scen_wlsmv_fit,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_summ                = ma_cfa_thr_1f_ce_neg_scen_wlsmv_summ,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_std_all             = ma_cfa_thr_1f_ce_neg_scen_wlsmv_std_all,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std           = ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig_p     = ma_cfa_thr_1f_ce_neg_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi                  = ma_cfa_thr_1f_ce_neg_scen_wlsmv_mi,
       ma_cfa_thr_1f_ce_neg_scen_wlsmv_p                   = ma_cfa_thr_1f_ce_neg_scen_wlsmv_p,
       ma_cfa_thr_1f_ce_val                                = ma_cfa_thr_1f_ce_val,
       ma_cfa_thr_1f_ce_val_wlsmv_fit                      = ma_cfa_thr_1f_ce_val_wlsmv_fit,
       ma_cfa_thr_1f_ce_val_wlsmv_summ                     = ma_cfa_thr_1f_ce_val_wlsmv_summ,
       ma_cfa_thr_1f_ce_val_wlsmv_std_all                  = ma_cfa_thr_1f_ce_val_wlsmv_std_all,
       ma_cfa_thr_1f_ce_val_wlsmv_p                        = ma_cfa_thr_1f_ce_val_wlsmv_p,
       ma_cfa_thr_1f_2cmf                                  = ma_cfa_thr_1f_2cmf,
       ma_cfa_thr_1f_2cmf_wlsmv_fit                        = ma_cfa_thr_1f_2cmf_wlsmv_fit,
       ma_cfa_thr_1f_2cmf_wlsmv_summ                       = ma_cfa_thr_1f_2cmf_wlsmv_summ,
       ma_cfa_thr_1f_2cmf_wlsmv_std_all                    = ma_cfa_thr_1f_2cmf_wlsmv_std_all,
       ma_cfa_thr_1f_2cmf_wlsmv_resid_std                  = ma_cfa_thr_1f_2cmf_wlsmv_resid_std,
       ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig_p            = ma_cfa_thr_1f_2cmf_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_2cmf_wlsmv_mi                         = ma_cfa_thr_1f_2cmf_wlsmv_mi,
       ma_cfa_thr_1f_2cmf_wlsmv_p                          = ma_cfa_thr_1f_2cmf_wlsmv_p,
       ma_cfa_thr_1f_2cmf_ce_scen                          = ma_cfa_thr_1f_2cmf_ce_scen,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit                = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_fit,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_summ               = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_summ,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_std_all            = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_std_all,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std          = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig_p    = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi                 = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_mi,
       ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p                  = ma_cfa_thr_1f_2cmf_ce_scen_wlsmv_p,
       ma_cfa_thr_1f_1mf_pos_ce_scen                       = ma_cfa_thr_1f_1mf_pos_ce_scen,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit             = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_fit,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_summ            = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_summ,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_std_all         = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_std_all,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std       = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig_p = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi              = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_mi,
       ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p               = ma_cfa_thr_1f_1mf_pos_ce_scen_wlsmv_p,
       ma_cfa_thr_1f_1mf_neg_ce_scen                       = ma_cfa_thr_1f_1mf_neg_ce_scen,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit             = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_fit,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_summ            = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_summ,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_std_all         = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_std_all,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std       = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig_p = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi              = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_mi,
       ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p               = ma_cfa_thr_1f_1mf_neg_ce_scen_wlsmv_p,
       ma_cfa_thr_bf_ce_scen                               = ma_cfa_thr_bf_ce_scen,
       ma_cfa_thr_bf_ce_scen_wlsmv_fit                     = ma_cfa_thr_bf_ce_scen_wlsmv_fit,
       ma_cfa_thr_bf_ce_scen_wlsmv_summ                    = ma_cfa_thr_bf_ce_scen_wlsmv_summ,
       ma_cfa_thr_bf_ce_scen_wlsmv_std_all                 = ma_cfa_thr_bf_ce_scen_wlsmv_std_all,
       ma_cfa_thr_bf_ce_scen_wlsmv_resid_std               = ma_cfa_thr_bf_ce_scen_wlsmv_resid_std,
       ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig_p         = ma_cfa_thr_bf_ce_scen_wlsmv_resid_std_sig_p,
       ma_cfa_thr_bf_ce_scen_wlsmv_mi                      = ma_cfa_thr_bf_ce_scen_wlsmv_mi,
       ma_cfa_thr_bf_ce_scen_wlsmv_p                       = ma_cfa_thr_bf_ce_scen_wlsmv_p,
       ma_cfa_thr_ho_ce                                    = ma_cfa_thr_ho_ce,
       ma_cfa_thr_ho_ce_wlsmv_fit                          = ma_cfa_thr_ho_ce_wlsmv_fit,
       ma_cfa_thr_ho_ce_wlsmv_summ                         = ma_cfa_thr_ho_ce_wlsmv_summ,
       ma_cfa_thr_ho_ce_wlsmv_std_all                      = ma_cfa_thr_ho_ce_wlsmv_std_all,
       ma_cfa_thr_ho_ce_wlsmv_p                            = ma_cfa_thr_ho_ce_wlsmv_p,
       ma_cfa_thr_2cf_ce_rev                               = ma_cfa_thr_2cf_ce_rev,
       ma_cfa_thr_2cf_ce_rev_wlsmv_fit                     = ma_cfa_thr_2cf_ce_rev_wlsmv_fit,
       ma_cfa_thr_2cf_ce_rev_wlsmv_summ                    = ma_cfa_thr_2cf_ce_rev_wlsmv_summ,
       ma_cfa_thr_2cf_ce_rev_wlsmv_std_all                 = ma_cfa_thr_2cf_ce_rev_wlsmv_std_all,
       ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std               = ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std,
       ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig_p         = ma_cfa_thr_2cf_ce_rev_wlsmv_resid_std_sig_p,
       ma_cfa_thr_2cf_ce_rev_wlsmv_mi                      = ma_cfa_thr_2cf_ce_rev_wlsmv_mi,
       ma_cfa_thr_2cf_ce_rev_wlsmv_p                       = ma_cfa_thr_2cf_ce_rev_wlsmv_p,
       ma_cfa_thr_1f_ce_scen_rev                           = ma_cfa_thr_1f_ce_scen_rev,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit                 = ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_summ                = ma_cfa_thr_1f_ce_scen_rev_wlsmv_summ,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_std_all             = ma_cfa_thr_1f_ce_scen_rev_wlsmv_std_all,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std           = ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig_p     = ma_cfa_thr_1f_ce_scen_rev_wlsmv_resid_std_sig_p,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi                  = ma_cfa_thr_1f_ce_scen_rev_wlsmv_mi,
       ma_cfa_thr_1f_ce_scen_rev_wlsmv_p                   = ma_cfa_thr_1f_ce_scen_rev_wlsmv_p)

# Save results

save(res_ma_cfa_28,  file = paste0(cfa_revised_path, "res_ma_cfa_28.RData"))
save(res_ma_cfa_thr, file = paste0(cfa_revised_path, "res_ma_cfa_thr.RData"))