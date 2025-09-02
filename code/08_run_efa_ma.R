# ---------------------------------------------------------------------------- #
# Run EFA in Managing Anxiety Study
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# Note: Random seed must be set for each EFA analysis for reproducible results

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

pkgs <- c("psych", "lavaan", "lavaanExtra")
groundhog.library(pkgs, groundhog_day)

# Set seed

set.seed(1234)

# ---------------------------------------------------------------------------- #
# Define functions used in script ----
# ---------------------------------------------------------------------------- #

# Define function to export basic EFA results and details to TXT and loadings to CSV

export_efa_res <- function(fit, path, filename_stem) {
  sink(paste0(path, paste0(filename_stem, ".txt")))
  print(summary(fit))
  sink()
  
  sink(paste0(path, paste0(filename_stem, "_detail.txt")))
  print(summary(fit, se = TRUE, zstat = TRUE, pvalue = TRUE))
  sink()
  
  sink(paste0(path, paste0(filename_stem, ".csv")))
  print(fit$loadings)
  sink()
}

# ---------------------------------------------------------------------------- #
# Import data and create EFA results path ----
# ---------------------------------------------------------------------------- #

load("./data/ma/final/dat_ma_rr_only.RData")

load("./data/helper/dat_items.RData")
load("./data/helper/rr_item_map.RData")

efa_path <- "./results/efa/"
dir.create(efa_path)

# ---------------------------------------------------------------------------- #
# Prepare data ----
# ---------------------------------------------------------------------------- #

# Order columns in same order as prior CFAs

target_order <- c(dat_items$rr_pos_thr_items, dat_items$rr_pos_non_items,
                  dat_items$rr_neg_thr_items, dat_items$rr_neg_non_items)

dat_ma_rr_only <- dat_ma_rr_only[match(target_order, names(dat_ma_rr_only))]

# ---------------------------------------------------------------------------- #
# Inspect item distributions ----
# ---------------------------------------------------------------------------- #

# Define function to plot histograms for six items at a time

plot_item_hists <- function(df, cols) {
  par(mfrow = c(3, 2))
  
  for (i in cols) {
    col_name <- names(df[i])
    hist(df[, i], main = col_name, xlab = "")
  }
}

# Run function and export plots to PDF

hists_path <- paste0(efa_path, "hists/")

dir.create(hists_path, recursive = TRUE)

pdf(paste0(hists_path, "ma_rr_hists.pdf"), height = 6, width = 6)
plot_item_hists(dat_ma_rr_only, 1:6)
plot_item_hists(dat_ma_rr_only, 7:12)
plot_item_hists(dat_ma_rr_only, 13:18)
plot_item_hists(dat_ma_rr_only, 19:24)
plot_item_hists(dat_ma_rr_only, 25:30)
plot_item_hists(dat_ma_rr_only, 31:36)
dev.off()

# Note: Some item distributions appear heavily skewed

# ---------------------------------------------------------------------------- #
# Determine number of factors for all items ----
# ---------------------------------------------------------------------------- #

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

all_items_path <- paste0(efa_path,       "all_items/")
scree_path     <- paste0(all_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only, factors = FALSE, pc = TRUE)
dev.off()

# Also consider parallel analysis based on principal component analysis (PA-PCA-m),
# which Lim and Jahng (2019) found is better than principal axis factoring (PA-PAF-m)
# across a wide variety of situations (inc. ordinal data) and recommend. PA-PCA-m
# using pairwise complete data suggests an upper bound (see Montoya & Edwards, 2021, 
# p. 416) of 4 components.

pdf(paste0(scree_path, "ma_rr_scree_pa_ml.pdf"), height = 6, width = 6)
(result_ml <- fa.parallel(dat_ma_rr_only, fm = "ml", n.iter = 100))
dev.off()

result_ml$ncomp == 4

sum(result_ml$pc.values > result_ml$pc.sim)  == 4
sum(result_ml$pc.values > result_ml$pc.simr) == 4

# Also do parallel analysis of polychoric correlations (given that RR items are ordinal; 
# and some are heavily skewed--see above). (Note: Per personal correspondence on 9/17/2023, 
# William Revelle suggested setting "correct = 0" to avoid error/warnings when setting 'cor 
# = "poly"' [though no warnings here so far here].) Unclear what estimation method to use for 
# polychoric correlations, so try both "ML" (as above, but with only 10 iterations as it 
# takes long to run) and "minres" (default, incl. for "fa.parallel.poly()"). Same results:
# PA-PCA-m using pairwise complete data suggests an upper bound of 5 components.

pdf(paste0(scree_path, "ma_rr_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_poly_ml <- fa.parallel(dat_ma_rr_only, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_poly_ml$ncomp == 5

sum(result_poly_ml$pc.values > result_poly_ml$pc.sim)  == 5 # Based on simulated data
sum(result_poly_ml$pc.values > result_poly_ml$pc.simr) == 3 # Based on resampled data

pdf(paste0(scree_path, "ma_rr_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_poly_minres <- fa.parallel(dat_ma_rr_only, correct = 0, cor = "poly"))
dev.off()

result_poly_minres$ncomp == 5

sum(result_poly_minres$pc.values > result_poly_minres$pc.sim)  == 5
sum(result_poly_minres$pc.values > result_poly_minres$pc.simr) == 4

# Based on scree plot, considered 4 factors, but parallel analysis suggests up to 5 (though 
# looking at +/- 1 is recommended; Lim & Jahng, 2019). Thus, consider 2, 3, 4, 5, or 6.

# ---------------------------------------------------------------------------- #
# Run EFA based on all items using PML-AC ----
# ---------------------------------------------------------------------------- #

# TODO: Negative chi-square and df and CFI of 1 (for oblimin, geomin, promax; 
# except for 3 factors). Otherwise no Heywood cases. Email lavaan group and authors
# of PML-AC paper, and may need to try a different estimator (e.g., WLSMV).





set.seed(1234)
ma_efa_oblimin_pml_ac_fit <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "oblimin", 
                                 ordered = TRUE, estimator = "PML", missing = "available.cases")
set.seed(1234)
ma_efa_geomin_pml_ac_fit  <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "geomin",  
                                 ordered = TRUE, estimator = "PML", missing = "available.cases")
set.seed(1234)
ma_efa_promax_pml_ac_fit  <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "promax",  
                                 ordered = TRUE, estimator = "PML", missing = "available.cases")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_oblimin_pml_ac_fit, all_items_path, "ma_efa_oblimin_pml_ac")
export_efa_res(ma_efa_geomin_pml_ac_fit,  all_items_path, "ma_efa_geomin_pml_ac")
export_efa_res(ma_efa_promax_pml_ac_fit,  all_items_path, "ma_efa_promax_pml_ac")

# ---------------------------------------------------------------------------- #
# Run EFA based on all items using FIML ----
# ---------------------------------------------------------------------------- #

# TODO: Some SEs seem too large (for oblimin, huge for 2 factor, some large for 4 
# factor, some gigantic for 6 factor; for geomin, some large for 4 factor), possibly
# due to using continuous estimator on ordinal data (or model misspecification; see
# Brown, 2015, p. 114). Otherwise no Heywood cases.





set.seed(1234)
ma_efa_oblimin_fiml_fit <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "oblimin", 
                               estimator = "MLR", missing = "fiml")
set.seed(1234)
ma_efa_geomin_fiml_fit  <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "geomin",  
                               estimator = "MLR", missing = "fiml")
set.seed(1234)
ma_efa_promax_fiml_fit  <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "promax",  
                               estimator = "MLR", missing = "fiml")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_oblimin_fiml_fit, all_items_path, "ma_efa_oblimin_fiml")
export_efa_res(ma_efa_geomin_fiml_fit,  all_items_path, "ma_efa_geomin_fiml")
export_efa_res(ma_efa_promax_fiml_fit,  all_items_path, "ma_efa_promax_fiml")

# ---------------------------------------------------------------------------- #
# Run EFA based on all items using WLSMV ----
# ---------------------------------------------------------------------------- #

# Given fit statistics issues with PML-AC, try WLSMV. No Heywood cases.

set.seed(1234)
ma_efa_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "oblimin", 
                                ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "geomin", 
                                ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only, nfactors = 2:6, rotation = "promax", 
                                ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_oblimin_wlsmv_fit, all_items_path, "ma_efa_oblimin_wlsmv")
export_efa_res(ma_efa_geomin_wlsmv_fit,  all_items_path, "ma_efa_geomin_wlsmv")
export_efa_res(ma_efa_promax_wlsmv_fit,  all_items_path, "ma_efa_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 35 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding item that consistently does not load with any factors in the 
# 2-4-factor solutions (which are more interpretable than the 5-6-factor solutions)

dat_ma_rr_only_35 <- dat_ma_rr_only[, names(dat_ma_rr_only)[names(dat_ma_rr_only) != 
                                                              "pos_thr_wedding_2a"]]

length(dat_ma_rr_only_35) == 35

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_35, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_35_items_path <- paste0(efa_path, "red_35_items/")
scree_path        <- paste0(red_35_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_35_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_35, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_35_scree_pa_ml.pdf"), height = 6, width = 6)
(result_35_ml <- fa.parallel(dat_ma_rr_only_35, fm = "ml", n.iter = 100))
dev.off()

result_35_ml$ncomp == 3

sum(result_35_ml$pc.values > result_35_ml$pc.sim)  == 3
sum(result_35_ml$pc.values > result_35_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 5 components

pdf(paste0(scree_path, "ma_rr_35_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_35_poly_ml <- fa.parallel(dat_ma_rr_only_35, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_35_poly_ml$ncomp == 4

sum(result_35_poly_ml$pc.values > result_35_poly_ml$pc.sim)  == 4 # Based on simulated data
sum(result_35_poly_ml$pc.values > result_35_poly_ml$pc.simr) == 3 # Based on resampled data

pdf(paste0(scree_path, "ma_rr_35_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_35_poly_minres <- fa.parallel(dat_ma_rr_only_35, correct = 0, cor = "poly"))
dev.off()

result_35_poly_minres$ncomp == 3

sum(result_35_poly_minres$pc.values > result_35_poly_minres$pc.sim)  == 4
sum(result_35_poly_minres$pc.values > result_35_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on 35 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_35_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_35, nfactors = 2:5, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_35_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_35, nfactors = 2:5, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_35_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_35, nfactors = 2:5, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_35_oblimin_wlsmv_fit, red_35_items_path, "ma_efa_35_oblimin_wlsmv")
export_efa_res(ma_efa_35_geomin_wlsmv_fit,  red_35_items_path, "ma_efa_35_geomin_wlsmv")
export_efa_res(ma_efa_35_promax_wlsmv_fit,  red_35_items_path, "ma_efa_35_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 34 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding negative threat item that consistently cross-loads onto nonthreat
# or positive threat factors

dat_ma_rr_only_34 <- dat_ma_rr_only_35[, names(dat_ma_rr_only_35)[names(dat_ma_rr_only_35) != 
                                                                    "neg_thr_scrape_7a"]]

length(dat_ma_rr_only_34) == 34

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_34, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_34_items_path <- paste0(efa_path, "red_34_items/")
scree_path        <- paste0(red_34_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_34_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_34, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_34_scree_pa_ml.pdf"), height = 6, width = 6)
(result_34_ml <- fa.parallel(dat_ma_rr_only_34, fm = "ml", n.iter = 100))
dev.off()

result_34_ml$ncomp == 3

sum(result_34_ml$pc.values > result_34_ml$pc.sim)  == 3
sum(result_34_ml$pc.values > result_34_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 4 components

pdf(paste0(scree_path, "ma_rr_34_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_34_poly_ml <- fa.parallel(dat_ma_rr_only_34, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_34_poly_ml$ncomp == 4

sum(result_34_poly_ml$pc.values > result_34_poly_ml$pc.sim)  == 4 # Based on simulated data
sum(result_34_poly_ml$pc.values > result_34_poly_ml$pc.simr) == 3 # Based on resampled data

# Warnings (two kinds) for "minres":
# In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   The estimated weights for the factor scores are probably incorrect.  Try a different factor score estimation method.
# In fac(r = r, nfactors = nfactors, n.obs = n.obs, rotate = rotate,  :
#   An ultra-Heywood case was detected.  Examine the results carefully

pdf(paste0(scree_path, "ma_rr_34_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_34_poly_minres <- fa.parallel(dat_ma_rr_only_34, correct = 0, cor = "poly"))
dev.off()

result_34_poly_minres$ncomp == 3

sum(result_34_poly_minres$pc.values > result_34_poly_minres$pc.sim)  == 4
sum(result_34_poly_minres$pc.values > result_34_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on 34 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_34_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_34, nfactors = 2:5, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_34_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_34, nfactors = 2:5, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_34_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_34, nfactors = 2:5, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_34_oblimin_wlsmv_fit, red_34_items_path, "ma_efa_34_oblimin_wlsmv")
export_efa_res(ma_efa_34_geomin_wlsmv_fit,  red_34_items_path, "ma_efa_34_geomin_wlsmv")
export_efa_res(ma_efa_34_promax_wlsmv_fit,  red_34_items_path, "ma_efa_34_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 33 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding positive threat item that consistently has weak loadings on
# positive threat factor and loads (negatively) onto negative threat item factor

dat_ma_rr_only_33 <- dat_ma_rr_only_34[, names(dat_ma_rr_only_34)[names(dat_ma_rr_only_34) != 
                                                                    "pos_thr_elevator_1d"]]

length(dat_ma_rr_only_33) == 33

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_33, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_33_items_path <- paste0(efa_path, "red_33_items/")
scree_path        <- paste0(red_33_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_33_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_33, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_33_scree_pa_ml.pdf"), height = 6, width = 6)
(result_33_ml <- fa.parallel(dat_ma_rr_only_33, fm = "ml", n.iter = 100))
dev.off()

result_33_ml$ncomp == 3

sum(result_33_ml$pc.values > result_33_ml$pc.sim)  == 3
sum(result_33_ml$pc.values > result_33_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 4 components

pdf(paste0(scree_path, "ma_rr_33_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_33_poly_ml <- fa.parallel(dat_ma_rr_only_33, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_33_poly_ml$ncomp == 3

sum(result_33_poly_ml$pc.values > result_33_poly_ml$pc.sim)  == 3 # Based on simulated data
sum(result_33_poly_ml$pc.values > result_33_poly_ml$pc.simr) == 3 # Based on resampled data

# Warnings (one kind) for "minres":
# In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   The estimated weights for the factor scores are probably incorrect.  Try a different factor score estimation method.

pdf(paste0(scree_path, "ma_rr_33_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_33_poly_minres <- fa.parallel(dat_ma_rr_only_33, correct = 0, cor = "poly"))
dev.off()

result_33_poly_minres$ncomp == 3

sum(result_33_poly_minres$pc.values > result_33_poly_minres$pc.sim)  == 4
sum(result_33_poly_minres$pc.values > result_33_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on 33 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_33_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_33, nfactors = 2:5, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_33_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_33, nfactors = 2:5, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_33_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_33, nfactors = 2:5, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_33_oblimin_wlsmv_fit, red_33_items_path, "ma_efa_33_oblimin_wlsmv")
export_efa_res(ma_efa_33_geomin_wlsmv_fit,  red_33_items_path, "ma_efa_33_geomin_wlsmv")
export_efa_res(ma_efa_33_promax_wlsmv_fit,  red_33_items_path, "ma_efa_33_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 32 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding positive threat item that consistently has weak loadings on
# positive threat factor and loads onto nonthreat factor

dat_ma_rr_only_32 <- dat_ma_rr_only_33[, names(dat_ma_rr_only_33)[names(dat_ma_rr_only_33) != 
                                                                    "pos_thr_job_3a"]]

length(dat_ma_rr_only_32) == 32

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_32, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_32_items_path <- paste0(efa_path, "red_32_items/")
scree_path        <- paste0(red_32_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_32_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_32, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_32_scree_pa_ml.pdf"), height = 6, width = 6)
(result_32_ml <- fa.parallel(dat_ma_rr_only_32, fm = "ml", n.iter = 100))
dev.off()

result_32_ml$ncomp == 3

sum(result_32_ml$pc.values > result_32_ml$pc.sim)  == 3
sum(result_32_ml$pc.values > result_32_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 4 components

pdf(paste0(scree_path, "ma_rr_32_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_32_poly_ml <- fa.parallel(dat_ma_rr_only_32, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_32_poly_ml$ncomp == 3

sum(result_32_poly_ml$pc.values > result_32_poly_ml$pc.sim)  == 4 # Based on simulated data
sum(result_32_poly_ml$pc.values > result_32_poly_ml$pc.simr) == 3 # Based on resampled data

# Warnings (one kind) for "minres":
# In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   The estimated weights for the factor scores are probably incorrect.  Try a different factor score estimation method.

pdf(paste0(scree_path, "ma_rr_32_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_32_poly_minres <- fa.parallel(dat_ma_rr_only_32, correct = 0, cor = "poly"))
dev.off()

result_32_poly_minres$ncomp == 3

sum(result_32_poly_minres$pc.values > result_32_poly_minres$pc.sim)  == 3
sum(result_32_poly_minres$pc.values > result_32_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on 32 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_32_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_32, nfactors = 2:5, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_32_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_32, nfactors = 2:5, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_32_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_32, nfactors = 2:5, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_32_oblimin_wlsmv_fit, red_32_items_path, "ma_efa_32_oblimin_wlsmv")
export_efa_res(ma_efa_32_geomin_wlsmv_fit,  red_32_items_path, "ma_efa_32_geomin_wlsmv")
export_efa_res(ma_efa_32_promax_wlsmv_fit,  red_32_items_path, "ma_efa_32_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 31 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding positive threat item that consistently has weak loadings on
# positive threat factor and loads onto nonthreat factor

dat_ma_rr_only_31 <- dat_ma_rr_only_32[, names(dat_ma_rr_only_32)[names(dat_ma_rr_only_32) != 
                                                                    "pos_thr_meeting_friend_5c"]]

length(dat_ma_rr_only_31) == 31

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_31, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_31_items_path <- paste0(efa_path, "red_31_items/")
scree_path        <- paste0(red_31_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_31_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_31, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_31_scree_pa_ml.pdf"), height = 6, width = 6)
(result_31_ml <- fa.parallel(dat_ma_rr_only_31, fm = "ml", n.iter = 100))
dev.off()

result_31_ml$ncomp == 3

sum(result_31_ml$pc.values > result_31_ml$pc.sim)  == 3
sum(result_31_ml$pc.values > result_31_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 4 components

pdf(paste0(scree_path, "ma_rr_31_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_31_poly_ml <- fa.parallel(dat_ma_rr_only_31, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_31_poly_ml$ncomp == 3

sum(result_31_poly_ml$pc.values > result_31_poly_ml$pc.sim)  == 4 # Based on simulated data
sum(result_31_poly_ml$pc.values > result_31_poly_ml$pc.simr) == 3 # Based on resampled data

# Warnings (one kind) for "minres":
# In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   The estimated weights for the factor scores are probably incorrect.  Try a different factor score estimation method.

pdf(paste0(scree_path, "ma_rr_31_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_31_poly_minres <- fa.parallel(dat_ma_rr_only_31, correct = 0, cor = "poly"))
dev.off()

result_31_poly_minres$ncomp == 3

sum(result_31_poly_minres$pc.values > result_31_poly_minres$pc.sim)  == 4
sum(result_31_poly_minres$pc.values > result_31_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on 31 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_31_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_31, nfactors = 2:5, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_31_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_31, nfactors = 2:5, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_31_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_31, nfactors = 2:5, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_31_oblimin_wlsmv_fit, red_31_items_path, "ma_efa_31_oblimin_wlsmv")
export_efa_res(ma_efa_31_geomin_wlsmv_fit,  red_31_items_path, "ma_efa_31_geomin_wlsmv")
export_efa_res(ma_efa_31_promax_wlsmv_fit,  red_31_items_path, "ma_efa_31_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 30 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding positive nonthreat item that consistently cross-loads onto
# positive threat factor

dat_ma_rr_only_30 <- dat_ma_rr_only_31[, names(dat_ma_rr_only_31)[names(dat_ma_rr_only_31) != 
                                                                    "pos_non_scrape_7d"]]

length(dat_ma_rr_only_30) == 30

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_30, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_30_items_path <- paste0(efa_path, "red_30_items/")
scree_path        <- paste0(red_30_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_30_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_30, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_30_scree_pa_ml.pdf"), height = 6, width = 6)
(result_30_ml <- fa.parallel(dat_ma_rr_only_30, fm = "ml", n.iter = 100))
dev.off()

result_30_ml$ncomp == 3

sum(result_30_ml$pc.values > result_30_ml$pc.sim)  == 3
sum(result_30_ml$pc.values > result_30_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 4 components

pdf(paste0(scree_path, "ma_rr_30_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_30_poly_ml <- fa.parallel(dat_ma_rr_only_30, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_30_poly_ml$ncomp == 3

sum(result_30_poly_ml$pc.values > result_30_poly_ml$pc.sim)  == 4 # Based on simulated data
sum(result_30_poly_ml$pc.values > result_30_poly_ml$pc.simr) == 3 # Based on resampled data

# Warnings (two kinds) for "minres":
# In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   The estimated weights for the factor scores are probably incorrect.  Try a different factor score estimation method.
# In fac(r = r, nfactors = nfactors, n.obs = n.obs, rotate = rotate,  :
#   An ultra-Heywood case was detected.  Examine the results carefully

pdf(paste0(scree_path, "ma_rr_30_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_30_poly_minres <- fa.parallel(dat_ma_rr_only_30, correct = 0, cor = "poly"))
dev.off()

result_30_poly_minres$ncomp == 3

sum(result_30_poly_minres$pc.values > result_30_poly_minres$pc.sim)  == 4
sum(result_30_poly_minres$pc.values > result_30_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on 30 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# Note: Standard errors of 4-factor solution for oblimin rotation seem implausibly
# large (does not resolve when changing seed to 12345, but does resolve when using
# listwise deletion below). Otherwise no Heywood cases.

set.seed(1234)
ma_efa_30_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_30, nfactors = 2:5, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_30_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_30, nfactors = 2:5, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_30_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_30, nfactors = 2:5, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# No Heywood cases for 4-factor oblimin solution using listwise deletion, but fit statistics 
# and magnitudes of factor loadings are very similar to those using pairwise deletion

set.seed(1234)
ma_efa_30_oblimin_wlsmv_fit_nf4_lw <- efa(data = dat_ma_rr_only_30, nfactors = 4, rotation = "oblimin", 
                                          ordered = TRUE, estimator = "WLSMV", missing = "listwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_30_oblimin_wlsmv_fit,        red_30_items_path, "ma_efa_30_oblimin_wlsmv")
export_efa_res(ma_efa_30_geomin_wlsmv_fit,         red_30_items_path, "ma_efa_30_geomin_wlsmv")
export_efa_res(ma_efa_30_promax_wlsmv_fit,         red_30_items_path, "ma_efa_30_promax_wlsmv")

export_efa_res(ma_efa_30_oblimin_wlsmv_fit_nf4_lw, red_30_items_path, "ma_efa_30_oblimin_wlsmv_nf4_lw")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 29 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding positive nonthreat item that consistently cross-loads onto
# positive threat factor

dat_ma_rr_only_29 <- dat_ma_rr_only_30[, names(dat_ma_rr_only_30)[names(dat_ma_rr_only_30) != 
                                                                    "pos_non_lunch_6d"]]

length(dat_ma_rr_only_29) == 29

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_29, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_29_items_path <- paste0(efa_path, "red_29_items/")
scree_path        <- paste0(red_29_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_29_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_29, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_29_scree_pa_ml.pdf"), height = 6, width = 6)
(result_29_ml <- fa.parallel(dat_ma_rr_only_29, fm = "ml", n.iter = 100))
dev.off()

result_29_ml$ncomp == 3

sum(result_29_ml$pc.values > result_29_ml$pc.sim)  == 3
sum(result_29_ml$pc.values > result_29_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_29_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_29_poly_ml <- fa.parallel(dat_ma_rr_only_29, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_29_poly_ml$ncomp == 3

sum(result_29_poly_ml$pc.values > result_29_poly_ml$pc.sim)  == 3 # Based on simulated data
sum(result_29_poly_ml$pc.values > result_29_poly_ml$pc.simr) == 3 # Based on resampled data

# No warnings

pdf(paste0(scree_path, "ma_rr_29_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_29_poly_minres <- fa.parallel(dat_ma_rr_only_29, correct = 0, cor = "poly"))
dev.off()

result_29_poly_minres$ncomp == 3

sum(result_29_poly_minres$pc.values > result_29_poly_minres$pc.sim)  == 3
sum(result_29_poly_minres$pc.values > result_29_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 3 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, or 4.

# ---------------------------------------------------------------------------- #
# Run EFA based on 29 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# Note: Standard errors of 2-factor solution for geomin rotation seem implausibly
# large (does not resolve when changing seed to 12345, but does resolve when using
# listwise deletion below). Otherwise no Heywood cases.

set.seed(1234)
ma_efa_29_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_29, nfactors = 2:4, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_29_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_29, nfactors = 2:4, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_29_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_29, nfactors = 2:4, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# No Heywood cases for 2-factor geomin solution using listwise deletion, but fit statistics 
# and magnitudes of factor loadings are very similar to those using pairwise deletion

set.seed(1234)
ma_efa_29_geomin_wlsmv_fit_nf2_lw  <- efa(data = dat_ma_rr_only_29, nfactors = 2, rotation = "geomin", 
                                          ordered = TRUE, estimator = "WLSMV", missing = "listwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_29_oblimin_wlsmv_fit,       red_29_items_path, "ma_efa_29_oblimin_wlsmv")
export_efa_res(ma_efa_29_geomin_wlsmv_fit,        red_29_items_path, "ma_efa_29_geomin_wlsmv")
export_efa_res(ma_efa_29_promax_wlsmv_fit,        red_29_items_path, "ma_efa_29_promax_wlsmv")

export_efa_res(ma_efa_29_geomin_wlsmv_fit_nf2_lw, red_29_items_path, "ma_efa_29_geomin_wlsmv_nf2_lw")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on 28 reduced items ----
# ---------------------------------------------------------------------------- #

# Consider excluding positive nonthreat item that consistently cross-loads onto
# positive threat factor

dat_ma_rr_only_28 <- dat_ma_rr_only_29[, names(dat_ma_rr_only_29)[names(dat_ma_rr_only_29) != 
                                                                    "pos_non_blood_test_9a"]]

length(dat_ma_rr_only_28) == 28

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_only_28, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

red_28_items_path <- paste0(efa_path, "red_28_items/")
scree_path        <- paste0(red_28_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_28_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_only_28, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_28_scree_pa_ml.pdf"), height = 6, width = 6)
(result_28_ml <- fa.parallel(dat_ma_rr_only_28, fm = "ml", n.iter = 100))
dev.off()

result_28_ml$ncomp == 3

sum(result_28_ml$pc.values > result_28_ml$pc.sim)  == 3
sum(result_28_ml$pc.values > result_28_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_28_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_28_poly_ml <- fa.parallel(dat_ma_rr_only_28, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_28_poly_ml$ncomp == 3

sum(result_28_poly_ml$pc.values > result_28_poly_ml$pc.sim)  == 3 # Based on simulated data
sum(result_28_poly_ml$pc.values > result_28_poly_ml$pc.simr) == 3 # Based on resampled data

# Warnings (two kinds) for "minres":
# In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   The estimated weights for the factor scores are probably incorrect.  Try a different factor score estimation method.
# In fac(r = r, nfactors = nfactors, n.obs = n.obs, rotate = rotate,  :
#   An ultra-Heywood case was detected.  Examine the results carefully

pdf(paste0(scree_path, "ma_rr_28_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_28_poly_minres <- fa.parallel(dat_ma_rr_only_28, correct = 0, cor = "poly"))
dev.off()

result_28_poly_minres$ncomp == 3

sum(result_28_poly_minres$pc.values > result_28_poly_minres$pc.sim)  == 3
sum(result_28_poly_minres$pc.values > result_28_poly_minres$pc.simr) == 3

# Based on scree plot, considered 4 factors, and PA-PCA-m suggests up to 3 (though 
# looking at +/- 1 is recommended). Thus, consider 2, 3, or 4.

# ---------------------------------------------------------------------------- #
# Run EFA based on 28 items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_28_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_only_28, nfactors = 2:4, rotation = "oblimin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_28_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_only_28, nfactors = 2:4, rotation = "geomin", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_28_promax_wlsmv_fit  <- efa(data = dat_ma_rr_only_28, nfactors = 2:4, rotation = "promax", 
                                   ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_28_oblimin_wlsmv_fit, red_28_items_path, "ma_efa_28_oblimin_wlsmv")
export_efa_res(ma_efa_28_geomin_wlsmv_fit,  red_28_items_path, "ma_efa_28_geomin_wlsmv")
export_efa_res(ma_efa_28_promax_wlsmv_fit,  red_28_items_path, "ma_efa_28_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Determine number of factors based on only threat items ----
# ---------------------------------------------------------------------------- #

# After not obtaining good absolute fit for revised CFA models based on 28 items 
# (see "run_revised_cfa_ma.R"), try focusing on only threat items. Start over with
# all items so that all of a given scenario's items can be removed at once.

dat_ma_rr_thr_only <- dat_ma_rr_only[names(dat_ma_rr_only) %in%
                                       c(dat_items$rr_pos_thr_items, 
                                         dat_items$rr_neg_thr_items)]

length(dat_ma_rr_thr_only) == 18

# Obtain eigenvalues of correlation matrix

eigen(cor(dat_ma_rr_thr_only, use = "pairwise.complete.obs"))$values

# Plot eigenvalues as scree plot to help decide how many factors to retain, which
# suggests a break point between cliff and scree at 4 factors

thr_items_path <- paste0(efa_path, "thr_items/")
scree_path     <- paste0(thr_items_path, "scree/")

dir.create(scree_path, recursive = TRUE)

pdf(paste0(scree_path, "ma_rr_thr_scree.pdf"), height = 6, width = 6)
scree(dat_ma_rr_thr_only, factors = FALSE, pc = TRUE)
dev.off()

# Also consider PA-PCA-m, which suggests an upper bound of 3 components

pdf(paste0(scree_path, "ma_rr_thr_scree_pa_ml.pdf"), height = 6, width = 6)
(result_thr_ml <- fa.parallel(dat_ma_rr_thr_only, fm = "ml", n.iter = 100))
dev.off()

result_thr_ml$ncomp == 3

sum(result_thr_ml$pc.values > result_thr_ml$pc.sim)  == 3
sum(result_thr_ml$pc.values > result_thr_ml$pc.simr) == 3

# PA-PCA-m on polychoric correlations suggests an upper bound of 4 components

pdf(paste0(scree_path, "ma_rr_thr_scree_pa_poly_ml.pdf"), height = 6, width = 6)
(result_thr_poly_ml <- fa.parallel(dat_ma_rr_thr_only, fm = "ml", n.iter = 10, correct = 0, cor = "poly"))
dev.off()

result_thr_poly_ml$ncomp == 3

sum(result_thr_poly_ml$pc.values > result_thr_poly_ml$pc.sim)  == 4 # Based on simulated data
sum(result_thr_poly_ml$pc.values > result_thr_poly_ml$pc.simr) == 3 # Based on resampled data

pdf(paste0(scree_path, "ma_rr_thr_scree_pa_poly_minres.pdf"), height = 6, width = 6)
(result_thr_poly_minres <- fa.parallel(dat_ma_rr_thr_only, correct = 0, cor = "poly"))
dev.off()

result_thr_poly_minres$ncomp == 3

sum(result_thr_poly_minres$pc.values > result_thr_poly_minres$pc.sim)  == 4
sum(result_thr_poly_minres$pc.values > result_thr_poly_minres$pc.simr) == 3

# Based on scree plot, considered 5 factors, and PA-PCA-m suggests up to 4 (though 
# looking at +/- 1 is recommended). Thus, consider 1, 2, 3, 4, or 5.

# ---------------------------------------------------------------------------- #
# Run EFA based on threat items using WLSMV ----
# ---------------------------------------------------------------------------- #

# No Heywood cases

set.seed(1234)
ma_efa_thr_oblimin_wlsmv_fit <- efa(data = dat_ma_rr_thr_only, nfactors = 1:5, rotation = "oblimin", 
                                    ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_thr_geomin_wlsmv_fit  <- efa(data = dat_ma_rr_thr_only, nfactors = 1:5, rotation = "geomin", 
                                    ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
set.seed(1234)
ma_efa_thr_promax_wlsmv_fit  <- efa(data = dat_ma_rr_thr_only, nfactors = 1:5, rotation = "promax", 
                                    ordered = TRUE, estimator = "WLSMV", missing = "pairwise")

# Export basic results and details to TXT and loadings to CSV

export_efa_res(ma_efa_thr_oblimin_wlsmv_fit, thr_items_path, "ma_efa_thr_oblimin_wlsmv")
export_efa_res(ma_efa_thr_geomin_wlsmv_fit,  thr_items_path, "ma_efa_thr_geomin_wlsmv")
export_efa_res(ma_efa_thr_promax_wlsmv_fit,  thr_items_path, "ma_efa_thr_promax_wlsmv")

# ---------------------------------------------------------------------------- #
# Export results ----
# ---------------------------------------------------------------------------- #

# Collect parallel analysis and EFA results in list

res_ma_pa <- list(ma_pa_ml               = result_ml,
                  ma_pa_poly_ml          = result_poly_ml,
                  ma_pa_poly_minres      = result_poly_minres,
                  ma_pa_35_ml            = result_35_ml,
                  ma_pa_35_poly_ml       = result_35_poly_ml,
                  ma_pa_35_poly_minres   = result_35_poly_minres,
                  ma_pa_34_ml            = result_34_ml,
                  ma_pa_34_poly_ml       = result_34_poly_ml,
                  ma_pa_34_poly_minres   = result_34_poly_minres,
                  ma_pa_33_ml            = result_33_ml,
                  ma_pa_33_poly_ml       = result_33_poly_ml,
                  ma_pa_33_poly_minres   = result_33_poly_minres,
                  ma_pa_32_ml            = result_32_ml,
                  ma_pa_32_poly_ml       = result_32_poly_ml,
                  ma_pa_32_poly_minres   = result_32_poly_minres,
                  ma_pa_31_ml            = result_31_ml,
                  ma_pa_31_poly_ml       = result_31_poly_ml,
                  ma_pa_31_poly_minres   = result_31_poly_minres,
                  ma_pa_30_ml            = result_30_ml,
                  ma_pa_30_poly_ml       = result_30_poly_ml,
                  ma_pa_30_poly_minres   = result_30_poly_minres,
                  ma_pa_29_ml            = result_29_ml,
                  ma_pa_29_poly_ml       = result_29_poly_ml,
                  ma_pa_29_poly_minres   = result_29_poly_minres,
                  ma_pa_28_ml            = result_28_ml,
                  ma_pa_28_poly_ml       = result_28_poly_ml,
                  ma_pa_28_poly_minres   = result_28_poly_minres,
                  result_thr_ml          = result_thr_ml,
                  result_thr_poly_ml     = result_thr_poly_ml,
                  result_thr_poly_minres = result_thr_poly_minres)

res_ma_efa <- list(ma_efa_oblimin_pml_ac_fit          = ma_efa_oblimin_pml_ac_fit,
                   ma_efa_geomin_pml_ac_fit           = ma_efa_geomin_pml_ac_fit,
                   ma_efa_promax_pml_ac_fit           = ma_efa_promax_pml_ac_fit,
                   ma_efa_oblimin_fiml_fit            = ma_efa_oblimin_fiml_fit,
                   ma_efa_geomin_fiml_fit             = ma_efa_geomin_fiml_fit,
                   ma_efa_promax_fiml_fit             = ma_efa_promax_fiml_fit,
                   ma_efa_oblimin_wlsmv_fit           = ma_efa_oblimin_wlsmv_fit,
                   ma_efa_geomin_wlsmv_fit            = ma_efa_geomin_wlsmv_fit,
                   ma_efa_promax_wlsmv_fit            = ma_efa_promax_wlsmv_fit,
                   ma_efa_35_oblimin_wlsmv_fit        = ma_efa_35_oblimin_wlsmv_fit,
                   ma_efa_35_geomin_wlsmv_fit         = ma_efa_35_geomin_wlsmv_fit,
                   ma_efa_35_promax_wlsmv_fit         = ma_efa_35_promax_wlsmv_fit,
                   ma_efa_34_oblimin_wlsmv_fit        = ma_efa_34_oblimin_wlsmv_fit,
                   ma_efa_34_geomin_wlsmv_fit         = ma_efa_34_geomin_wlsmv_fit,
                   ma_efa_34_promax_wlsmv_fit         = ma_efa_34_promax_wlsmv_fit,
                   ma_efa_33_oblimin_wlsmv_fit        = ma_efa_33_oblimin_wlsmv_fit,
                   ma_efa_33_geomin_wlsmv_fit         = ma_efa_33_geomin_wlsmv_fit,
                   ma_efa_33_promax_wlsmv_fit         = ma_efa_33_promax_wlsmv_fit,
                   ma_efa_32_oblimin_wlsmv_fit        = ma_efa_32_oblimin_wlsmv_fit,
                   ma_efa_32_geomin_wlsmv_fit         = ma_efa_32_geomin_wlsmv_fit,
                   ma_efa_32_promax_wlsmv_fit         = ma_efa_32_promax_wlsmv_fit,
                   ma_efa_31_oblimin_wlsmv_fit        = ma_efa_31_oblimin_wlsmv_fit,
                   ma_efa_31_geomin_wlsmv_fit         = ma_efa_31_geomin_wlsmv_fit,
                   ma_efa_31_promax_wlsmv_fit         = ma_efa_31_promax_wlsmv_fit,
                   ma_efa_30_oblimin_wlsmv_fit        = ma_efa_30_oblimin_wlsmv_fit,
                   ma_efa_30_oblimin_wlsmv_fit_nf4_lw = ma_efa_30_oblimin_wlsmv_fit_nf4_lw,
                   ma_efa_30_geomin_wlsmv_fit         = ma_efa_30_geomin_wlsmv_fit,
                   ma_efa_30_promax_wlsmv_fit         = ma_efa_30_promax_wlsmv_fit,
                   ma_efa_29_oblimin_wlsmv_fit        = ma_efa_29_oblimin_wlsmv_fit,
                   ma_efa_29_geomin_wlsmv_fit         = ma_efa_29_geomin_wlsmv_fit,
                   ma_efa_29_geomin_wlsmv_fit_nf2_lw  = ma_efa_29_geomin_wlsmv_fit_nf2_lw,
                   ma_efa_29_promax_wlsmv_fit         = ma_efa_29_promax_wlsmv_fit,
                   ma_efa_28_oblimin_wlsmv_fit        = ma_efa_28_oblimin_wlsmv_fit,
                   ma_efa_28_geomin_wlsmv_fit         = ma_efa_28_geomin_wlsmv_fit,
                   ma_efa_28_promax_wlsmv_fit         = ma_efa_28_promax_wlsmv_fit,
                   ma_efa_thr_oblimin_wlsmv_fit       = ma_efa_thr_oblimin_wlsmv_fit,
                   ma_efa_thr_geomin_wlsmv_fit        = ma_efa_thr_geomin_wlsmv_fit,
                   ma_efa_thr_promax_wlsmv_fit        = ma_efa_thr_promax_wlsmv_fit)

# Save results

save(res_ma_pa,  file = paste0(efa_path, "res_ma_pa.RData"))
save(res_ma_efa, file = paste0(efa_path, "res_ma_efa.RData"))