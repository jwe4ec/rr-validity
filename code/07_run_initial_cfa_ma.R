# ---------------------------------------------------------------------------- #
# Run Initial CFA in Managing Anxiety Study
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

# Load packages

pkgs <- c("lavaan", "lavaanExtra", "semPlot", "qgraph")
groundhog.library(pkgs, groundhog_day)

# Allow printing of more lines

getOption("max.print")
options(max.print=10000)

# ---------------------------------------------------------------------------- #
# Import data and create CFA initial results path ----
# ---------------------------------------------------------------------------- #

load("./data/ma/intermediate/dat_ma_rr.RData")

load("./data/helper/dat_items.RData")
load("./data/helper/rr_item_map.RData")

cfa_initial_path <- "./results/cfa/initial/"
dir.create(cfa_initial_path)

# ---------------------------------------------------------------------------- #
# Restrict to RR items ----
# ---------------------------------------------------------------------------- #

dat_ma_rr_only <- dat_ma_rr[, names(dat_ma_rr) %in% dat_items$rr_items]

# ---------------------------------------------------------------------------- #
# Run CFA for correlated-factors model ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_pos_thr_items,
            pos_non = dat_items$rr_pos_non_items,
            neg_thr = dat_items$rr_neg_thr_items,
            neg_non = dat_items$rr_neg_non_items)

ma_cfa_cf <- write_lavaan(latent = lat)
cat(ma_cfa_cf)

# ---------------------------------------------------------------------------- #
# Run using PML-AC ----
# ---------------------------------------------------------------------------- #

# Note: Negative thresholds are OK

filename <- file(paste0(cfa_initial_path, "ma_cfa_cf_pml_ac.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_cf_pml_ac_fit <- cfa(ma_cfa_cf, data = dat_ma_rr_only, std.lv = TRUE,
                            ordered = TRUE, estimator = "PML", missing = "available.cases")
(ma_cfa_cf_pml_ac_summ <- summary(ma_cfa_cf_pml_ac_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_cf_pml_ac_std_all <- standardizedsolution(ma_cfa_cf_pml_ac_fit))
sink(type = "message")
sink()

# TODO: Unsure why scaled CFI is 0 (worst possible value), whereas TLI is 0.877 
# (closer to 1 is better fit) and SRMR is .09 (close to 0, showing good fit)





ma_cfa_cf_pml_ac_p <- 
  semPaths(ma_cfa_cf_pml_ac_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_non_items), rev(dat_items$rr_neg_thr_items),
                         rev(dat_items$rr_pos_non_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 30, 1, 12),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_cf_pml_ac_plot")

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_initial_path, "ma_cfa_cf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_cf_wlsmv_fit <- cfa(ma_cfa_cf, data = dat_ma_rr_only, std.lv = TRUE,
                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
(ma_cfa_cf_wlsmv_summ <- summary(ma_cfa_cf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_cf_wlsmv_std_all <- standardizedsolution(ma_cfa_cf_wlsmv_fit))
sink(type = "message")
sink()

ma_cfa_cf_wlsmv_p <- 
  semPaths(ma_cfa_cf_wlsmv_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_non_items), rev(dat_items$rr_neg_thr_items),
                         rev(dat_items$rr_pos_non_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 30, 1, 12),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_cf_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run using FIML ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_initial_path, "ma_cfa_cf_fiml.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_cf_fiml_fit <- cfa(ma_cfa_cf, data = dat_ma_rr_only, std.lv = TRUE,
                          estimator = "MLR", missing = "fiml")
(ma_cfa_cf_fiml_summ <- summary(ma_cfa_cf_fiml_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_cf_fiml_std_all <- standardizedsolution(ma_cfa_cf_fiml_fit))
sink(type = "message")
sink()

ma_cfa_cf_fiml_p <- 
  semPaths(ma_cfa_cf_fiml_fit, what = "path", whatLabels = "stand", 
           intercepts = FALSE, residuals = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr"),
           manifests = c(rev(dat_items$rr_neg_non_items), rev(dat_items$rr_neg_thr_items),
                         rev(dat_items$rr_pos_non_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 30, 1, 12),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_cf_fiml_plot")

# ---------------------------------------------------------------------------- #
# Run CFA for bifactor model ----
# ---------------------------------------------------------------------------- #

# Define model (note: in addition to examples cited in preregistration, Morin et al., 
# 2020, p. 1055, says the two general factors should be allowed to correlate)

lat <- list(pos     = dat_items$rr_pos_items,
            neg     = dat_items$rr_neg_items,
            pos_thr = dat_items$rr_pos_thr_items,
            pos_non = dat_items$rr_pos_non_items,
            neg_thr = dat_items$rr_neg_thr_items,
            neg_non = dat_items$rr_neg_non_items)

cov <- list(pos = "neg")

ma_cfa_bf <- write_lavaan(latent = lat, covariance = cov)
cat(ma_cfa_bf)

# ---------------------------------------------------------------------------- #
# Run using PML-AC ----
# ---------------------------------------------------------------------------- #

# Note: Negative thresholds are OK

filename <- file(paste0(cfa_initial_path, "ma_cfa_bf_pml_ac.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_bf_pml_ac_fit <- cfa(ma_cfa_bf, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                            ordered = TRUE, estimator = "PML", missing = "available.cases")
(ma_cfa_bf_pml_ac_summ <- summary(ma_cfa_bf_pml_ac_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_bf_pml_ac_std_all <- standardizedsolution(ma_cfa_bf_pml_ac_fit))
sink(type = "message")
sink()

# TODO: Scaled CFI of 1 (perfect fit) and scaled TLI of 1 (good fit) consistent 
# with SRMR of .09 (good fit), but if scaled CFI is exactly 1, seems implausible





ma_cfa_bf_pml_ac_p <- 
  semPaths(ma_cfa_bf_pml_ac_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = c("neg", "pos"),
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr", "neg", "pos"),
           manifests = c(rev(dat_items$rr_neg_items), rev(dat_items$rr_pos_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 3, 1, 3),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_bf_pml_ac_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_bf_pml_ac_p$graphAttributes$Edges$edgeConnectPoints[1:36, 2]  <- 1.5*pi
ma_cfa_bf_pml_ac_p$graphAttributes$Edges$edgeConnectPoints[37:72, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_bf_pml_ac_p$graphAttributes$Edges$label.color[ma_cfa_bf_pml_ac_p$Arguments$edge.color == 
                                                       "transparent"] <- "transparent"

# Make edges black

ma_cfa_bf_pml_ac_p$graphAttributes$Edges$color[ma_cfa_bf_pml_ac_p$graphAttributes$Edges$color == 
                                                 "#808080FF"] <- "black"

plot(ma_cfa_bf_pml_ac_p)

# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_initial_path, "ma_cfa_bf_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_bf_wlsmv_fit <- cfa(ma_cfa_bf, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
(ma_cfa_bf_wlsmv_summ <- summary(ma_cfa_bf_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_bf_wlsmv_std_all <- standardizedsolution(ma_cfa_bf_wlsmv_fit))
sink(type = "message")
sink()

ma_cfa_bf_wlsmv_p <- 
  semPaths(ma_cfa_bf_wlsmv_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = c("neg", "pos"),
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr", "neg", "pos"),
           manifests = c(rev(dat_items$rr_neg_items), rev(dat_items$rr_pos_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 3, 1, 3),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_bf_wlsmv_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_bf_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[1:36, 2]  <- 1.5*pi
ma_cfa_bf_wlsmv_p$graphAttributes$Edges$edgeConnectPoints[37:72, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_bf_wlsmv_p$graphAttributes$Edges$label.color[ma_cfa_bf_wlsmv_p$Arguments$edge.color == 
                                                      "transparent"] <- "transparent"

# Make edges black

ma_cfa_bf_wlsmv_p$graphAttributes$Edges$color[ma_cfa_bf_wlsmv_p$graphAttributes$Edges$color == 
                                                "#808080FF"] <- "black"

plot(ma_cfa_bf_wlsmv_p)

# ---------------------------------------------------------------------------- #
# Run using FIML ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_initial_path, "ma_cfa_bf_fiml.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_bf_fiml_fit <- cfa(ma_cfa_bf, data = dat_ma_rr_only, std.lv = TRUE, orthogonal = TRUE,
                          estimator = "MLR", missing = "fiml")
(ma_cfa_bf_fiml_summ <- summary(ma_cfa_bf_fiml_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_bf_fiml_std_all <- standardizedsolution(ma_cfa_bf_fiml_fit))
sink(type = "message")
sink()

ma_cfa_bf_fiml_p <- 
  semPaths(ma_cfa_bf_fiml_fit, what = "path", whatLabels = "stand", layout = "tree2",
           intercepts = FALSE, residuals = FALSE,
           fixedStyle = c("transparent", "blank"),
           rotation = 2, curve = 2, curvature = 2.5, reorder = FALSE,
           bifactor = c("neg", "pos"),
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr", "neg", "pos"),
           manifests = c(rev(dat_items$rr_neg_items), rev(dat_items$rr_pos_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 3, 1, 3),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_bf_fiml_plot")

# Connect edges to left and right sides of nodes rather than to their centers

ma_cfa_bf_fiml_p$graphAttributes$Edges$edgeConnectPoints[1:36, 2]  <- 1.5*pi
ma_cfa_bf_fiml_p$graphAttributes$Edges$edgeConnectPoints[37:72, 2] <- 0.5*pi

# Hide labels for fixed edges (whose edge color was made transparent above)

ma_cfa_bf_fiml_p$graphAttributes$Edges$label.color[ma_cfa_bf_fiml_p$Arguments$edge.color == 
                                                     "transparent"] <- "transparent"

# Make edges black

ma_cfa_bf_fiml_p$graphAttributes$Edges$color[ma_cfa_bf_fiml_p$graphAttributes$Edges$color == 
                                               "#808080FF"] <- "black"

plot(ma_cfa_bf_fiml_p)

# ---------------------------------------------------------------------------- #
# Run CFA for higher-order model ----
# ---------------------------------------------------------------------------- #

# Define model

lat <- list(pos_thr = dat_items$rr_pos_thr_items,
            pos_non = dat_items$rr_pos_non_items,
            neg_thr = dat_items$rr_neg_thr_items,
            neg_non = dat_items$rr_neg_non_items,
            pos     = c("pos_thr", "pos_non"),
            neg     = c("neg_thr", "neg_non"))

ma_cfa_ho <- write_lavaan(latent = lat)
cat(ma_cfa_ho)

# ---------------------------------------------------------------------------- #
# Run using PML-AC ----
# ---------------------------------------------------------------------------- #

# Note: Negative thresholds are OK

filename <- file(paste0(cfa_initial_path, "ma_cfa_ho_pml_ac.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_ho_pml_ac_fit <- cfa(ma_cfa_ho, data = dat_ma_rr_only, std.lv = TRUE,
                            ordered = TRUE, estimator = "PML", missing = "available.cases",
                            check.vcov = FALSE)
(ma_cfa_ho_pml_ac_summ <- summary(ma_cfa_ho_pml_ac_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_ho_pml_ac_std_all <- standardizedsolution(ma_cfa_ho_pml_ac_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below (and negative factor covariance)





# Warning message:
#   In lavaan::lavaan(model = ma_cfa_ho, data = dat_ma_rr_only, ordered = TRUE,  :
#     lavaan WARNING:
#       the optimizer warns that a solution has NOT been found!
# lavaan 0.6.17 did NOT end normally after 4 iterations
# ** WARNING ** Estimates below are most likely unreliable

# Warning message:
#   In lav_object_summary(object = object, header = header, fit.measures = fit.measures,  :
#     lavaan WARNING: fit measures not available if model did not converge

# Warning message:
#   In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = object@SampleStats,  :
#     lavaan WARNING:
#       The variance-covariance matrix of the estimated parameters (vcov)
#       does not appear to be positive definite! The smallest eigenvalue
#       (= -4.248977e-32) is smaller than zero. This may be a symptom that
#       the model is not identified.

# However, no eigenvalues of "vcov" below are actually below 0, and -4.248977e-32 in warning 
# above is just barely negative and may simply be due to a machine precision issue (see 
# https://groups.google.com/g/lavaan/c/4y5pmqRz4nk/m/PXSq9VEdBwAJ). Tried adding "check.vcov = 
# FALSE" (see https://search.r-project.org/CRAN/refmans/lavaan/html/lavOptions.html), which
# removes the warning about nonpositive definite "vcov", but the other two warnings remain.

raw_cov                       <- cov(dat_ma_rr_only, use = "pairwise.complete.obs")
ma_cfa_ho_pml_ac_fit_samp_cov <- lavInspect(ma_cfa_ho_pml_ac_fit, "sampstat")$cov
ma_cfa_ho_pml_ac_fit_vcov     <- lavInspect(ma_cfa_ho_pml_ac_fit, "vcov")
ma_cfa_ho_pml_ac_fit_cov_ov   <- lavInspect(ma_cfa_ho_pml_ac_fit, "cov.ov")
ma_cfa_ho_pml_ac_fit_cov_lv   <- lavInspect(ma_cfa_ho_pml_ac_fit, "cov.lv")

eigen(raw_cov)$values
eigen(ma_cfa_ho_pml_ac_fit_samp_cov)$values
eigen(ma_cfa_ho_pml_ac_fit_vcov)$values
all(eigen(ma_cfa_ho_pml_ac_fit_vcov)$values > 0) # TODO: In contrast to warning, all are > 0
eigen(ma_cfa_ho_pml_ac_fit_cov_ov)$values
eigen(ma_cfa_ho_pml_ac_fit_cov_lv)$values

# TODO: Create figure if model runs

ma_cfa_ho_pml_ac_p <- NULL





# ---------------------------------------------------------------------------- #
# Run using WLSMV with pairwise deletion ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_initial_path, "ma_cfa_ho_wlsmv.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_ho_wlsmv_fit <- cfa(ma_cfa_ho, data = dat_ma_rr_only, std.lv = TRUE,
                           ordered = TRUE, estimator = "WLSMV", missing = "pairwise")
(ma_cfa_ho_wlsmv_summ <- summary(ma_cfa_ho_wlsmv_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_ho_wlsmv_std_all <- standardizedsolution(ma_cfa_ho_wlsmv_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below





# Warning message:
# In lavaan::lavaan(model = ma_cfa_ho, data = dat_ma_rr_only, ordered = TRUE,  :
#   lavaan WARNING:
#     the optimizer (NLMINB) claimed the model converged, but not all
#     elements of the gradient are (near) zero; the optimizer may not
#     have found a local solution use check.gradient = FALSE to skip
#     this check.
# lavaan 0.6.17 did NOT end normally after 2010 iterations
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

# Create figure without parameter estimates just to show structure

ma_cfa_ho_wlsmv_p <- 
  semPaths(ma_cfa_ho_wlsmv_fit, what = "path", whatLabels = "no", 
           intercepts = FALSE, residuals = FALSE, thresholds = FALSE,
           rotation = 2, curve = 2, curvature = 2.5, cardinal = TRUE, reorder = FALSE,
           latents = c("neg_non", "neg_thr", "pos_non", "pos_thr", "neg", "pos"),
           manifests = c(rev(dat_items$rr_neg_non_items), rev(dat_items$rr_neg_thr_items),
                         rev(dat_items$rr_pos_non_items), rev(dat_items$rr_pos_thr_items)),
           sizeLat = 8, sizeMan = 20, sizeMan2 = 2,
           nCharNodes = 0, label.scale = FALSE, label.cex = 0.7,
           edge.label.position = .6, edge.label.margin = .01,
           edge.color="black", edge.label.color="black", 
           mar = c(1, 6, 1, 6),
           filetype = "pdf", filename = "./results/cfa/initial/ma_cfa_ho_wlsmv_plot")

# ---------------------------------------------------------------------------- #
# Run using FIML ----
# ---------------------------------------------------------------------------- #

filename <- file(paste0(cfa_initial_path, "ma_cfa_ho_fiml.txt"), open = "wt")
sink(file = filename, type = "output")
sink(file = filename, type = "message")
set.seed(1234)
ma_cfa_ho_fiml_fit <- cfa(ma_cfa_ho, data = dat_ma_rr_only, std.lv = TRUE,
                          estimator = "MLR", missing = "fiml")
(ma_cfa_ho_fiml_summ <- summary(ma_cfa_ho_fiml_fit, fit.measures = TRUE, standardized = TRUE))
(ma_cfa_ho_fiml_std_all <- standardizedsolution(ma_cfa_ho_fiml_fit))
sink(type = "message")
sink()

# TODO: Resolve warnings below (and some negative residual variances of first-order 
# factors). Also, paths from second-order to first-order factors seem very large.





# Warning message:
# In lavaan::lavaan(model = ma_cfa_ho, data = dat_ma_rr_only, std.lv = TRUE,  :
#   lavaan WARNING:
#     the optimizer warns that a solution has NOT been found!
# lavaan 0.6.17 did NOT end normally after 1004 iterations
# ** WARNING ** Estimates below are most likely unreliable

# Warning message:
# In lav_object_summary(object = object, header = header, fit.measures = fit.measures,  :
#   lavaan WARNING: fit measures not available if model did not converge

# Warning message:
# In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = object@SampleStats,  :
#   lavaan WARNING:
#     Could not compute standard errors! The information matrix could
#     not be inverted. This may be a symptom that the model is not identified.

raw_cov                     <- cov(dat_ma_rr_only, use = "pairwise.complete.obs")
ma_cfa_ho_fiml_fit_samp_cov <- lavInspect(ma_cfa_ho_fiml_fit, "sampstat")$cov
ma_cfa_ho_fiml_fit_cov_ov   <- lavInspect(ma_cfa_ho_fiml_fit, "cov.ov")
ma_cfa_ho_fiml_fit_cov_lv   <- lavInspect(ma_cfa_ho_fiml_fit, "cov.lv")

eigen(raw_cov)$values
eigen(ma_cfa_ho_fiml_fit_samp_cov)$values
eigen(ma_cfa_ho_fiml_fit_cov_ov)$values
eigen(ma_cfa_ho_fiml_fit_cov_lv)$values # TODO: Small eigenvalues approaching 0 (semidefinite?)

ma_cfa_ho_fiml_fit_information <- lavInspect(ma_cfa_ho_fiml_fit, "information")
eigen(ma_cfa_ho_fiml_fit_information)$values # TODO: Small eigenvalues approaching 0 (semidefinite?)
det(ma_cfa_ho_fiml_fit_information)   # TODO: Determinant is not 0 (i.e., matrix not singular)
solve(ma_cfa_ho_fiml_fit_information) # TODO: But error states system is computationally singular

# TODO: Create figure if model runs

ma_cfa_ho_fiml_p <- NULL





# ---------------------------------------------------------------------------- #
# Export results and data ----
# ---------------------------------------------------------------------------- #

# Collect results in list

res_ma_cfa_initial <- list(ma_cfa_cf                = ma_cfa_cf,
                           ma_cfa_cf_pml_ac_fit     = ma_cfa_cf_pml_ac_fit,
                           ma_cfa_cf_pml_ac_summ    = ma_cfa_cf_pml_ac_summ,
                           ma_cfa_cf_pml_ac_std_all = ma_cfa_cf_pml_ac_std_all,
                           ma_cfa_cf_pml_ac_p       = ma_cfa_cf_pml_ac_p,
                           ma_cfa_cf_wlsmv_fit      = ma_cfa_cf_wlsmv_fit,
                           ma_cfa_cf_wlsmv_summ     = ma_cfa_cf_wlsmv_summ,
                           ma_cfa_cf_wlsmv_std_all  = ma_cfa_cf_wlsmv_std_all,
                           ma_cfa_cf_wlsmv_p        = ma_cfa_cf_wlsmv_p,
                           ma_cfa_cf_fiml_fit       = ma_cfa_cf_fiml_fit,
                           ma_cfa_cf_fiml_summ      = ma_cfa_cf_fiml_summ,
                           ma_cfa_cf_fiml_std_all   = ma_cfa_cf_fiml_std_all,
                           ma_cfa_cf_fiml_p         = ma_cfa_cf_fiml_p,
                           ma_cfa_bf                = ma_cfa_bf,
                           ma_cfa_bf_pml_ac_fit     = ma_cfa_bf_pml_ac_fit,
                           ma_cfa_bf_pml_ac_summ    = ma_cfa_bf_pml_ac_summ,
                           ma_cfa_bf_pml_ac_std_all = ma_cfa_bf_pml_ac_std_all,
                           ma_cfa_bf_pml_ac_p       = ma_cfa_bf_pml_ac_p,
                           ma_cfa_bf_wlsmv_fit      = ma_cfa_bf_wlsmv_fit,
                           ma_cfa_bf_wlsmv_summ     = ma_cfa_bf_wlsmv_summ,
                           ma_cfa_bf_wlsmv_std_all  = ma_cfa_bf_wlsmv_std_all,
                           ma_cfa_bf_wlsmv_p        = ma_cfa_bf_wlsmv_p,
                           ma_cfa_bf_fiml_fit       = ma_cfa_bf_fiml_fit,
                           ma_cfa_bf_fiml_summ      = ma_cfa_bf_fiml_summ,
                           ma_cfa_bf_fiml_std_all   = ma_cfa_bf_fiml_std_all,
                           ma_cfa_bf_fiml_p         = ma_cfa_bf_fiml_p,
                           ma_cfa_ho                = ma_cfa_ho,
                           ma_cfa_ho_pml_ac_fit     = ma_cfa_ho_pml_ac_fit,
                           ma_cfa_ho_pml_ac_summ    = ma_cfa_ho_pml_ac_summ,
                           ma_cfa_ho_pml_ac_std_all = ma_cfa_ho_pml_ac_std_all,
                           ma_cfa_ho_pml_ac_p       = ma_cfa_ho_pml_ac_p,
                           ma_cfa_ho_wlsmv_fit      = ma_cfa_ho_wlsmv_fit,
                           ma_cfa_ho_wlsmv_summ     = ma_cfa_ho_wlsmv_summ,
                           ma_cfa_ho_wlsmv_std_all  = ma_cfa_ho_wlsmv_std_all,
                           ma_cfa_ho_wlsmv_p        = ma_cfa_ho_wlsmv_p,
                           ma_cfa_ho_fiml_fit       = ma_cfa_ho_fiml_fit,
                           ma_cfa_ho_fiml_summ      = ma_cfa_ho_fiml_summ,
                           ma_cfa_ho_fiml_std_all   = ma_cfa_ho_fiml_std_all,
                           ma_cfa_ho_fiml_p         = ma_cfa_ho_fiml_p)

# Save results and data

save(res_ma_cfa_initial, file = paste0(cfa_initial_path, "res_ma_cfa_initial.RData"))

dir.create("./data/ma/final")
save(dat_ma_rr_only, file = "./data/ma/final/dat_ma_rr_only.RData")