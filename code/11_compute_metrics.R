# ---------------------------------------------------------------------------- #
# Compute Metrics in Managing Anxiety Study
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

pkgs <- c("lavaan", "semTools", "MBESS")
groundhog.library(pkgs, groundhog_day)

# Set seed

set.seed(1234)

# ---------------------------------------------------------------------------- #
# Define functions used in script ----
# ---------------------------------------------------------------------------- #

# Define function to export internal consistency results

export_omega_res <- function(reliability, omega, path, filename) {
  sink(paste0(path, paste0(filename, ".txt")))
  cat("Output from reliability():", "\n", "\n")
  print(reliability)
  cat("\n")
  cat("Output from compRelSEM():", "\n", "\n")
  print(omega)
  sink()
}

# Define function to export factor determinacy results

export_fd_res <- function(fd, path, filename) {
  sink(paste0(path, paste0(filename, ".txt")))
  cat("Correlation of factor scores with their respective factors:", "\n", "\n")
  print(fd$cor_fs_f)
  cat("\n")
  cat("Minimum possible correlation between two sets of factor scores:", "\n", "\n")
  print(fd$min_cor_fs)
  sink()
}

# Define function to export construct reliability results

export_H_res <- function(H, path, filename) {
  sink(paste0(path, paste0(filename, ".txt")))
  cat("Reliability of optimally weighted composites:", "\n", "\n")
  print(H)
  sink()
}

# ---------------------------------------------------------------------------- #
# Import data and results ----
# ---------------------------------------------------------------------------- #

load("./data/ma/final/dat_ma_rr_only.RData")

load("./data/helper/dat_items.RData")
load("./data/helper/rr_item_map.RData")

load("./results/cfa/initial/res_ma_cfa_initial.RData")
load("./results/cfa/revised/res_ma_cfa_28.RData")
load("./results/cfa/revised/res_ma_cfa_thr.RData")

# ---------------------------------------------------------------------------- #
# Restrict data to threat items ----
# ---------------------------------------------------------------------------- #

dat_ma_rr_thr_only <- dat_ma_rr_only[c(dat_items$rr_pos_thr_items, dat_items$rr_neg_thr_items)]

# ---------------------------------------------------------------------------- #
# Compute reliability of unit-weighted composites ----
# ---------------------------------------------------------------------------- #

# Note: "reliability" function is now deprecated. "compRelSEM" yields identical
# results as "omega2" (based on model-implied variance) from "reliability".

# For initial correlated-factor model

ma_cfa_cf_wlsmvv_reliability        <- reliability(res_ma_cfa_initial$ma_cfa_cf_wlsmv_fit)

ma_cfa_cf_wlsmvv_omega        <- compRelSEM(res_ma_cfa_initial$ma_cfa_cf_wlsmv_fit, 
                                            obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)

# For initial bifactor model

ma_cfa_bf_wlsmv_reliability         <- reliability(res_ma_cfa_initial$ma_cfa_bf_wlsmv_fit)

ma_cfa_bf_wlsmv_omega        <- compRelSEM(res_ma_cfa_initial$ma_cfa_bf_wlsmv_fit, 
                                            obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)

# For key 28-item correlated-factor models

ma_cfa_28_2cf_wlsmv_reliability        <- reliability(res_ma_cfa_28$ma_cfa_28_2cf_wlsmv_fit)
ma_cfa_28_3cf_wlsmv_reliability        <- reliability(res_ma_cfa_28$ma_cfa_28_3cf_wlsmv_fit)
ma_cfa_28_3cf_ce_thr_wlsmv_reliability <- reliability(res_ma_cfa_28$ma_cfa_28_3cf_ce_thr_wlsmv_fit)
ma_cfa_28_3cf_ce_all_wlsmv_reliability <- reliability(res_ma_cfa_28$ma_cfa_28_3cf_ce_all_wlsmv_fit)

ma_cfa_28_2cf_wlsmv_omega        <- compRelSEM(res_ma_cfa_28$ma_cfa_28_2cf_wlsmv_fit, 
                                               obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)
ma_cfa_28_3cf_wlsmv_omega        <- compRelSEM(res_ma_cfa_28$ma_cfa_28_3cf_wlsmv_fit, 
                                               obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)
ma_cfa_28_3cf_ce_thr_wlsmv_omega <- compRelSEM(res_ma_cfa_28$ma_cfa_28_3cf_ce_thr_wlsmv_fit, 
                                               obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)
ma_cfa_28_3cf_ce_all_wlsmv_omega <- compRelSEM(res_ma_cfa_28$ma_cfa_28_3cf_ce_all_wlsmv_fit, 
                                               obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)

# For key all threat-items correlated-factor models (note: for one-factor model, use 
# reverse-scored negative items to ensure positive factor loadings; Flora, 2020, p. 489)

ma_cfa_thr_2cf_wlsmv_reliability            <- reliability(res_ma_cfa_thr$ma_cfa_thr_2cf_wlsmv_fit)
ma_cfa_thr_2cf_ce_wlsmv_reliability         <- reliability(res_ma_cfa_thr$ma_cfa_thr_2cf_ce_wlsmv_fit)
ma_cfa_thr_1f_ce_scen_rev_wlsmv_reliability <- reliability(res_ma_cfa_thr$ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit)

ma_cfa_thr_2cf_wlsmv_omega            <- compRelSEM(res_ma_cfa_thr$ma_cfa_thr_2cf_wlsmv_fit, 
                                                    obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)
ma_cfa_thr_2cf_ce_wlsmv_omega         <- compRelSEM(res_ma_cfa_thr$ma_cfa_thr_2cf_ce_wlsmv_fit, 
                                                    obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)
ma_cfa_thr_1f_ce_scen_rev_wlsmv_omega <- compRelSEM(res_ma_cfa_thr$ma_cfa_thr_1f_ce_scen_rev_wlsmv_fit, 
                                                    obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)

# TODO: Consider computing for single-trait correlated-methods model and bifactor 
# model (note: need to reverse-score negative items first; see above)





# TODO: "ci.reliability" function is unable to account for error covariances. Thus,
# skip confidence intervals for now.





# Export results

internal_consistency_path <- "./results/internal_consistency/"

dir.create(internal_consistency_path, recursive = TRUE)

export_omega_res(ma_cfa_cf_wlsmvv_reliability,                ma_cfa_cf_wlsmvv_omega,
                 internal_consistency_path, "ma_cfa_cf_wlsmvv")
export_omega_res(ma_cfa_bf_wlsmv_reliability,                 ma_cfa_bf_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_bf_wlsmv")

export_omega_res(ma_cfa_28_2cf_wlsmv_reliability,             ma_cfa_28_2cf_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_28_2cf_wlsmv")
export_omega_res(ma_cfa_28_3cf_wlsmv_reliability,             ma_cfa_28_3cf_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_28_3cf_wlsmv")
export_omega_res(ma_cfa_28_3cf_ce_thr_wlsmv_reliability,      ma_cfa_28_3cf_ce_thr_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_28_3cf_ce_thr_wlsmv")
export_omega_res(ma_cfa_28_3cf_ce_all_wlsmv_reliability,      ma_cfa_28_3cf_ce_all_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_28_3cf_ce_all_wlsmv")

export_omega_res(ma_cfa_thr_2cf_wlsmv_reliability,            ma_cfa_thr_2cf_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_thr_2cf_wlsmv")
export_omega_res(ma_cfa_thr_2cf_ce_wlsmv_reliability,         ma_cfa_thr_2cf_ce_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_thr_2cf_ce_wlsmv")
export_omega_res(ma_cfa_thr_1f_ce_scen_rev_wlsmv_reliability, ma_cfa_thr_1f_ce_scen_rev_wlsmv_omega,
                 internal_consistency_path, "ma_cfa_thr_1f_ce_scen_rev_wlsmv")

# ---------------------------------------------------------------------------- #
# Compute factor determinacy ----
# ---------------------------------------------------------------------------- #

# Define function to compute correlation of factor scores with their respective
# factors (using Eq. 8 on p. 142 of Rodriguez et al., 2016; https://doi.org/gfrkj9)
# and minimum possible correlation between two sets of factor scores (see pp. 142
# and 150 of Rodriguez et al., 2016)

compute_factor_determinacy <- function(fit) {
  # Compute correlation of factor scores with their respective factors
  
  phi    <- lavInspect(fit, "cov.lv")
  lambda <- lavInspect(fit, "std")$lambda
  sigma  <- lavInspect(fit, "cor.ov")
  
  rho <- diag(phi %*% t(lambda) %*% solve(sigma) %*% lambda %*% phi) ^ .5
  
  # Compute minimum possible correlation between two sets of factor scores
  
  min_cor_fs <- 2 * rho ^ 2 - 1
  
  # Collect results in list
  
  fd = list(cor_fs_f   = rho,
            min_cor_fs = min_cor_fs)
  
  return(fd)
}

# Run function

ma_cfa_cf_wlsmv_fd            <- compute_factor_determinacy(res_ma_cfa_initial$ma_cfa_cf_wlsmv_fit)
ma_cfa_bf_wlsmv_fd            <- compute_factor_determinacy(res_ma_cfa_initial$ma_cfa_bf_wlsmv_fit)

ma_cfa_28_3cf_ce_thr_wlsmv_fd <- compute_factor_determinacy(res_ma_cfa_28$ma_cfa_28_3cf_ce_thr_wlsmv_fit)
ma_cfa_28_3cf_ce_all_wlsmv_fd <- compute_factor_determinacy(res_ma_cfa_28$ma_cfa_28_3cf_ce_all_wlsmv_fit)

ma_cfa_thr_2cf_ce_wlsmv_fd    <- compute_factor_determinacy(res_ma_cfa_thr$ma_cfa_thr_2cf_ce_wlsmv_fit)

# Export results

factor_determinacy_path <- "./results/factor_determinacy/"

dir.create(factor_determinacy_path, recursive = TRUE)

export_fd_res(ma_cfa_cf_wlsmv_fd,            factor_determinacy_path, "ma_cfa_cf_wlsmv_fd")
export_fd_res(ma_cfa_bf_wlsmv_fd,            factor_determinacy_path, "ma_cfa_bf_wlsmv_fd")

export_fd_res(ma_cfa_28_3cf_ce_thr_wlsmv_fd, factor_determinacy_path, "ma_cfa_28_3cf_ce_thr_wlsmv_fd")
export_fd_res(ma_cfa_28_3cf_ce_all_wlsmv_fd, factor_determinacy_path, "ma_cfa_28_3cf_ce_all_wlsmv_fd")

export_fd_res(ma_cfa_thr_2cf_ce_wlsmv_fd,    factor_determinacy_path, "ma_cfa_thr_2cf_ce_wlsmv_fd")

# ---------------------------------------------------------------------------- #
# Compute reliability of optimally weighted composites ----
# ---------------------------------------------------------------------------- #

# Define function to compute construct reliability for a given factor based on its
# completely standardized loadings (see Eq. 9 on p. 143 of Rodriguez et al., 2016)

compute_construct_reliability <- function(fit) {
  lambda <- lavInspect(fit, "std")$lambda
  
  lambda <- as.data.frame(lambda)
  
  H <- apply(lambda, 2, function(x) {
    1 / (1 + 1 / sum(x ^ 2 / (1 - x ^ 2)))
    })
}

# Run function

ma_cfa_cf_wlsmvv_H           <- compute_construct_reliability(res_ma_cfa_initial$ma_cfa_cf_wlsmv_fit)
ma_cfa_bf_wlsmvv_H           <- compute_construct_reliability(res_ma_cfa_initial$ma_cfa_bf_wlsmv_fit)

ma_cfa_28_3cf_ce_thr_wlsmv_H <- compute_construct_reliability(res_ma_cfa_28$ma_cfa_28_3cf_ce_thr_wlsmv_fit)
ma_cfa_28_3cf_ce_all_wlsmv_H <- compute_construct_reliability(res_ma_cfa_28$ma_cfa_28_3cf_ce_all_wlsmv_fit)

ma_cfa_thr_2cf_ce_wlsmv_H    <- compute_construct_reliability(res_ma_cfa_thr$ma_cfa_thr_2cf_ce_wlsmv_fit)

# Export results

construct_reliability_path <- "./results/construct_reliability/"

dir.create(construct_reliability_path, recursive = TRUE)

export_H_res(ma_cfa_cf_wlsmvv_H,           construct_reliability_path, "ma_cfa_cf_wlsmvv_H")
export_H_res(ma_cfa_bf_wlsmvv_H,           construct_reliability_path, "ma_cfa_bf_wlsmvv_H")

export_H_res(ma_cfa_28_3cf_ce_thr_wlsmv_H, construct_reliability_path, "ma_cfa_28_3cf_ce_thr_wlsmv_H")
export_H_res(ma_cfa_28_3cf_ce_all_wlsmv_H, construct_reliability_path, "ma_cfa_28_3cf_ce_all_wlsmv_H")

export_H_res(ma_cfa_thr_2cf_ce_wlsmv_H,    construct_reliability_path, "ma_cfa_thr_2cf_ce_wlsmv_H")