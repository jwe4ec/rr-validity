# ---------------------------------------------------------------------------- #
# Explore Uniquenesses
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
# Import EFA results ----
# ---------------------------------------------------------------------------- #

load("./results/efa/res_ma_efa.RData")

load("./data/helper/dat_items.RData")
load("./data/helper/rr_item_map.RData")

# ---------------------------------------------------------------------------- #
# Explore uniquenesses by scenario ----
# ---------------------------------------------------------------------------- #

# Define function to explore uniquenesses by scenario

explore_uniquenesses <- function(fit) {
  # Extract uniquenesses
  
  theta <- unlist(summary(fit)$efa$theta)
  
  theta_df <- data.frame(items_rename = names(theta),
                         uniqueness   = theta)
  
  row.names(theta_df) <- 1:nrow(theta_df)
  
  theta_df <- merge(theta_df, rr_item_map, "items_rename", all.x = TRUE, sort = FALSE)
  
  # Compute mean uniqueness of each scenario's items
  
  m_uniqueness_df <- aggregate(uniqueness ~ scenario, theta_df, mean)
  names(m_uniqueness_df)[names(m_uniqueness_df) == "uniqueness"] <- "m_uniqueness"
  
  theta_df <- merge(theta_df, m_uniqueness_df, "scenario", all.x = TRUE, sort = FALSE)
  
  # Sort items by mean uniqueness of their scenario
  
  theta_df <- theta_df[order(theta_df$m_uniqueness,
                             theta_df$valence, theta_df$threat_relevance,
                             decreasing = TRUE), ]
  
  # Visualize mean uniquenesses
  
  hist(m_uniqueness_df$m_uniqueness, breaks = seq(0, 1, .01))
  
  # Identify any scenarios with outlying mean uniquenesses > 2.5 median absolute 
  # deviations (per Leys et al., 2013, https://doi.org/f42jwd) above median
  
  med         <- median(unique(theta_df$m_uniqueness))
  med_abs_dev <- mad(unique(theta_df$m_uniqueness))
  
  upper_thres <- med + 2.5*med_abs_dev

  theta_df$mad_outlier <- NA
  theta_df$mad_outlier <- ifelse(theta_df$m_uniqueness > upper_thres, 1, 0)
  
  return(theta_df)
}

# Run function

unq_df_all_items_pml_ac_nf2   <- explore_uniquenesses(res_ma_efa$ma_efa_oblimin_pml_ac_fit$nf2)
unq_df_all_items_pml_ac_nf3   <- explore_uniquenesses(res_ma_efa$ma_efa_oblimin_pml_ac_fit$nf3)
unq_df_all_items_pml_ac_nf4   <- explore_uniquenesses(res_ma_efa$ma_efa_oblimin_pml_ac_fit$nf4)

unq_df_all_items_wlsmv_nf2    <- explore_uniquenesses(res_ma_efa$ma_efa_oblimin_wlsmv_fit$nf2)
unq_df_all_items_wlsmv_nf3    <- explore_uniquenesses(res_ma_efa$ma_efa_oblimin_wlsmv_fit$nf3)
unq_df_all_items_wlsmv_nf4    <- explore_uniquenesses(res_ma_efa$ma_efa_oblimin_wlsmv_fit$nf4)

unq_df_red_35_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_35_oblimin_wlsmv_fit$nf2)
unq_df_red_35_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_35_oblimin_wlsmv_fit$nf3)
unq_df_red_35_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_35_oblimin_wlsmv_fit$nf4)

unq_df_red_34_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_34_oblimin_wlsmv_fit$nf2)
unq_df_red_34_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_34_oblimin_wlsmv_fit$nf3)
unq_df_red_34_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_34_oblimin_wlsmv_fit$nf4)

unq_df_red_33_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_33_oblimin_wlsmv_fit$nf2)
unq_df_red_33_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_33_oblimin_wlsmv_fit$nf3)
unq_df_red_33_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_33_oblimin_wlsmv_fit$nf4)

unq_df_red_32_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_32_oblimin_wlsmv_fit$nf2)
unq_df_red_32_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_32_oblimin_wlsmv_fit$nf3)
unq_df_red_32_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_32_oblimin_wlsmv_fit$nf4)

unq_df_red_31_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_31_oblimin_wlsmv_fit$nf2)
unq_df_red_31_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_31_oblimin_wlsmv_fit$nf3)
unq_df_red_31_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_31_oblimin_wlsmv_fit$nf4)

unq_df_red_30_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_30_oblimin_wlsmv_fit$nf2)
unq_df_red_30_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_30_oblimin_wlsmv_fit$nf3)
unq_df_red_30_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_30_oblimin_wlsmv_fit$nf4)

unq_df_red_29_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_29_oblimin_wlsmv_fit$nf2)
unq_df_red_29_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_29_oblimin_wlsmv_fit$nf3)
unq_df_red_29_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_29_oblimin_wlsmv_fit$nf4)

unq_df_red_28_items_wlsmv_nf2 <- explore_uniquenesses(res_ma_efa$ma_efa_28_oblimin_wlsmv_fit$nf2)
unq_df_red_28_items_wlsmv_nf3 <- explore_uniquenesses(res_ma_efa$ma_efa_28_oblimin_wlsmv_fit$nf3)
unq_df_red_28_items_wlsmv_nf4 <- explore_uniquenesses(res_ma_efa$ma_efa_28_oblimin_wlsmv_fit$nf4)

unq_df_thr_items_wlsmv_nf1    <- explore_uniquenesses(res_ma_efa$ma_efa_thr_oblimin_wlsmv_fit$nf1)
unq_df_thr_items_wlsmv_nf2    <- explore_uniquenesses(res_ma_efa$ma_efa_thr_oblimin_wlsmv_fit$nf2)

# View any scenarios with outlying mean uniquenesses

  # Although "elevator" scenario is outlying for some factor solutions, across factor solutions and 
  # rotations only the item "pos_thr_elevator_1d" appears to have low loadings. Thus, try removing
  # this item but retaining the other "elevator" items.

unique(unq_df_all_items_pml_ac_nf2$scenario[unq_df_all_items_pml_ac_nf2$mad_outlier == 1])     == "elevator"
unique(unq_df_all_items_pml_ac_nf3$scenario[unq_df_all_items_pml_ac_nf3$mad_outlier == 1])     # None
unique(unq_df_all_items_pml_ac_nf4$scenario[unq_df_all_items_pml_ac_nf4$mad_outlier == 1])     # None

unique(unq_df_all_items_wlsmv_nf2$scenario[unq_df_all_items_wlsmv_nf2$mad_outlier == 1])       == "elevator"
unique(unq_df_all_items_wlsmv_nf3$scenario[unq_df_all_items_wlsmv_nf3$mad_outlier == 1])       # None
unique(unq_df_all_items_wlsmv_nf4$scenario[unq_df_all_items_wlsmv_nf4$mad_outlier == 1])       # None

unique(unq_df_red_35_items_wlsmv_nf2$scenario[unq_df_red_35_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_35_items_wlsmv_nf3$scenario[unq_df_red_35_items_wlsmv_nf3$mad_outlier == 1]) == "elevator"
unique(unq_df_red_35_items_wlsmv_nf4$scenario[unq_df_red_35_items_wlsmv_nf4$mad_outlier == 1]) == c("elevator", "shopping")

unique(unq_df_red_34_items_wlsmv_nf2$scenario[unq_df_red_34_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_34_items_wlsmv_nf3$scenario[unq_df_red_34_items_wlsmv_nf3$mad_outlier == 1]) == "elevator"
unique(unq_df_red_34_items_wlsmv_nf4$scenario[unq_df_red_34_items_wlsmv_nf4$mad_outlier == 1]) == "elevator"

  # Although "elevator" scenario remains outlying for 4-factor solution, its remaining 3 items
  # load strongly or fairly strongly on reasonable factors with no cross-loadings

unique(unq_df_red_33_items_wlsmv_nf2$scenario[unq_df_red_33_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_33_items_wlsmv_nf3$scenario[unq_df_red_33_items_wlsmv_nf3$mad_outlier == 1]) # None
unique(unq_df_red_33_items_wlsmv_nf4$scenario[unq_df_red_33_items_wlsmv_nf4$mad_outlier == 1]) == "elevator"

  # No outlying scenarios

unique(unq_df_red_32_items_wlsmv_nf2$scenario[unq_df_red_32_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_32_items_wlsmv_nf3$scenario[unq_df_red_32_items_wlsmv_nf3$mad_outlier == 1]) # None
unique(unq_df_red_32_items_wlsmv_nf4$scenario[unq_df_red_32_items_wlsmv_nf4$mad_outlier == 1]) # None

unique(unq_df_red_31_items_wlsmv_nf2$scenario[unq_df_red_31_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_31_items_wlsmv_nf3$scenario[unq_df_red_31_items_wlsmv_nf3$mad_outlier == 1]) # None
unique(unq_df_red_31_items_wlsmv_nf4$scenario[unq_df_red_31_items_wlsmv_nf4$mad_outlier == 1]) # None

unique(unq_df_red_30_items_wlsmv_nf2$scenario[unq_df_red_30_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_30_items_wlsmv_nf3$scenario[unq_df_red_30_items_wlsmv_nf3$mad_outlier == 1]) # None
unique(unq_df_red_30_items_wlsmv_nf4$scenario[unq_df_red_30_items_wlsmv_nf4$mad_outlier == 1]) # None

  # Although "elevator" scenario is now outlying for 3-factor solution, its remaining 3 items
  # still load strongly or fairly strongly on reasonable factors with no cross-loadings

unique(unq_df_red_29_items_wlsmv_nf2$scenario[unq_df_red_29_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_29_items_wlsmv_nf3$scenario[unq_df_red_29_items_wlsmv_nf3$mad_outlier == 1]) == "elevator"
unique(unq_df_red_29_items_wlsmv_nf4$scenario[unq_df_red_29_items_wlsmv_nf4$mad_outlier == 1]) # None

  # Although "elevator" and "shopping" scenarios are outlying for 3-factor solution, their (3
  # and 4 items, respectively) load strongly or fairly strongly on reasonable factors with no
  # consistent cross-loadings

unique(unq_df_red_28_items_wlsmv_nf2$scenario[unq_df_red_28_items_wlsmv_nf2$mad_outlier == 1]) # None
unique(unq_df_red_28_items_wlsmv_nf3$scenario[unq_df_red_28_items_wlsmv_nf3$mad_outlier == 1]) == c("elevator", "shopping")
unique(unq_df_red_28_items_wlsmv_nf4$scenario[unq_df_red_28_items_wlsmv_nf4$mad_outlier == 1]) # None

  # No outlying scenarios

unique(unq_df_thr_items_wlsmv_nf1$scenario[unq_df_thr_items_wlsmv_nf1$mad_outlier == 1]) # None
unique(unq_df_thr_items_wlsmv_nf2$scenario[unq_df_thr_items_wlsmv_nf2$mad_outlier == 1]) # None