# ---------------------------------------------------------------------------- #
# Write Tables to MS Word
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Although the present script writes a demographics table for Calm Thinking, we 
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

# Load "officer" package properties

source("./code/01c_set_officer_properties.R")

# ---------------------------------------------------------------------------- #
# Import flextables ----
# ---------------------------------------------------------------------------- #

dem_path <- "./results/demographics/"

load(paste0(dem_path, "dem_tbl_ma_ft.RData"))
load(paste0(dem_path, "dem_tbl_ct_ft.RData"))

# ---------------------------------------------------------------------------- #
# Write tables to MS Word ----
# ---------------------------------------------------------------------------- #

# Note: "flextable" seems to have a bug in which blank page is at end of doc

tbls <- list(dem_tbl_ma_ft, dem_tbl_ct_ft)

tbl_orientations <- c("p", "p")
tbl_numbers      <- c("1", "2")

doc <- read_docx()
doc <- body_set_default_section(doc, psect_prop)

for (i in 1:length(tbls)) {
  doc <- body_add_fpar(doc, fpar(ftext(paste0("Table ", tbl_numbers[[i]]),
                                       prop = text_prop_bold)))
  doc <- body_add_par(doc, "")
  
  doc <- body_add_flextable(doc, tbls[[i]], align = "left")
  
  if (tbl_orientations[[i]] == "p") {
    doc <- body_end_block_section(doc, block_section(psect_prop))
  } else if (tbl_orientations[[i]] == "l") {
    doc <- body_end_block_section(doc, block_section(lsect_prop))
  }
}

dem_tbl_path <- paste0(dem_path, "tables/")

dir.create(dem_tbl_path)

print(doc, target = paste0(dem_tbl_path, "dem_tbls.docx"))