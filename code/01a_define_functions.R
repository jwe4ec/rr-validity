# ---------------------------------------------------------------------------- #
# Define Functions
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Define version_control() ----
# ---------------------------------------------------------------------------- #

# Define function to check R version, load groundhog package, and return groundhog_day

version_control <- function() {
  # Ensure you are using the same version of R used at the time the script was 
  # written. To install a previous version, go to 
  # https://cran.r-project.org/bin/windows/base/old/
  
  script_R_version <- "R version 4.3.2 (2023-10-31 ucrt)"
  current_R_version <- R.Version()$version.string
  
  if(current_R_version != script_R_version) {
    warning(paste0("This script is based on ", script_R_version,
                   ". You are running ", current_R_version, "."))
  }
  
  # Load packages using "groundhog", which installs and loads the most recent
  # versions of packages available on the specified date ("groundhog_day"). This 
  # is important for reproducibility so that everyone running the script is using
  # the same versions of packages used at the time the script was written.
  
  # Note that packages may take longer to load the first time you load them with
  # "groundhog.library". This is because you may not have the correct versions of 
  # the packages installed based on the "groundhog_day". After "groundhog.library"
  # automatically installs the correct versions alongside other versions you may 
  # have installed, it will load the packages more quickly.
  
  # If in the process of loading packages with "groundhog.library" for the first 
  # time the console states that you first need to install "Rtools", follow steps 
  # here (https://cran.r-project.org/bin/windows/Rtools/) for installing "Rtools" 
  # and putting "Rtools" on the PATH. Then try loading the packages again.
  
  library(groundhog)
  meta.groundhog("2024-02-15")
  groundhog_day <- "2024-02-15"
  
  return(groundhog_day)
}

# ---------------------------------------------------------------------------- #
# Define identify_cols() ----
# ---------------------------------------------------------------------------- #

# Define function to identify columns matching a grep pattern in a data frame.
# Ignore case. Use lapply to apply the function to all data frames in a list.

identify_cols <- function(df, grep_pattern) {
  df_colnames <- colnames(df)
  
  selected_cols <- grep(grep_pattern, df_colnames, ignore.case = TRUE)
  if (length(selected_cols) != 0) {
    df_colnames[selected_cols]
  }
}

# ---------------------------------------------------------------------------- #
# Define check_relevant_files() ----
# ---------------------------------------------------------------------------- #

# Note: Although this function defines tables relevant to Calm Thinking analyses,
# we decided not to analyze Calm Thinking data as part of the dissertation's scope

# Define function to check that selected intermediate clean CSV data files from
# Calm Thinking study contain those relevant to present paper (for full set 
# of intermediate clean CSV data files, see https://github.com/TeachmanLab/MT-Data-CalmThinkingStudy)

check_relevant_files <- function(filenames) {
  relevant_files <- c("bbsiq.csv", "demographics.csv", "mechanisms.csv", "oa.csv", 
                      "rr.csv", "task_log.csv", "dass21_as.csv", "demographics_race.csv", 
                      "participant.csv", "study.csv")
  
  if (all(relevant_files %in% filenames) == FALSE) {
    missing_files <- setdiff(relevant_files, filenames)
    
    warning(paste0(c("You are missing these files:", paste0(" ", missing_files))))
  }
}

# ---------------------------------------------------------------------------- #
# Define convert_POSIXct() ----
# ---------------------------------------------------------------------------- #

# Define function to convert system-generated timestamps from Calm Thinking study 
# to POSIXct data types (with "tz = 'UTC'" for user-provided "return_date_as_POSIXct" 
# of "return_intention" table and "tz = 'EST'" for all system-generated timestamps)

convert_POSIXct <- function(dat) {
  for (i in 1:length(dat)) {
    POSIXct_colnames <- c(names(dat[[i]])[grep("as_POSIXct", names(dat[[i]]))],
                          "system_date_time_earliest",
                          "system_date_time_latest")
    
    for (j in 1:length(POSIXct_colnames)) {
      # Strip timezone from character vector
      
      dat[[i]][, POSIXct_colnames[j]] <- sub(" UTC| EST", "", 
                                             dat[[i]][, POSIXct_colnames[j]])
      
      # Convert character vector to POSIXct, specifying timezone
      
      if (names(dat[i]) == "return_intention" & 
          POSIXct_colnames[j] == "return_date_as_POSIXct") {
        dat[[i]][, POSIXct_colnames[j]] <- as.POSIXct(dat[[i]][, POSIXct_colnames[j]],
                                                      format = "%Y-%m-%d %H:%M:%S",
                                                      tz = "UTC")
      } else {
        dat[[i]][, POSIXct_colnames[j]] <- as.POSIXct(dat[[i]][, POSIXct_colnames[j]],
                                                      format = "%Y-%m-%d %H:%M:%S",
                                                      tz = "EST")
      }
    }
  }
  
  return(dat)
}

# ---------------------------------------------------------------------------- #
# Define format_dem_tbl() ----
# ---------------------------------------------------------------------------- #

# Define function to format demographics table

format_dem_tbl <- function(dem_tbl, footnotes, title, country_other_row_idx) {
  # Format "label" column using Markdown
  
  dem_tbl$label_md <- dem_tbl$label
  
  rows_no_indent <- dem_tbl$label_md == "n" | 
    grepl("\\b(Age|Gender|Race|Ethnicity|Country|Education|Employment|Annual|Marital)\\b",
          dem_tbl$label_md)
  rows_indent <- !rows_no_indent
  
  indent_spaces <- "\\ \\ \\ \\ \\ "
  
  dem_tbl$label_md[rows_indent] <- paste0(indent_spaces, dem_tbl$label_md[rows_indent])
  
  dem_tbl$label_md[dem_tbl$label_md == "n"] <- "*n*"
  dem_tbl$label_md <- gsub("n \\(%\\)", "*n* \\(%\\)", dem_tbl$label_md)
  
  dem_tbl$label_md <- gsub("M \\(SD\\)", "*M* \\(*SD*\\)", dem_tbl$label_md)
  
  dem_tbl <- dem_tbl[c("label_md", names(dem_tbl)[names(dem_tbl) != "label_md"])]
  
  # Identify rows for footnotes
  
  if (dem_tbl$label[country_other_row_idx] != "Other") {
    stop("Row index for 'country' value of 'Other' is incorrect. Update 'country_other_row_idx'.")
  }
  
  # Define columns
  
  left_align_body_cols <- "label_md"
  target_cols <- names(dem_tbl)[names(dem_tbl) != "label"]
  
  # Create flextable
  
  dem_tbl_ft <- flextable(dem_tbl[, target_cols]) |>
    set_table_properties(align = "left") |>
    
    set_caption(as_paragraph(as_i(title)), word_stylename = "heading 1",
                fp_p = fp_par(padding.left = 0, padding.right = 0),
                align_with_table = FALSE) |>
    
    align(align = "center", part = "header") |>
    align(align = "center", part = "body") |>
    align(j = left_align_body_cols, align = "left", part = "body") |>
    align(align = "left", part = "footer") |>
    
    valign(valign = "bottom", part = "header") |>
    
    set_header_labels(label_md = "Characteristic",
                      value    = "Value") |>
    
    colformat_md(j = "label_md", part = "body") |>
    
    add_footer_lines(gen_note) |>
    
    footnote(i = country_other_row_idx,
             j = 1,
             value = as_paragraph_md(footnotes$country_other),
             ref_symbols = " a",
             part = "body") |>
    
    autofit()
}