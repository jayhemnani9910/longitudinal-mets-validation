# Utility for loading and stacking NHANES tables across cycles.
# Tables are downloaded by scripts/01_download_nhanes.R into data/raw/nhanes/.

suppressMessages({
  library(dplyr)
  library(purrr)
})

#' Load and stack one NHANES table across all cycles
#'
#' @param table Table name without cycle suffix (e.g., "DEMO", "BMX")
#' @return tibble with all available cycles row-stacked plus a `cycle` column
#'   identifying the survey years (e.g., "1999-2000", "2001-2002").
load_nhanes_table <- function(table) {
  cycles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  cycle_years <- c(
    "1999-2000", "2001-2002", "2003-2004", "2005-2006",
    "2007-2008", "2009-2010", "2011-2012", "2013-2014",
    "2015-2016", "2017-2018"
  )
  paths <- file.path("data/raw/nhanes", sprintf("%s_%s.rds", table, cycles))
  existing <- file.exists(paths)
  if (!any(existing)) {
    stop(sprintf("No files found for table '%s' in data/raw/nhanes/", table))
  }
  # Load each cycle, stripping factors/labels to base types
  dfs <- map2(paths[existing], cycle_years[existing], function(p, yr) {
    df <- readRDS(p)
    df$cycle <- yr
    for (col in names(df)) {
      v <- df[[col]]
      if (inherits(v, "haven_labelled")) v <- haven::zap_labels(v)
      if (is.factor(v))                 v <- as.character(v)
      df[[col]] <- v
    }
    df
  })

  # Unify column types across cycles so bind_rows succeeds.
  # If any cycle stores a column as character, coerce all cycles to character.
  all_cols <- unique(unlist(lapply(dfs, names)))
  for (col in all_cols) {
    types <- sapply(dfs, function(d) if (col %in% names(d)) class(d[[col]])[1] else NA_character_)
    types <- types[!is.na(types)]
    if (length(unique(types)) > 1 && "character" %in% types) {
      dfs <- lapply(dfs, function(d) {
        if (col %in% names(d)) d[[col]] <- as.character(d[[col]])
        d
      })
    }
  }

  dplyr::bind_rows(dfs)
}

#' Safe column accessor that returns NA when a column is missing
#'
#' Some NHANES columns are renamed across cycles. Use this in transmute() to
#' tolerate missing columns gracefully.
get_or_na <- function(df, colname) {
  if (colname %in% names(df)) df[[colname]] else NA
}
