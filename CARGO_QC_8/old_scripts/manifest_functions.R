library(tidyverse)
library(tibble)
library(shiny)

col_names <- c("Row", "Group ID", "Sample ID", "Plate Barcode or Number", "Plate", "Well", "Sample Well", "Species",
               "Gender (M/F/U)", "Volume (ul)", "Concentration (ng/ul)", "OD 260/280", "Tissue Source",
               "Extraction Method", "Ethnicity", "Parent 1 ID", "Parent 2 ID", "Replicate(s) ID", "Cancer Sample (Y/N)")

controls <- c("Hypo-methylated Control", "Hyper-methylated control")

add_column_na <- function(d, col_names) {
  add_cols <- col_names[!col_names %in% colnames(d)]
  
  if(length(add_cols) != 0) d[add_cols] <- NA
  d
}

get_plate_wells <- function(n = 1, r = 12, c = 8, trans = FALSE) {
  if (n == 0) { return(NULL); }
  else if (n < 1) { return( get_plate_wells(n = 1, r = r * n, c) ) }
  else { return( c(
    paste0(rep(LETTERS[1:8], each = r), formatC(rep(1:r, times = c), width = 2, flag = "0")),
    get_plate_wells(n-1, r, c)
  ) ) }
}

get_plate_wells_trans <- function(n = 1, transpose = FALSE, r = 12, c = 8) {
  if (n == 0) { return(NULL); }
  else if (n < 1) { return( get_plate_wells(n = 1, r = r * n, c) ) }
  else { 
    well_rows <- if_else(rep(transpose, times = r * c),
                         rep(LETTERS[1:8], times = r),
                         rep(LETTERS[1:8], each = r))
    well_columns <- if_else(rep(transpose, times = r * c),
                            rep(1:r, each = c),
                            rep(1:r, times = c))
    return( c( str_c(well_rows, formatC(well_columns, width = 2, flag = "0")),
    get_plate_wells_trans(n - 1, transpose, r, c) ) )
  }
}

get_info <- function(samples, controls, plate_size, chip_size) {
  n_samples <- nrow(samples)
  n_controls <- length(controls)
  total_plates <- ceiling( n_samples / ( plate_size - n_controls ) )
  used_wells <- n_samples + (total_plates * n_controls)
  
  total_chips <- ceiling( used_wells / chip_size)
  empty_wells <- ( total_chips * chip_size ) - used_wells
  samples_per_plate <- plate_size - n_controls
  total_controls <- total_plates * n_controls
  
  lst(plate_size, chip_size, n_samples, n_controls, samples_per_plate,
      total_plates, total_chips, total_controls, empty_wells)
}

create_manifest_2 <- function(fromQC, controls, id_col, by_cols, col_vals, add_cols) {
  info <- get_info(fromQC, controls, 96, 8)
  # Fill empty wells
  #### TODO: Add UI to decide what to do with empty wells
  samples <- get_data() %>%
    mutate("Sample ID" = as.character(!!! syms(id_col))) %>%
    bind_rows(tibble("Sample ID" = rep(controls, length.out = info$empty_wells)))

  # Distribute samples across plates and add plate controls
  samples_w_controls <- samples %>%
    sample_n(n()) %>%
    mutate(Plate = rep(1:info$total_plates, each = info$samples_per_plate, length.out = n())) %>%
    bind_rows(tibble("Sample ID" = rep(controls, length.out = info$total_controls),
                   Plate = rep(1:info$total_plates, each = info$n_controls, length.out = info$total_controls)))

  # Arrange by plates, randomize and add well information
  plated_samples <- samples_w_controls %>%
    group_split(Plate) %>%
    map(~ sample_n(., n())) %>%
    map(~ mutate(., Well = get_plate_wells_trans(1, transpose = TRUE) %>% head(n()))) %>%
    reduce(rbind)
  
  # Format for manifest output
  formatted_samples <- plated_samples %>%
    mutate(!!! col_vals) %>%
    add_column_na(col_names) %>%
    select(union(col_names, c(by_cols, add_cols)))
  
  return(formatted_samples)
}

