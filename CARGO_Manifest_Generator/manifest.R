#SHINY Get this from UI, use template radio button.
col_names <- c("Row", "Group ID", "Sample ID", "Plate Barcode or Number", "Plate", "Well", "Sample Well", "Species",
               "Gender (M/F/U)", "Volume (ul)", "Concentration (ng/ul)", "OD 260/280", "Tissue Source",
               "Extraction Method", "Ethnicity", "Parent 1 ID", "Parent 2 ID", "Replicate(s) ID", "Cancer Sample (Y/N)")

get_n_plates <- function (n_samples, n_controls, n_plate_wells, n_chip_wells) {
  total_plates <- ceiling( n_samples / ( n_plate_wells - n_controls ) )
  total_controls <- total_plates * n_controls
  used_wells <- n_samples + total_controls
  
  total_chips <- ceiling( used_wells / n_chip_wells )
  empty_wells <- total_chips * n_chip_wells - used_wells
  
  ( used_wells + empty_wells ) / n_plate_wells
}

add_column_na <- function(d, col_names) {
  add_cols <- col_names[!col_names %in% colnames(d)]
  
  if(length(add_cols) != 0) d[add_cols] <- NA
  d
}

get_plates <- function(n_plates, n_wells = 96) {
  rep(1:ceiling(n_plates), each = n_wells, length.out = n_wells * n_plates)
}

get_plate_wells <- function(n = 1, r = 12, c = 8) {
  if (n == 0) { return(NULL); }
  else if (n < 1) { return( get_plate_wells(n = 1, r = r * n, c) ) }
  else { return( c(
    paste0(rep(LETTERS[1:8], each = r), formatC(rep(1:r, times = c), width = 2, flag = "0")),
    get_plate_wells(n-1, r, c)
  ) ) }
}

get_plate_chips <- function(wells) {
  as.numeric(substring(wells, 2))
}

get_controls <- function(control_list, whole_plates, fill = 0) {
  n_controls <- whole_plates * length(control_list) + fill
  
  tibble("Sample ID" = rep(control_list, length.out = n_controls)) %>%
    mutate(r = sample(rep(sample(LETTERS[1:8]), length.out = n_controls)),
           c = sample(1:12, size = n_controls),
           n = rep(1:whole_plates, each = n_controls),
           Well = paste0(r, formatC(c, width = 2, flag = "0"))
    ) %>%
    select (-r, -c) %>% group_split(n, keep = FALSE) %>%
    imap(~ mutate(.x, Plate = .y, Chip = as.numeric(substring(Well, 2))))
}

get_open_wells <- function(all_controls, n_plates) {
  tibble(Plate = get_plates(n_plates), Well = get_plate_wells(n_plates), Chip = get_plate_chips(Well)) %>%
    group_split(Plate) %>%
    map2(all_controls, ~ filter(.x, ! Well %in% .y$Well)) %>%
    reduce(rbind) %>%
    arrange(Well)
}

plate_samples <- function(pheno, open_wells, by_cols, col_vals) {
  pheno %>% slice(sample(1:n(), n())) %>%
    mutate(!!! col_vals) %>%
    arrange(!!! quos(!!! syms(by_cols))) %>%
    cbind(open_wells) %>%
    group_split(Plate, Chip) %>%
    map(~ slice(., sample(1:n(), n()))) %>%
    reduce(rbind)
}

create_manifest <- function(pheno, controls, by_cols, col_vals) {
  set.seed(42)  # RSHINY Add control to allow changing seed
  
  wells_per_plate <- 96
  samples_per_chip <- 8
  
  n_plates <- get_n_plates(nrow(pheno), length(controls), 96, 8)
  whole_plates <- ceiling (n_plates)
  fill_wells <- (96 - length(controls)) * whole_plates - nrow(pheno)
  
  all_controls <- get_controls(controls, whole_plates, fill_wells)
  open_wells <- get_open_wells(all_controls, n_plates)
  plated_samples <- plate_samples(pheno, open_wells, by_cols, col_vals)
  formatted_plates <- plated_samples %>% add_column_na(col_names) %>% select(union(col_names, by_cols))
  return(formatted_plates)
}