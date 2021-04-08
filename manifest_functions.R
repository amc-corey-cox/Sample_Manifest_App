col_names <- c("Row", "Group ID", "Sample ID", "Plate Barcode or Number", "Plate", "Well", "Sample Well", "Species",
               "Gender (M/F/U)", "Volume (ul)", "Concentration (ng/ul)", "OD 260/280", "Tissue Source",
               "Extraction Method", "Ethnicity", "Parent 1 ID", "Parent 2 ID", "Replicate(s) ID", "Cancer Sample (Y/N)")

add_column_na <- function(d, col_names) {
  add_cols <- col_names %>% setdiff(colnames(d))
  if(length(add_cols) != 0) d[add_cols] <- NA
  d
}

disperse <- function(df_a, df_b) {
  if(is.null(df_a) | is.null(df_b)) { return(bind_rows(df_a, df_b)) }
  if(nrow(df_b) > nrow(df_a)) { return(disperse(df_b, df_a)) }
  a <- nrow(df_a)
  b <- nrow(df_b)
  
  bind_rows(df_a, df_b) %>% arrange(c(seq(a), a * seq(b)/(b + 1)))
}

get_wells <- function(n = 1, transpose = FALSE, r = 12, c = 8) {
  if (n == 0) { return(NULL); }
  else if (n < 1) { return( get_wells(n = 1, r = r * n, c) ) }
  else {
    well_rows <- if_else(rep(transpose, times = r * c),
                         rep(LETTERS[1:8], times = r),
                         rep(LETTERS[1:8], each = r))
    well_columns <- if_else(rep(transpose, times = r * c),
                            rep(1:r, each = c),
                            rep(1:r, times = c))
    return( c( str_c(well_rows, formatC(well_columns, width = 2, flag = "0")),
    get_wells(n - 1, transpose, r, c) ) )
  }
}

get_info <- function(samples, controls, plate_size, chip_size) {
  lst( plate_size, chip_size,
    n_samples = nrow(samples),
    n_controls = length(controls),
    total_plates = ceiling( n_samples / ( plate_size - n_controls ) ),
    used_wells = n_samples + (total_plates * n_controls),
    
    total_chips = ceiling( used_wells / chip_size ),
    empty_wells = ( total_chips * chip_size ) - used_wells,
    
    samples_per_plate = plate_size - n_controls,
    total_controls = total_plates * n_controls)
}

format_manifest <- function(samples, by_cols, add_cols, col_vals = NULL) {
  samples_w_wells <- samples %>%
    group_split(Plate) %>%
    map(~ mutate(., Well = get_wells(1, transpose = TRUE) %>% head(n()))) %>%
    bind_rows
  
  samples_w_wells %>%
    mutate(!!! col_vals) %>%
    mutate("Gender (M/F/U)" = Gender) %>%
    add_column_na(col_names) %>%
    select(union(col_names, c(all_of(by_cols), add_cols)))
}

simple_disperse <- function(samples, controls, seed, id_col, by_cols, empty_wells) {
  set.seed(seed)
  ### TODO: get plate dimensions from UI
  info <- get_info(samples, controls, 96, 8)

  randomized_samples <- samples %>% sample_n(n()) %>%
    mutate("Sample ID" = as.character(!!! syms(id_col))) %>%
    group_split(!!! syms(by_cols)) %>% sample()

  if (info$empty_wells > 0 & empty_wells == "Use Controls") {
    empty <- tibble("Sample ID" = rep(controls, length.out = info$empty_wells))
    dispersed_samples <- randomized_samples %>%
      list_modify(empty = empty) %>% reduce(disperse)
  }
  else {
    empty <- tibble("Sample ID" = rep("Empty", length.out = info$empty_wells))
    dispersed_samples <- randomized_samples %>%
      reduce(disperse) %>% bind_rows(empty)
  }

  plated_samples <- dispersed_samples %>%
    mutate(Plate = rep(1:info$total_plates, each = info$samples_per_plate, length.out = n())) %>%
    group_split(Plate) %>% imap(~ disperse(.x, tibble("Sample ID" = controls, Plate = .y))) %>%
    bind_rows

  plated_samples %>%
    mutate(Chip = rep(1:info$total_chips, each = info$chip_size, length.out = n())) %>%
    group_split(Chip) %>% map(~ sample_n(., n())) %>%
    bind_rows
}

col_split <- function(.data, split_on) {
  len <- length(split_on)
  if(len == 0) { return(.data) }
  map_depth(col_split(.data, split_on[-len]), .depth = len - 1, ~ group_split(., !! sym(split_on[len])))
}

multi_reduce <- function(.x, .f) {
  if(is_tibble(.x[[1]])) { return(reduce(.x, .f)) }
  map(.x, ~ multi_reduce(.x, .f)) %>% reduce(.f)
}

# Testing new dispersal mechanism
grouped_disperse <- function(samples, controls, seed, id_col, by_cols, empty_wells) {
  set.seed(seed)
  ### TODO: get plate dimensions from UI
  info <- get_info(samples, controls, 96, 8)
  
  randomized_samples <- samples %>% sample_n(n()) %>%
    mutate("Sample ID" = as.character(!!! syms(id_col))) %>%
    col_split(by_cols) %>% multi_reduce(disperse)
  
  if (info$empty_wells > 0 & empty_wells == "Use Controls") {
    empty <- tibble("Sample ID" = rep(controls, length.out = info$empty_wells))
    dispersed_samples <- disperse(randomized_samples, empty)
  } else {
    empty <- tibble("Sample ID" = rep("Empty", length.out = info$empty_wells))
    dispersed_samples <-bind_rows(randomized_samples, empty)
  }
  
  plated_samples <- dispersed_samples %>%
    mutate(Plate = rep(1:info$total_plates, each = info$samples_per_plate, length.out = n())) %>%
    group_split(Plate) %>% imap(~ disperse(.x, tibble("Sample ID" = controls, Plate = .y))) %>%
    bind_rows
  
  plated_samples %>%
    mutate(Chip = rep(1:info$total_chips, each = info$chip_size, length.out = n())) %>%
    group_split(Chip) %>% map(~ sample_n(., n())) %>%
    bind_rows
}

plate_randomize <- function(samples, controls, seed, id_col, by_cols, empty_wells) {
  set.seed(seed)
  
  ### TODO: get plate dimensions from UI
  info <- get_info(samples, controls, 96, 8)
  
  randomized_samples <- samples %>% sample_n(n()) %>%
    mutate("Sample ID" = as.character(!!! syms(id_col)))
  
  if (info$empty_wells > 0 & empty_wells == "Use Controls") {
    empty <- tibble("Sample ID" = rep(controls, length.out = info$empty_wells))
    dispersed_samples <- disperse(randomized_samples, empty)
  } else {
    empty <- tibble("Sample ID" = rep("Empty", length.out = info$empty_wells))
    dispersed_samples <-bind_rows(randomized_samples, empty)
  }
  
  dispersed_samples %>% mutate(Plate = rep(1:info$total_plates, each = info$samples_per_plate, length.out = n())) %>%
    group_split(Plate) %>% imap(~ disperse(.x, tibble("Sample ID" = controls, Plate = .y))) %>%
    bind_rows
}