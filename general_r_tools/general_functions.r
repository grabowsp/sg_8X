# General purpose functions used by many scripts

add_slash <- function(dir_string){
  # Add '/' to end of a string, particularly directory names that will be
  #  used with paste() to generate file names
  # INPUTS
  # dir_string = character string; ex: '/FULL/DIR/STRING'
  # OUPUT
  # character string
  #   if dir_string already ends with '/', then returns dir_string
  #   if dir_string does not end with '/', returns dir_string ending with '/'
  #########
  final_string <- dir_string
  last_char <- rev(unlist(strsplit(dir_string, split = '')))[1]
  if(last_char != '/'){
    final_string <- paste(final_string, '/', sep = '')
  }
  return(final_string)
}

gen_unified_df <- function(df_list){
  # Generate a combined data.frame from a list of data.frames with the
  #   same formats
  #  Example: df_list has 4 elements, each with a data.frame with same
  #    colnames generated for a different sample or comparison; this
  #    generates a combined data.frame with the same number of columns,
  #    but column 1 from all elements are concatenated together, same for
  #    column 2, etc.
  # INPUTS
  # df_list = list where each element is a data.frame and all the data.frames
  #            have the same number of columns and same column names
  # OUTPUT
  # data.frame with same number of columns and same column names as the 
  #  data.frames in the elements of df_list; column contains the
  #  concatenated vector of all the column 1's from each element, continuing
  #  for all the columns
  ###########################
  n_cols <- ncol(df_list[[1]])
  df_colnames <- colnames(df_list[[1]])
  tot_list <- list()
  for(i in seq(n_cols)){
    tot_list[[i]] <- unlist(lapply(df_list, function(x) x[,i]))
  }
  tot_df <- data.frame(tot_list, stringsAsFactors = F)
  colnames(tot_df) <- colnames(df_list[[1]])
  return(tot_df)
}

generate_dist_window_df <- function(dist_vec, value_vec, window_size, 
  avg_type = mean){
  # Get average of all values in value vec in each window interval
  #  Written to generate window averages of r^2, but can be used
  #  for any type of pairwise comparisons across SNPs
  # INPUTS
  # dist_vec = vector of distances between SNPs used for comparison
  # value_vec = values (ex: r^2); must have same length to and correspond with
  #              dist_vec
  # window_size = bp size of the window
  # ave_type = the type of average, though could also probably use 
  #               'sum' or 'max'
  # OUTPUT
  # data.frame with max distance in a window, the average value in the 
  #  window, and the number of values in the window
  ##############
  max_dist <- max(dist_vec)
  dist_mult_vec <- seq(ceiling(max_dist/window_size))
  wind_dists <- lapply(dist_mult_vec, function(x) 
    ((x-1)*window_size) + c(1:window_size))
  window_inds <- lapply(wind_dists, function(x)
    which(dist_vec %in% x))
  num_window_inds <- unlist(lapply(window_inds, length))
  window_avg <- unlist(lapply(window_inds, function(x)
    avg_type(value_vec[x], na.rm = T)))
  tot_df <- data.frame(window_pos = dist_mult_vec * window_size,
    window_avg = window_avg, n_vals = num_window_inds, stringsAsFactors = F)
  return(tot_df)
}

generate_perHighVal_window_df <- function(dist_vec, value_vec, window_size,
  max_val){
  # Get the percentage of values in value vec above max_val in each window;
  #  Written to generate the % of r2 values above a high r2 in windows, but
  #   can be used for any tuype of pairwise comparisons across SNPs
  # INPUTS
  # dist_vec = vector of distances between SNPs used for comparison
  # value_vec = values (ex: r^2); must have same length to and correspond with
  #              dist_vec
  # window_size = bp size of the window
  # max_val = the maximum value above which a value is considered "High"
  # OUTPUT
  # data.frame with max distance in a window, the percentage of "High Values"
  #  in the window, the number of comparisons in the window, and the number
  #  of NA's
  ###############
  max_dist <- max(dist_vec)
  dist_mult_vec <- seq(ceiling(max_dist/window_size))
  wind_dists <- lapply(dist_mult_vec, function(x)
    ((x-1)*window_size) + c(1:window_size))
  window_inds <- lapply(wind_dists, function(x)
    which(dist_vec %in% x))
  num_window_inds <- unlist(lapply(window_inds, length))
  n_na_inds <- unlist(lapply(window_inds, function(x)
    sum(is.na(value_vec[x]))))
  n_high_ld_inds <- unlist(lapply(window_inds, function(x)
    sum(value_vec[x] >= max_val, na.rm = T)))
  per_high_ld <- n_high_ld_inds/(num_window_inds-n_na_inds)
  tot_df <- data.frame(window_pos = dist_mult_vec * window_size,
    per_high_ld = per_high_ld, n_comps = num_window_inds, n_nas = n_na_inds,
    stringsAsFactors = F)
  return(tot_df)
}


