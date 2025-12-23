# General helper functions --------------------------------------------------------

find_ladder_peaks <- function(ladder_df,
                              n_reference_sizes,
                              minimum_ladder_signal,
                              sample_id) {

  if(!is.na(minimum_ladder_signal)){
    peaks <- pracma::findpeaks(ladder_df$smoothed_signal,
      peakpat = "[+]{3,}[0]*[-]{3,}" # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
      )

    ladder_subset_df <- ladder_df[peaks[, 2], ]
    ladder_peaks <- ladder_subset_df[which(ladder_subset_df$signal > minimum_ladder_signal), "scan"]

  } else{
    median_signal <- median(ladder_df$smoothed_signal, na.rm = TRUE)
    sd_signal <- stats::sd(ladder_df$smoothed_signal, na.rm = TRUE)
  
    ladder_peaks <- vector("numeric")
    ladder_peak_threshold <- 1

    while (length(ladder_peaks) < n_reference_sizes) {
      peaks <- pracma::findpeaks(ladder_df$smoothed_signal,
                                 peakpat = "[+]{3,}[0]*[-]{3,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                                 minpeakheight = median_signal + sd_signal * ladder_peak_threshold
      )
  
      ladder_peaks <- ladder_df$scan[peaks[, 2]]
  
      # lower the threshold for the next cycle
      ladder_peak_threshold <- ladder_peak_threshold - 0.05
  
      # provide an exit if there are not enough peaks found
      if (sd_signal * ladder_peak_threshold <= 0) {
        break
      }
    }
  }

  # go through raw signal and make sure that the identified scan in the smoothed signal is still the highest
  # it will also deal with cases where the scans have the same signal (which.max will chose first)
  n_scans <- length(ladder_df$scan)
  window_width <- 3
  peak_position <- numeric(length(ladder_peaks))
  for (i in seq_along(peak_position)) {
    if (ladder_peaks[i] + window_width > 1 & ladder_peaks[i] + window_width < n_scans) { # make sure that the subsetting would be in bounds when taking window into account
      max_peak <- which.max(ladder_df$signal[(ladder_peaks[i] - window_width):(ladder_peaks[i] + window_width)])

      peak_position[i] <- ladder_peaks[i] - window_width - 1 + max_peak
    } else {
      peak_position[i] <- ladder_peaks[i]
    }
  }

  return(ladder_peaks)
}

mean_rsq <- function(scan, size, choose){

  # set cor window size to be one less than selected window size for ladder chunking
  window_size = choose -1

  cors <- vector("numeric", length = length(scan) - window_size)

  for (i in seq_along(cors)) {
    xi <- scan[i:(i + window_size)]
    yi <- size[i:(i + window_size)]
    cors[i] <- stats::cor(yi, xi)^2
  }

  mean_rsq <- mean(cors)

  return(mean_rsq)

}

ladder_iteration <- function(reference_sizes, observed_sizes, choose = 5, 
                              max_combinations = 10000, r_squared_threshold = 0.9, top_n_branching = 3) {

  # Initialize counter for total combinations tested in parent environment
  total_combinations_tested <- 0

  find_best_combinations <- function(recombinations, reference_sizes, top_n) {
    rsq_vector <- vector("numeric", ncol(recombinations))
    for (i in 1:ncol(recombinations)) {
      rsq_vector[[i]] <- stats::cor(reference_sizes, recombinations[, i])^2
    }
    
    # Return top N combinations
    top_indices <- order(rsq_vector, decreasing = TRUE)[1:min(top_n, length(rsq_vector))]
    return(list(
      recombinations = recombinations[, top_indices, drop = FALSE],
      rsq = rsq_vector[top_indices]
    ))
  }
  
  # Recursive function to explore assignments with backtracking
  explore_assignments <- function(remaining_ref, remaining_obs, 
                                 current_assigned_ref, current_assigned_obs, 
                                 current_rsq, choose, depth = 1) {
    
    if(length(remaining_ref) != 0 && length(remaining_ref) == length(remaining_obs)){
      current_assigned_obs <- c(current_assigned_obs, remaining_obs)
      current_assigned_ref <- c(current_assigned_ref, remaining_ref)
      # set to empty so recursive function ends
      remaining_ref <- vector("numeric")
      remaining_obs <- vector("numeric")
    }

    
    # Base case: all reference sizes assigned
    if (length(remaining_ref) == 0) {
      output_list <- list(
        assigned_observed = current_assigned_obs,
        assigned_reference = current_assigned_ref,
        final_rsq = current_rsq
      )

      return(output_list)
    }
    
    # Calculate window size for this iteration
    n_obs <- length(remaining_obs)
    n_ref <- length(remaining_ref)

    # If few reference sizes remain (less than 4 + choose), take all remaining at once
    if (n_ref < 4 + choose && n_ref <= n_obs) {
      current_choose <- n_ref
    } else {
      current_choose <- min(choose, n_ref)
    }
    
    # Calculate start window
    remainder <- n_ref - current_choose
    start_window <- min(n_obs - remainder, n_obs)
    start_window <- max(start_window, current_choose) # Ensure we have enough points
    
    # Generate combinations
    n_recombinations <- choose(length(remaining_obs[1:start_window]), current_choose)

    if (n_recombinations > max_combinations) {
      stop(
        call. = FALSE,
        paste0("Too many combinations to test (", n_recombinations, "). Adjust parameters or max_combinations.")
      )
    }

    # Increment the total number of combinations tested
    total_combinations_tested <<- total_combinations_tested + n_recombinations
   
    recombinations <- utils::combn(remaining_obs[1:start_window], current_choose)
   
    # Get top combinations
    top_combinations <- find_best_combinations(
      recombinations = recombinations,
      reference_sizes = remaining_ref[1:current_choose],
      top_n = top_n_branching
    )
    
    best_result <- NULL
    best_rsq <- -Inf

    # filter for rsq but keep best
    if (depth > 1) {
      recombinations_to_keep <- top_combinations$recombinations[, which(top_combinations$rsq > r_squared_threshold | top_combinations$rsq == max(top_combinations$rsq)), drop = FALSE]
    } else{
      recombinations_to_keep <- top_combinations$recombinations
    }
    
    # Try each top combination
    for (i in 1:ncol(recombinations_to_keep)) {
      selected_comb <- recombinations_to_keep[, i]

      # Find positions in remaining observations
      positions <- sapply(selected_comb, function(x) which(remaining_obs == x)[1])
      
      # Set up for next iteration
      new_assigned_obs <- c(current_assigned_obs, selected_comb)
      new_assigned_ref <- c(current_assigned_ref, remaining_ref[1:current_choose])
      
      # Calculate overall R² so far
      overall_rsq <- if (length(new_assigned_obs) > 2) {
        mean_rsq(new_assigned_ref, new_assigned_obs, choose)
      } else {
        1.0 # Not enough points for correlation
      }

      last_selected_obs_position <- which(remaining_obs == selected_comb[length(selected_comb)])
      if(last_selected_obs_position == length(remaining_obs)){
        new_remaining_obs <- numeric()
      } else{
        new_remaining_obs <- remaining_obs[(last_selected_obs_position + 1):length(remaining_obs)]
      }
  
      last_selected_reference_position <- which(remaining_ref == remaining_ref[current_choose])
      if(last_selected_reference_position == length(remaining_ref)){ 
        new_remaining_ref <- numeric()
      } else{
        new_remaining_ref <- remaining_ref[(last_selected_reference_position + 1):length(remaining_ref)]
      }
  
      reference_sizes_left <- length(remaining_ref) - last_selected_reference_position

      # deal with situation of just a couple of reference peaks left over, but you need a good amount to make correlation
      # therefore increase chose to the max
      if (reference_sizes_left - choose < 3) {
        new_choose <- reference_sizes_left
      } else{
        new_choose <- choose
      }
      
      result <- explore_assignments(
        new_remaining_ref, new_remaining_obs,
        new_assigned_ref, new_assigned_obs,
        overall_rsq, new_choose, depth + 1
      )
      
      if (!is.null(result) && result$final_rsq > best_rsq) {
        best_result <- result
        best_rsq <- result$final_rsq
      }
    }
   
    
    return(best_result)
  }

  if(length(observed_sizes) > length(reference_sizes)){

    result <- explore_assignments(
      reference_sizes, observed_sizes,
      vector("numeric"), vector("numeric"),
      1.0, choose = choose
    )

  } else{

    result <- explore_assignments(
      observed_sizes, reference_sizes,
      vector("numeric"), vector("numeric"),
      1.0, choose = choose
    )

    #need to switch around sizes and scans
    result[c("assigned_reference", "assigned_observed")] <- result[c("assigned_observed", "assigned_reference")]

  }
  
  if (is.null(result)) {
    stop("Could not find a valid assignment")
  }

  # Add the total combinations tested to the result
  final_result <- data.frame(scan = result$assigned_observed, size = result$assigned_reference)
  attr(final_result, "total_combinations_tested") <- total_combinations_tested
  
  return(final_result)

}


exhaustive_ladder_matching <- function(reference_sizes, observed_sizes, max_combinations = 1e6) {

  find_best_combination <- function(ref, obs){
    # Check that we're not exceeding computational limits
    n_combinations <- choose(length(obs), length(ref))
    if(n_combinations > max_combinations) {
      stop(paste0("Too many combinations to test (", n_combinations, 
                  "). Adjust inputs or increase max_combinations."))
    }
    
    # Generate all possible combinations
    combination_indices <- utils::combn(length(obs), length(ref))
    
    # Evaluate each combination
    best_rsq <- -Inf
    best_combo <- NULL
    
    for(i in 1:ncol(combination_indices)) {
      # Get the indices for this combination
      indices <- combination_indices[, i]
      
      # Extract corresponding values
      selected_values <- obs[indices]
      
      # Calculate R-squared for this assignment
      current_rsq <- mean_rsq(selected_values, ref, length(ref))
      
      # Track best combination
      if(current_rsq > best_rsq) {
        best_rsq <- current_rsq
        best_combo <- indices
      }

    }

    best_comb <- data.frame(
      scan = obs[best_combo],
      size = ref
    )
    attr(best_comb, "total_combinations_tested") <- n_combinations

    return(best_comb)
  }

  if(length(observed_sizes) > length(reference_sizes)){
    result <- find_best_combination(
      ref = reference_sizes, 
      obs = observed_sizes
    )    
  } else{
    result <- find_best_combination(
      ref = observed_sizes, 
      obs = reference_sizes
    )    
    colnames(result) <- rev(colnames(result))
  }
  
  return(result)
}
  



# predict bp size

predict_bp_size <- function(
    ladder_df,
    scans) {
  ladder_df <- ladder_df[which(!is.na(ladder_df$size)), ]
  ladder_df <- ladder_df[which(!is.na(ladder_df$scan)), ]

  n_knots <- ifelse(
    nrow(ladder_df) > 10, 
    -1, #the default setting
    floor(nrow(ladder_df) / 2)
)

  p_spline_model <- mgcv::gam(size ~ s(scan, bs = "cr", k = n_knots), data = ladder_df)
  predicted_size <- predict(p_spline_model, newdata = data.frame(scan = scans))
  
  return(predicted_size)
}



ladder_fit_cor <- function(fragments){
  ladder_df <- fragments$ladder_df[order(fragments$ladder_df$size),]
  ladder_df <- ladder_df[which(!is.na(ladder_df$size)), ]

  # Function to calculate the fitting constants for each group of three neighboring points
  cor_list <- vector("list", length = nrow(ladder_df) - 2)

  for (i in seq_along(cor_list)) {
    xi <- ladder_df$scan[i:(i + 2)]
    yi <- ladder_df$size[i:(i + 2)]
    cor_list[[i]] <- list(
      rsq = stats::cor(yi, xi)^2,
      size_ranges = yi
    )
  }

  return(cor_list)
}


ladder_rsq_warning_helper <- function(
    fragments,
    rsq_threshold) {
  
  cor_list <- ladder_fit_cor(fragments)
  rsq <- sapply(cor_list, function(x) x$rsq)

  if (any(rsq < rsq_threshold)) {
    size_ranges <- sapply(cor_list, function(x) x$size_ranges)
    size_ranges <- size_ranges[, which(rsq < rsq_threshold), drop = FALSE]
    size_ranges_vector <- vector("numeric", ncol(size_ranges))
    for (j in seq_along(size_ranges_vector)) {
      size_ranges_vector[j] <- paste0(size_ranges[1, j], "-", size_ranges[3, j])
    }
    warning(
      call. = FALSE,
      paste(
        "sample", fragments$unique_id, "has badly fitting ladder for bp sizes:",
        paste0(size_ranges_vector, collapse = ", ")
      )
    )
  }
}

# ladder ------------------------------------------------------------------


#' Ladder and bp sizing
#'
#' Find the ladder peaks in and use that to call bp size
#'
#' @param fragments_list list from 'read_fsa' function
#' @param config A trace_config object generated using [load_config()].
#' @param ... additional parameters from any of the functions in the pipeline detailed below may be passed to this function. This overwrites values in the `config`. These parameters include:
#'   \itemize{
#'     \item `ladder_channel`: string, which channel in the fsa file contains the ladder signal. Default: `"DATA.105"`.
#'     \item `signal_channel`: string, which channel in the fsa file contains the data signal. Default: `"DATA.1"`.
#'     \item `ladder_sizes`: numeric vector, bp sizes of ladder used in fragment analysis. Default: `c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)`.
#'     \item `ladder_start_scan`: single numeric indicating the scan number to start looking for ladder peaks (only required when ladder signal does not have large spike at start). Usually this can be automatically found (when set to NA) through the detection of the large spike at the start of the signal. Default: `NA`.
#'     \item `minimum_ladder_signal`: single numeric for minimum signal of peak from smoothed signal. Default: `NA`.
#'     \item `ladder_assign_left_to_right`: single logical for if the ladder should be assigned from the smallest base pair size to largest (TRUE), or if the order should be reversed and assigned from largest to smallest (FALSE), which can be helpful since the end often has cleaner signal than the start. Default: `TRUE`.
#'     \item `ladder_selection_window`: single numeric for the ladder assigning algorithm. We iterate through the scans in blocks and test their linear fit (We can assume that the ladder is linear over a short distance). This value defines how large that block of peaks should be. Larger values should be better because the fit is tested in greater context, but larger numbers will make the fit increasingly slower. Default: `5`.
#'     \item `ladder_top_n_branching`: single numeric. The ladder assigning algorithm branches as it tests the various combinations. This value defines how many branches should be created. If the correct combination is not found, you could try increasing this value, but it will make it increasingly slower. Default: `5`.
#'     \item `ladder_branching_r_squared_threshold`: single numeric. The branches of the ladder assigning algorithm are pruned by R-squared values above this threshold to discard fits that are not promising. If the correct combination is not found, you could try decreasing this value, but it will make it increasingly slower. Default: `0.99`. 
#'     \item `min_scan`: single numeric indicating the lower scan limit to filter out scans below. Default: `NA`.
#'     \item `max_scan`: single numeric indicating the upper scan limit to filter out scans above Default: `NA`.
#'     \item `max_combinations`: single numeric indicating what is the maximum number of ladder combinations that should be tested. Default: `2500000`.
#'     \item `warning_rsq_threshold`: single numeric for the value for which this function will warn you when parts of the ladder have R-squared values below the specified threshold. Default: `0.998`.
#'     \item `show_progress_bar`: single logical for showing progress bar. Default: `TRUE`.
#'   }
#'
#' @return This function modifies list of fragments objects in place with the ladder assigned and base pair calculated.
#' @keywords internal
#'
#' @details
#' This function takes a list of fragments files (the output from read_fsa) and identifies
#' the ladders in the ladder channel which is used to call the bp size. The output
#' is a list of fragments. 
#' 
#' In this package, base pair (bp) sizes are assigned using a generalized additive model (GAM) with cubic regression splines. The model is fit to known ladder fragment sizes and their corresponding scan positions, capturing the relationship between scan number and bp size. Once trained, the model predicts bp sizes for all scans by interpolating between the known ladder points. This approach provides a flexible and accurate assignment of bp sizes, accommodating the slightly non-linear relationship.
#'
#' Use [plot_data_channels()] to plot the raw data on the fsa file to identify which channel the ladder and data are in.
#'
#'
#' Each ladder should be manually inspected to make sure that is has been correctly assigned.
#'
#' @seealso [plot_data_channels()] to plot the raw data in all channels. [plot_ladders()] to plot the assigned ladder
#' peaks onto the raw ladder signal. [fix_ladders_interactive()] to fix ladders with
#' incorrectly assigned peaks.
#' 
#' @importFrom mgcv gam
#'
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#' config <- load_config()
#'
#' trace:::find_ladders(fsa_list, config, show_progress_bar = FALSE)
#'
#' # Manually inspect the ladders
#' plot_ladders(fsa_list[1])
#'
find_ladders <- function(
  fragments_list,
  config,
  ...) {

  fit_ladder <- function(
      ladder,
      scans,
      sample_id) {
    
    ladder_df <- data.frame(signal = ladder, scan = scans)
    ladder_df <- ladder_df[which(ladder_df$scan >= ladder_start_scan), ]
    ladder_df$smoothed_signal <- pracma::savgol(
      ladder_df$signal,
      21
    )

    ladder_peaks <- find_ladder_peaks(
      ladder_df = ladder_df,
      n_reference_sizes = length(config$ladder_sizes),
      minimum_ladder_signal = config$minimum_ladder_signal,
      sample_id = sample_id
    )

    if(config$ladder_assign_left_to_right){
      reference_sizes = config$ladder_sizes
      observed_sizes = ladder_peaks
    }else{
      reference_sizes = rev(config$ladder_sizes)
      observed_sizes = rev(ladder_peaks)
    }

    if(config$ladder_selection_window >= length(reference_sizes)){
      # if testing all combinations ladder_iteration is slower
      peaks_fit_df <- exhaustive_ladder_matching(
        reference_sizes = reference_sizes, 
        observed_sizes = observed_sizes,
        max_combinations = config$max_combinations)
    } else{
      peaks_fit_df <- ladder_iteration(
        reference_sizes = reference_sizes, 
        observed_sizes = observed_sizes,
        choose = config$ladder_selection_window,
        max_combinations = config$max_combinations,
        top_n_branching = config$ladder_top_n_branching,
        r_squared_threshold = config$ladder_branching_r_squared_threshold
      )
    }

    peaks_not_fit <- ladder_peaks[which(!ladder_peaks %in% peaks_fit_df$scan)]

    peaks_not_fit_df <- data.frame(
      scan = peaks_not_fit,
      size = rep(NA_real_, length(peaks_not_fit))
    )


    combined_ladder_peaks <- rbind(peaks_fit_df, peaks_not_fit_df)
    combined_ladder_peaks <- combined_ladder_peaks[order(combined_ladder_peaks$scan), ]

    # add back in attr
    attr(combined_ladder_peaks, "total_combinations_tested") <- attr(peaks_fit_df, "total_combinations_tested")

    return(combined_ladder_peaks)
  }

  # prepare output
  output <- trace_output$new("find_ladders")

  config <- tryCatch(
    update_config(config, list(...)),
    error = function(e) e
  )
  if("error" %in% class(config)){
    output$set_status(
      "error", 
      config$message
    )
    return(output)
  }

  if (config$show_progress_bar) {
    pb <- utils::txtProgressBar(min = 0, max = length(fragments_list), style = 3)
  }

  for (i in seq_along(fragments_list)) {
    # populate the ladder and data channels with the supplied channel name
    ## check first to make sure that name is actually in the object

    for (channel in c(config$ladder_channel, config$signal_channel)) {
      if(!channel %in% names(fragments_list[[i]]$fsa$Data)){
        output$set_status(
          "error", 
          paste0(channel, " not detected as a channel in fsa")
        )
        return(output)
      }
    }

    fragments_list[[i]]$raw_ladder <- fragments_list[[i]]$fsa$Data[[config$ladder_channel]]
    fragments_list[[i]]$raw_data <- fragments_list[[i]]$fsa$Data[[config$signal_channel]]
    fragments_list[[i]]$scan <- 0:(length(fragments_list[[i]]$fsa$Data[[config$signal_channel]]) - 1)
    fragments_list[[i]]$off_scale_scans <- fragments_list[[i]]$fsa$Data$OfSc.1

    # set ladder_start_scan
    if (is.na(config$ladder_start_scan)) {
      ladder_start_scan <- which.max(fragments_list[[i]]$raw_ladder) + 50
    } else{
      ladder_start_scan <- config$ladder_start_scan
    }

    # allow user to subset to particular scans
    if (!is.na(config$min_scan) | !is.na(config$max_scan)) {
      min_scan = ifelse(is.na(config$min_scan), min(fragments_list[[i]]$scan), config$min_scan)
      max_scan = ifelse(is.na(config$max_scan), max(fragments_list[[i]]$scan), config$max_scan)

      fragments_list[[i]]$raw_ladder <- fragments_list[[i]]$raw_ladder[min_scan:max_scan]
      fragments_list[[i]]$raw_data <- fragments_list[[i]]$raw_data[min_scan:max_scan]
      fragments_list[[i]]$scan <- fragments_list[[i]]$scan[min_scan:max_scan]

      # set spike location since it's automatically set usually, and user may select scans to start after
      if(is.na(config$ladder_start_scan)){
        ladder_start_scan <- config$min_scan
      }
      
    }

    # ladder
    ladder_df <- tryCatch(
      fit_ladder(
        ladder = fragments_list[[i]]$raw_ladder,
        scans = fragments_list[[i]]$scan,
        sample_id = fragments_list[[i]]$unique_id
      ),
      error = function(e) e
    )
    if("error" %in% class(ladder_df)){
      output$set_status(
        "error", 
        paste0("There was an issue fitting the ladder for ", fragments_list[[i]]$unique_id,":\n",
                ladder_df$message
      )
      )
      return(output)
    }

    fragments_list[[i]]$ladder_df <- ladder_df
    if(!is.null(ladder_df) || nrow(ladder_df) > 0){
      fragments_list[[i]]$ladder_total_combinations_tested <- attr(ladder_df, "total_combinations_tested")
    }

    # ladder correlation stats
    # make a warning if one of the ladder modes is bad
    rsq_warning <- tryCatch(
      ladder_rsq_warning_helper(fragments_list[[i]],
        rsq_threshold = config$warning_rsq_threshold
      ),
      warning = function(w) w
    )
    if("warning" %in% class(rsq_warning)){
      output$set_status(
        "warning", 
        rsq_warning$message
      )
    }

    predicted_size <- predict_bp_size(
      ladder_df = ladder_df,
      scans = fragments_list[[i]]$scan
    )

    fragments_list[[i]]$trace_bp_df <- data.frame(
      unique_id = rep(fragments_list[[i]]$unique_id, length(fragments_list[[i]]$scan)),
      scan = fragments_list[[i]]$scan,
      size = predicted_size,
      signal = fragments_list[[i]]$raw_data,
      ladder_signal = fragments_list[[i]]$raw_ladder,
      off_scale = fragments_list[[i]]$scan %in% fragments_list[[i]]$off_scale_scans
    )

    if (config$show_progress_bar) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  # if too many messages, change to just say many samples have bad ladders
  if(output$status == "warning" && length(output$warning_message) > 10){
    output$warning_message = NULL
    output$set_status(
      "warning", 
      "Many samples appear to have badly fitting ladders based on the warning_rsq_threshold. Use extract_ladder_summary() to get a summary of the ladder fits."
    )
  }

  if (config$show_progress_bar) {
    # make sure progress bar ends on new line
    cat("\n")
  }
  
  return(output)
}



#' Fix ladders manually
#'
#' Manually assign the ladder peaks for samples in a fragments_list
#'
#' @param fragments_list list of fragments objects
#' @param ladder_df_list a list of dataframes, with the names being the unique id
#' and the value being a dataframe. The dataframe has two columns, size (indicating
#' the bp of the standard) and scan (the scan value of the ladder peak). It's
#' critical that the element name in the list is the unique id of the sample.
#' @param warning_rsq_threshold The value for which this function will warn you when parts of the ladder have R-squared values below the specified threshold.
#'
#' @return This function modifies list of fragments objects in place with the selected ladders fixed.
#' @keywords internal
#'
#' @details
#' This function returns a fragments list the same length as was supplied.
#' It goes through each sample and either just returns the same fragments
#' if the unique id doesn't match the samples that need the ladder fixed, or if
#' it is one to fix, it will use the supplied dataframe in the ladder_df_list
#' as the ladder. It then reruns the bp sizing methods on those samples.
#'
#' This is best used with [fix_ladders_interactive()] that can generate a `ladder_df_list`.
#'
#'
#' @examples
#' 
#' config <- load_config()
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' trace:::find_ladders(fsa_list, config, show_progress_bar = FALSE)
#'
#' # first manually determine the real ladder peaks using your judgment
#' # the raw ladder signal can be extracted
#' raw_ladder <- fsa_list[1]$raw_ladder
#'
#' # or we can look at the "trace_bp_df" to see a dataframe that includes the scan and ladder signal
#' raw_ladder_df <- fsa_list[[1]]$trace_bp_df[, c("unique_id", "scan", "ladder_signal")]
#' plot(raw_ladder_df$scan, raw_ladder_df$ladder_signal)
#'
#' # once you have figured what sizes align with which peak, make a dataframe. The
#' # fix_ladders_manual() function takes a list as an input so that multiple ladders
#' # can be fixed. Each sample would have the the list element name as it's unique id.
#'
#' example_list <- list(
#'  "20230413_A07.fsa" = data.frame(
#'    size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#'    scan = c(1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)
#'  )
#' )
#'
#' trace:::fix_ladders_manual(
#'   fsa_list,
#'   example_list
#' )
#'
fix_ladders_manual <- function(fragments_list,
                               ladder_df_list,
                               warning_rsq_threshold = 0.998) {
  # prepare output file
  output <- trace_output$new("fix_ladders_manual")

  
  samples_to_fix <- names(ladder_df_list)
  for (i in seq_along(fragments_list)) {
    if (fragments_list[[i]]$unique_id %in% samples_to_fix) {
      message(paste("Fixing ladder for", fragments_list[[i]]$unique_id))

      tmp_ladder_df <- ladder_df_list[[which(names(ladder_df_list) == fragments_list[[i]]$unique_id)]]

      # do some quality control of the df user supplied
      for (col in c("scan", "size")) {
        if(!col %in% colnames(tmp_ladder_df)){
          output$set_status(
            "error", 
            "Dataframe must contain column names 'size' and 'scan'"
          )
          return(output)
        }
      }

      fragments_list[[i]]$ladder_df <- tmp_ladder_df
      predicted_size <- predict_bp_size(
        fragments_list[[i]]$ladder_df,
        fragments_list[[i]]$scan
      )
    
      fragments_list[[i]]$trace_bp_df <- data.frame(
        unique_id = rep(fragments_list[[i]]$unique_id, length(fragments_list[[i]]$scan)),
        scan = fragments_list[[i]]$scan,
        size = predicted_size,
        signal = fragments_list[[i]]$raw_data,
        ladder_signal = fragments_list[[i]]$raw_ladder,
        off_scale = fragments_list[[i]]$scan %in% fragments_list[[i]]$off_scale_scans
      )
    
      # make a warning if one of the ladder modes is bad
      rsq_warning <- tryCatch(
        ladder_rsq_warning_helper(fragments_list[[i]],
          rsq_threshold = warning_rsq_threshold
        ),
        warning = function(w) w
      )
      if("warning" %in% class(rsq_warning)){
        output$set_status(
          "warning", 
          rsq_warning$message
        )
      }
    }
  }

  return(output)
}
