######################## Helper functions ####################################


deshoulder <- function(peak_table_df, shoulder_window) {
  fragment_size <- peak_table_df$size
  heights <- peak_table_df$height
  fragment_size_deshoulder <- numeric()
  for (i in seq_along(fragment_size)) {
    if (i == 1) {
      if (fragment_size[i + 1] - fragment_size[i] > shoulder_window) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else if (heights[i] > heights[i + 1]) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else {
        next
      }
    } else if (i == length(fragment_size)) {
      if (fragment_size[i - 1] - fragment_size[i] > shoulder_window) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else if (heights[i] > heights[i - 1]) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else {
        next
      }
    } else if (fragment_size[i + 1] - fragment_size[i] < shoulder_window | fragment_size[i] - fragment_size[i - 1] < shoulder_window) {
      before <- fragment_size[i + 1] - fragment_size[i] < shoulder_window
      after <- fragment_size[i] - fragment_size[i - 1] < shoulder_window
      before_and_higher <- heights[i] > heights[i + 1]
      after_and_higher <- heights[i] > heights[i - 1]

      if (before & after) {
        if (before_and_higher & after_and_higher) {
          fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
        } else {
          next
        }
      } else if (before && before_and_higher) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else if (after && after_and_higher) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else {
        next
      }
    } else {
      fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
    }
  }

  new_peak_table_df <- peak_table_df[which(peak_table_df$size %in% fragment_size_deshoulder), ]

  return(new_peak_table_df)
}






# force whole repeat unit algorithm -----------------------------------------------------------
np_repeat <- function(size,
                      main_peak_size,
                      main_peak_repeat,
                      repeat_size) {
  # For loop to get values
  # first find the distance to main peak and the peak before
  size_delta_from_main_peak <- size - main_peak_size
  # create a numeric vector of length of peak to store values
  np_repeat <- vector(mode = "numeric", length = length(size))
  # note that the loop has to start in the middle at the main peak since it's looking for the previous value in the new vector
  for (i in seq_along(size)) {
    # set main peak
    if (i == which(size_delta_from_main_peak == 0)) {
      np_repeat[[i]] <- main_peak_repeat
    }
    # calculate for peaks greater than main peak
    if (i > which(size_delta_from_main_peak == 0)) {
      # calculate size distance to nearest peak in cag length,
      # add that on to previous cag, then round to whole cag
      np_repeat[[i]] <-
        np_repeat[[i - 1]] + round((size[[i]] - size[[i - 1]]) / repeat_size)
    }
  }
  # --- reverse order and find smaller repeats ---
  size_rev <- rev(size)
  np_repeat_rev <- rev(np_repeat)
  size_delta_from_main_peak_rev <- rev(size_delta_from_main_peak)
  # second loop that goes along the vector in the from larger to smaller bp peaks
  for (i in seq_along(size_rev)) {
    # calculate peaks greater than main peak (smaller bp peaks since reversed)
    if (i > which(size_delta_from_main_peak_rev == 0)) {
      # note this is a subtraction since the size_delta_peak_after is a negative value
      np_repeat_rev[[i]] <-
        np_repeat_rev[[i - 1]] + round((size_rev[[i]] - size_rev[[i - 1]]) / repeat_size)
    }
  }
  return(rev(np_repeat_rev))
}


# fft algo ----------------------------------------------------------------


####### fft helpers
fragments_fft <- function(df) {
  Fs <- 1 # Arbitrary sampling frequency
  Ts <- 1 / Fs # Arbitrary sampling time

  signal_periodic <- as.vector(pracma::detrend(df$signal))

  n <- length(signal_periodic)
  xf <- seq(0, 1 / (2 * Ts), length.out = n %/% 2) # int div in R is % / %
  yf <- stats::fft(signal_periodic)
  yf <- (2 / n) * abs(yf[seq(1, n %/% 2)])

  data.frame(
    "xf" = xf,
    "yf" = yf
  )
}


find_main_freq <- function(fft_df, skip_rows = 3) {
  # skip rows to avoid noise at start
  df2 <- fft_df[skip_rows + 1:nrow(fft_df), ]
  peaks <- pracma::findpeaks(df2$yf)
  fund_freq_position <- peaks[which(peaks[, 1] == max(peaks[, 1], na.rm = TRUE)), 2][1]
  freq <- df2[fund_freq_position, "xf"] * 3

  return(freq)
}

find_scan_period <- function(df,
                             main_peak_scan) {
  fft_df <- fragments_fft(df)

  # detrend signal
  fft_df$yf <- detrend_signal(fft_df$yf)

  main_freq <- find_main_freq(fft_df)
  pure_wave <- cos(2 * main_freq * (df$scan - main_peak_scan))
  cos_max <- pracma::findpeaks(pure_wave, nups = 3)
  scans_diffs <- diff(df[cos_max[, 2], "scan"])
  peak_scan_peroid <- round(median(scans_diffs))


  return(peak_scan_peroid)
}

find_peaks_by_scan_period <- function(df,
                                      main_peak_scan,
                                      peak_scan_period,
                                      direction,
                                      window) {
  if (direction == 1) {
    df_post_main <- df[which(df$scan > main_peak_scan), ]
  } else {
    df_post_main <- df[which(df$scan < main_peak_scan), ]
    df_post_main <- df_post_main[order(df_post_main$scan, decreasing = TRUE), ]
    peak_scan_period <- peak_scan_period * -1
  }

  called_peaks <- numeric()
  current_scan_position <- main_peak_scan + peak_scan_period
  while (TRUE) {
    window_range <- (current_scan_position - window):(current_scan_position + window)
    window_df <- df_post_main[df_post_main$scan %in% window_range, ]

    if (nrow(window_df) > 0) {
      tallest_in_window <- window_df[which.max(window_df$signal), "scan"]
      called_peaks <- c(called_peaks, tallest_in_window)

      # Update current scan position
      current_scan_position <- tallest_in_window + peak_scan_period
    } else {
      # If no more data points, terminate
      break
    }
  }

  return(called_peaks)
}


fft_repeat_caller <- function(fragments_repeat,
                              scan_peak_window = 3,
                              fragment_window = 3 * 5) {
  if (is.na(fragments_repeat$get_alleles()$allele_1_size)) {
    df <- data.frame(
      "unique_id" = character(),
      "scan" = numeric(),
      "size" = numeric(),
      "signal" = numeric()
    )

    return(df)
  }

  fragment_window_positions <- which(fragments_repeat$trace_bp_df$size > fragments_repeat$get_alleles()$allele_1_size - fragment_window & fragments_repeat$trace_bp_df$size < fragments_repeat$get_alleles()$allele_1_size + fragment_window)
  window_df <- fragments_repeat$trace_bp_df[fragment_window_positions, ]
  main_peak_scan <- window_df[which(window_df$size == fragments_repeat$get_alleles()$allele_1_size), "scan"]

  peak_scan_peroid <- find_scan_period(window_df, main_peak_scan)

  pos_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_peroid,
    direction = 1,
    window = scan_peak_window
  )

  neg_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_peroid,
    direction = -1,
    window = scan_peak_window
  )

  peak_table <- fragments_repeat$trace_bp_df
  peak_table <- peak_table[which(peak_table$scan %in% c(neg_peaks, main_peak_scan, pos_peaks)), ]
  peak_table <- peak_table[which(peak_table$size > fragments_repeat$.__enclos_env__$private$min_bp_size & peak_table$size < fragments_repeat$.__enclos_env__$private$max_bp_size), ]


  return(peak_table)
}






size_period_repeat_caller <- function(fragments_repeat,
                                      size_period,
                                      scan_peak_window = 3,
                                      fragment_window = 3 * 5) {
  if (is.na(fragments_repeat$get_alleles()$allele_1_size)) {
    df <- data.frame(
      "unique_id" = character(),
      "scan" = numeric(),
      "size" = numeric(),
      "signal" = numeric()
    )

    return(df)
  }

  fragment_window_positions <- which(fragments_repeat$trace_bp_df$size > fragments_repeat$get_alleles()$allele_1_size - fragment_window & fragments_repeat$trace_bp_df$size < fragments_repeat$get_alleles()$allele_1_size + fragment_window)
  window_df <- fragments_repeat$trace_bp_df[fragment_window_positions, ]
  main_peak_scan <- window_df[which(window_df$size == fragments_repeat$get_alleles()$allele_1_size), "scan"]


  # determine period

  peak_scan_peroid <- round(size_period / median(diff(window_df$size)))

  pos_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_peroid,
    direction = 1,
    window = scan_peak_window
  )

  neg_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_peroid,
    direction = -1,
    window = scan_peak_window
  )

  peak_table <- fragments_repeat$trace_bp_df
  peak_table <- peak_table[which(peak_table$scan %in% c(neg_peaks, main_peak_scan, pos_peaks)), ]
  peak_table <- peak_table[which(peak_table$size > fragments_repeat$.__enclos_env__$private$min_bp_size & peak_table$size < fragments_repeat$.__enclos_env__$private$max_bp_size), ]


  return(peak_table)
}



# find_size_batch_correction_factor

find_batch_correction_factor <- function(fragments_list, trace_window_size = 50, smoothing_window = 301){
  # make df for all samples of plate id, batch sample id
  metadata_list <- lapply(fragments_list, function(x){
    df <- data.frame(
      unique_id = x$unique_id,
      batch_run_id = x$batch_run_id,
      batch_sample_id = x$batch_sample_id
    )
    return(df)
  })

  metadata_df <- do.call(rbind, metadata_list)

  # make sure there is a least one run with common samples to all others
  size_std_df <- metadata_df[which(!is.na(metadata_df$batch_sample_id)), , drop = FALSE]
  if(nrow(size_std_df) == 0 ){
    stop(call. = FALSE, "Size correction requires 'batch_sample_id' samples indicated in the metadata. See ?add_metadata for more info.")
  }


  sample_run_table <- table(metadata_df$batch_sample_id, metadata_df$batch_run_id)
  runs_with_all_samples <- colnames(sample_run_table)[colSums(sample_run_table == 1) == length(unique(na.omit(metadata_df$batch_sample_id)))]

  if(length(runs_with_all_samples) == 0){
    stop(call. = FALSE, "There needs to be at least one batch_run_id with all batch_sample_id for batch correction.")
  }
  index_batch_run_id <- runs_with_all_samples[1]

  # first group samples across runs by sample id
  size_std_fragments <- fragments_list[size_std_df$unique_id]
  split_by_sample_id <- split(size_std_fragments, sapply(size_std_fragments, function(x) x$batch_sample_id ))

  correction_factor_by_sample_list <- lapply(split_by_sample_id, function(sample_list){
    sample_modal_size <- sapply(sample_list, function(fragment){
      df <- fragment$trace_bp_df[which(fragment$trace_bp_df$size < fragment$get_alleles()$allele_1_size + trace_window_size & fragment$trace_bp_df$size > fragment$get_alleles()$allele_1_size - trace_window_size ), ]
      df$smoothed <- pracma::savgol(df$signal, smoothing_window)
      smoothed_modal_size <- df[which.max(df$smoothed), "size"]
      return(smoothed_modal_size)
    })

    df <- data.frame(
      unique_id = sapply(sample_list, function(x) x$unique_id),
      batch_run_id = sapply(sample_list, function(x) x$batch_run_id),
      smoothed_modal_size = sample_modal_size
    )

    index_plate_smoothed_modal_size <- median(df[which(df$batch_run_id == index_batch_run_id), "smoothed_modal_size"])
    df$correction_factor <- df$smoothed_modal_size - index_plate_smoothed_modal_size

    return(df)

  })
  correction_factor_by_sample_df <- do.call(rbind, correction_factor_by_sample_list)

  #now group by plate and find the average correction factor
  correction_factor_by_plate_df_list <- split(correction_factor_by_sample_df, correction_factor_by_sample_df$batch_run_id)
  correction_factor_by_plate_list <- lapply(correction_factor_by_plate_df_list, function(x) median(x$correction_factor))

  # save correction factor for each class object
  for (i in seq_along(fragments_list)) {
    # Made the plate id explicity match the list name for cases when the plate name is a number. It could cause subsetting issues
    fragments_list[[i]]$.__enclos_env__$private$batch_correction_factor <- correction_factor_by_plate_list[[which(names(correction_factor_by_plate_list) == fragments_list[[i]]$batch_run_id)]]
  }
  
}




# call_repeats ------------------------------------------------------------

#' Call Repeats for Fragments
#'
#' This function calls the repeat lengths for a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#' @param assay_size_without_repeat An integer specifying the assay size without repeat for repeat calling. Default is 87.
#' @param repeat_size An integer specifying the repeat size for repeat calling. Default is 3.
#' @param force_whole_repeat_units A logical value specifying if the peaks should be forced to be whole repeat units apart. Usually the peaks are slightly under the whole repeat unit if left unchanged.
#' @param batch_correction A logical specifying if the size should be adjusted across fragment analysis runs. Requires metadata to be added to specify samples (\code{"batch_sample_id"}) common across runs (\code{"batch_run_id"})(see [add_metadata()]).
#' @param repeat_calling_algorithm A character specifying the repeat calling algorithm. Options: \code{"simple"}, \code{"fft"}, or \code{"size_period"} (see details section for more information on these).
#' @param repeat_calling_algorithm_size_window_around_allele A numeric value for how big of a window around the tallest peak should be used to find the peak periodicity. Used for both \code{"fft"} and \code{"size_period"}. For \code{"fft"}, you want to make sure that this window is limited to where there are clear peaks. For \code{"size_period"}, it will not make a big difference.
#' @param repeat_calling_algorithm_peak_assignment_scan_window A numeric value for the scan window when assigning the peak. This is used for both \code{"fft"} and \code{"size_period"}. When the scan period is determined, the algorithm jumps to the predicted scan for the next peak. This value opens a window of the neighboring scans to pick the tallest in.
#' @param repeat_calling_algorithm_size_period A numeric value \code{"size_period"} algorithm to set the peak periodicity by bp size. This is the key variable to change for \code{"size_period"}. In fragment analysis, the peaks are usually slightly below the actual repeat unit size.
#'
#' @return A list of \code{"fragments_repeats"} objects with repeat data added.
#'
#' @details
#' The calculated repeat lengths are assigned to the corresponding peaks in the provided `fragments_repeats` object. The repeat lengths can be used for downstream instability analysis.
#'
#' The `simple` algorithm is just the repeat size calculated directly from bp: (bp - assay_size_without_repeat) / repeat_size
#'
#' The `fft` or `size_period` algorithms both re-call the peaks based on empirically determined (`fft`) or specified (`size_period`) periodicity of the peaks. The main application of these algorithms is to solve the issue of contaminating peaks that are not expected in the expected regular pattern of peaks. The `fft` approach applies a fourier transform to the peak signal to determine the underlying periodicity of the signal. `size_period` is similar and simpler, where instead of automatically figuring out the periodicity we as users usually know the size distance between repeat units. We can use that known peroidicty to jump between peaks.
#'
#' The `force_whole_repeat_units` algorithm aims to correct for the systematic drift in fragment sizes that occurs. It calculates repeat lengths in a way that helps align peaks with the underlying repeat pattern, making the estimation of repeat lengths more reliable relative to the main peak. The calculated repeat lengths start from the main peak's repeat length and increases in increments of the specified `repeat_size`.
#'
#' @seealso [find_alleles()], [add_metadata()]
#'
#' @export
#'
#' @examples
#'
#' file_list <- trace::cell_line_fsa_list[c(90:92)]
#'
#' test_ladders <- find_ladders(file_list)
#'
#' fragments_list <- find_fragments(test_ladders,
#'   min_bp_size = 300
#' )
#'
#' test_alleles <- find_alleles(
#'   fragments_list = fragments_list
#' )
#'
#' # Simple conversion from bp size to repeat size
#' test_repeats <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats[1], xlim = c(120, 170))
#'
#'
#' # use different algorithms to call the repeats to ensure only periodic peaks are called
#'
#' # fft to automatically find peak period
#' test_repeats_fft <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "fft",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats_fft[1], xlim = c(120, 170))
#'
#' # size_period to manually supply the peak period
#' test_repeats_size_period <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "size_period",
#'   repeat_calling_algorithm_size_period = 2.75,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats_size_period[1], xlim = c(120, 170))
#'
#'
#' # Use force_whole_repeat_units algorithm to make sure called
#' # repeats are the exact number of bp apart
#'
#' test_repeats_whole_units <- call_repeats(
#'   fragments_list = test_alleles,
#'   force_whole_repeat_units = TRUE,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats_whole_units[1], xlim = c(120, 170))
#'
#' # correct repeat length from metadata
#'
#' test_alleles_metadata <- add_metadata(
#'   fragments_list = test_alleles,
#'   metadata_data.frame = trace::metadata
#' )
#'
#' test_repeats_corrected <- call_repeats(
#'   fragments_list = test_alleles_metadata,
#'   batch_correction = TRUE
#' )
#'
#'
#' plot_traces(test_repeats_corrected[1], xlim = c(120, 170))
#'
call_repeats <- function(
    fragments_list,
    assay_size_without_repeat = 87,
    repeat_size = 3,
    force_whole_repeat_units = FALSE,
    batch_correction = FALSE,
    repeat_calling_algorithm = "simple",
    repeat_calling_algorithm_size_window_around_allele = repeat_size * 5,
    repeat_calling_algorithm_peak_assignment_scan_window = 3,
    repeat_calling_algorithm_size_period = repeat_size * 0.93) {
  # Check to see if repeats are to be batch corrected
  # if so, find correction factor by looking across all samples before drilling down to a per sample level
  # the size correction then needs to be applied downstream
  if (batch_correction) {
    find_batch_correction_factor(fragments_list)
  }
  # call repeats for each sample
  added_repeats <- lapply(
    fragments_list,
    function(fragment) {
      ### in this function, we are doing three key things
      #### 1) correct bp size if required
      #### 2) calculate repeats by assay_size_without_repeat and repeat size
      #### 3) use a method to calculate repeats

      # check to make sure all the required inputs for the function have been given
      if (fragment$.__enclos_env__$private$find_main_peaks_used == FALSE) {
        stop(paste0(fragment$unique_id, " requires main alleles to be identified before repeats can be called. Find alleles using 'find_main_peaks()' whitin the class, or use the 'find_alleles()' accesesor to find the main peaks across a list of 'fragments_repeats' objects"),
          call. = FALSE
        )
      }

      # only continue from here if main peaks were successfully found, otherwise, don't return repeat data (ie it can be an empty df)
      if (is.na(fragment$get_alleles()$allele_1_size) | is.na(fragment$get_alleles()$allele_1_height)) {
        fragment$.__enclos_env__$private$repeats_not_called_reason <- "No main peaks"
        warning(paste0(fragment$unique_id, ": repeats were not called (no main peaks in sample)"),
          call. = FALSE
        )
        # populate with empty dataframe to help the rest of the pipeline
        fragment$repeat_table_df <- data.frame(
          unique_id = character(),
          size = numeric(),
          height = numeric(),
          repeats = numeric(),
          off_scale = logical()
        )
        # exit lapply early 
        return(fragment)
      } 
      # repeat calling algorithm
      if (repeat_calling_algorithm == "simple") {

        if (batch_correction) {
          peak_table_size <- fragment$peak_table_df$size - fragment$.__enclos_env__$private$batch_correction_factor
        } else{
          peak_table_size <- fragment$peak_table_df$size
        }

        repeat_table_df <- data.frame(
          unique_id = fragment$peak_table_df$unique_id,
          size = peak_table_size,
          height = fragment$peak_table_df$height,
          calculated_repeats = (peak_table_size- assay_size_without_repeat) / repeat_size,
          repeats = (peak_table_size - assay_size_without_repeat) / repeat_size,
          off_scale = ifelse(any(colnames(fragment$peak_table_df) == "off_scale"),
          fragment$peak_table_df$off_scale,
            rep(FALSE, nrow(fragment$peak_table_df))
          )
        )
      } else if (repeat_calling_algorithm == "fft") {
        # check to see that fragments repeats has trace data since that is required.
        if (is.null(fragment$trace_bp_df)) {
          stop("fft algorithim requires trace data. Use fsa samples rather than peak table is inputs into the pipeline.",
            call. = FALSE
          )
        }

        fft_peak_df <- fft_repeat_caller(fragment,
          fragment_window = repeat_calling_algorithm_size_window_around_allele,
          scan_peak_window = repeat_calling_algorithm_peak_assignment_scan_window
        )

        if (batch_correction) {
          fft_peak_size <- fft_peak_df$size - fragment$.__enclos_env__$private$batch_correction_factor
        } else{
          fft_peak_size <- fft_peak_df$size
        }

        repeat_table_df <- data.frame(
          unique_id = fft_peak_df$unique_id,
          size = fft_peak_size,
          height = fft_peak_df$signal,
          calculated_repeats = (fft_peak_size - assay_size_without_repeat) / repeat_size,
          repeats = (fft_peak_size - assay_size_without_repeat) / repeat_size,
          off_scale = fft_peak_df$off_scale
        )
      } else if (repeat_calling_algorithm == "size_period") {
        # check to see that fragments repeats has trace data since that is required.
        if (is.null(fragment$trace_bp_df)) {
          stop("size_period algorithim requires trace data. Use fsa samples rather than peak table is inputs into the pipeline.",
            call. = FALSE
          )
        }

        size_period_df <- size_period_repeat_caller(fragment,
          size_period = repeat_calling_algorithm_size_period,
          fragment_window = repeat_calling_algorithm_size_window_around_allele,
          scan_peak_window = repeat_calling_algorithm_peak_assignment_scan_window
        )

        if (batch_correction) {
          size_period_size <- size_period_df$size - fragment$.__enclos_env__$private$batch_correction_factor
        } else{
          size_period_size <- size_period_df$size
        }

        repeat_table_df <- data.frame(
          unique_id = size_period_df$unique_id,
          size = size_period_size,
          height = size_period_df$signal,
          calculated_repeats = (size_period_size - assay_size_without_repeat) / repeat_size,
          repeats = (size_period_size - assay_size_without_repeat) / repeat_size,
          off_scale = size_period_df$off_scale
        )
      } else {
        stop(
          call. = FALSE,
          "Invalid repeat calling algorithim selected"
        )
      }

      # Force the repeat units to be whole numbers
      if (force_whole_repeat_units == TRUE) {
        repeat_table_df$repeats <- np_repeat(
          size = repeat_table_df$size,
          main_peak_size = fragment$get_alleles()$allele_1_size,
          main_peak_repeat = repeat_table_df$calculated_repeats[which(repeat_table_df$size == fragment$get_alleles()$allele_1_size)],
          repeat_size = repeat_size
        )
      }

      # Finally save main peak repeat length and repeats data
      fragment$repeat_table_df <- repeat_table_df
      allele_1_subset <- repeat_table_df$repeats[which(repeat_table_df$size == fragment$get_alleles()$allele_1_size)]
      fragment$set_allele(allele = 1, unit = "repeats", value = ifelse(length(allele_1_subset) == 1, allele_1_subset, NA_real_))
      allele_2_subset <- repeat_table_df$repeats[which(repeat_table_df$size == fragment$get_alleles()$allele_2_size)]
      fragment$set_allele(allele = 2, unit = "repeats", value = ifelse(length(allele_2_subset) == 1, allele_2_subset, NA_real_))

      # also calculate repeat length for the trace-level data if it exists
      if (!is.null(fragment$trace_bp_df)) {
        if (batch_correction) {
          trace_bp_size <- fragment$trace_bp_df$size - fragment$.__enclos_env__$private$batch_correction_factor
        } else{
          trace_bp_size <- fragment$trace_bp_df$size
        }
        fragment$trace_bp_df$calculated_repeats <- (trace_bp_size - assay_size_without_repeat) / repeat_size
      }

      return(fragment)
    }
  )

  return(added_repeats)
}
