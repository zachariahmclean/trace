# General helper functions --------------------------------------------------------

detrend_signal <- function(x, bins = 50) {
  x_split <- split(x, ceiling(seq_along(x) / (length(x) / bins)))
  median_signal <- sapply(x_split, median)

  x_detrended_list <- vector("list", length = length(x_split))
  for (i in seq_along(x_detrended_list)) {
    x_detrended_list[[i]] <- x_split[[i]] - median_signal[i]
  }
  x_detrended <- unlist(x_detrended_list, use.names = FALSE)

  return(x_detrended)
}

process_ladder_signal <- function(ladder,
                                  scans,
                                  ladder_start_scan,
                                  smoothing_window) {

  ladder_df <- data.frame(signal = ladder, scan = scans)
  ladder_df <- ladder_df[which(ladder_df$scan >= ladder_start_scan), ]
  ladder_df$detrended_signal <- detrend_signal(ladder_df$signal)
  ladder_df$smoothed_signal <- pracma::savgol(
    ladder_df$detrended_signal,
    smoothing_window
  )
  return(ladder_df)
}

find_ladder_peaks <- function(ladder_df,
                              n_reference_sizes,
                              minimum_peak_signal,
                              sample_id) {
  median_signal <- median(ladder_df$smoothed_signal, na.rm = TRUE)
  sd_signal <- stats::sd(ladder_df$smoothed_signal, na.rm = TRUE)

  ladder_peaks <- vector("numeric")
  ladder_peak_threshold <- 1

  #allow user to set min height

  if(!is.null(minimum_peak_signal)){
    peaks <- pracma::findpeaks(ladder_df$smoothed_signal,
                               peakpat = "[+]{6,}[0]*[-]{6,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                               minpeakheight = minimum_peak_signal
    )

    ladder_peaks <- ladder_df$scan[peaks[, 2]]
  } else {
    while (length(ladder_peaks) < n_reference_sizes) {
      peaks <- pracma::findpeaks(ladder_df$smoothed_signal,
                                 peakpat = "[+]{6,}[0]*[-]{6,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                                 minpeakheight = median_signal + sd_signal * ladder_peak_threshold
      )

      ladder_peaks <- ladder_df$scan[peaks[, 2]]

      # lower the threshold for the next cycle
      ladder_peak_threshold <- ladder_peak_threshold - 0.01

      # provide an exit if there are not enough peaks found
      if (sd_signal * ladder_peak_threshold <= 0) {
        break
      }
    }
  }

  if(length(ladder_peaks) < n_reference_sizes){
    stop(call. = FALSE,
         paste0("Fewer ladder peaks than reference ladder sizes were identified for ",
                sample_id,
                ". Adjust settings/ladder sizes to ensure the expected number of peaks are found."))
  }

  # go through raw signal and make sure that the identified scan in the smoothed signal is still the highest
  # it will also deal with cases where the scans have the same height (which.max will chose first)
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

ladder_iteration <- function(reference_sizes,
                             observed_sizes,
                             choose,
                             max_combinations) {
  find_best_combination <- function(recombinations,
                                    reference_sizes) {
    rsq_vector <- vector("numeric", ncol(recombinations))
    for (i in 1:ncol(recombinations)) {
      rsq_vector[[i]] <- stats::cor(reference_sizes, recombinations[, i])^2
    }
    # select the values from the original recombinations matrix not the one with the already selected values added for the regression
    selected_recombinations <- recombinations[, which.max(rsq_vector)]

    return(selected_recombinations)
  }

  assigned_reference <- vector("numeric")
  assigned_observed <- vector("numeric")
  max_iterations <- 1000
  iteration_count <- 1

  # Keep going until the number of reference peaks left is small
  while (length(reference_sizes) > 0 && iteration_count < max_iterations) {
    # observed_sizes and reference_sizes reset each loop

    # calculate the window to look at for this loop
    n_observed <- length(observed_sizes)
    n_reference <- length(reference_sizes)
    remainder <- n_reference - choose # how many peaks will be left to be assigned after this selection
    start_window <- n_observed - remainder


    if (n_observed == n_reference) {
      assigned_observed <- append(assigned_observed, observed_sizes)
      assigned_reference <- append(assigned_reference, reference_sizes)
      break
    } else {

      n_recombinations <- factorial(length(observed_sizes[1:start_window])) / (factorial(choose) * factorial(length(observed_sizes[1:start_window]) - choose))


      if (n_recombinations > max_combinations) {
        stop(
          call. = FALSE,
          paste0("Too many combinations to test (", n_recombinations, ")")
        )
      }

      recombinations <- utils::combn(observed_sizes[1:start_window], choose)

      selected_recombinations <- find_best_combination(
        recombinations = recombinations,
        reference_sizes = reference_sizes[1:choose]
      )
    }

    # assign the selections
    assigned_observed <- append(assigned_observed, selected_recombinations)
    assigned_reference <- append(assigned_reference, reference_sizes[1:choose])

    # set up up the next iteration
    last_selected_scan <- selected_recombinations[length(selected_recombinations)]
    last_selected_scan_position <- which(observed_sizes == last_selected_scan)

    last_selected_reference <- reference_sizes[choose]
    last_selected_reference_position <- which(reference_sizes == last_selected_reference)

    reference_sizes_left <- length(reference_sizes) - last_selected_reference_position

    # fix introduction of NAs by choose being longer than remaining sizes
    if (reference_sizes_left < choose) {
      choose <- reference_sizes_left
    }

    # break for final iteration
    if (choose == 0) {
      break
    }

    # deal with situation of just a couple of reference peaks left over, but you need a good amount to make correlation
    # therefore increase chose to the max
    if (reference_sizes_left - choose < 4) {
      choose <- reference_sizes_left
    } else if (reference_sizes_left < choose) { # fix introduction of NAs by choose being longer than remaining sizes
      choose <- reference_sizes_left
    }

    reference_sizes <- reference_sizes[(last_selected_reference_position + 1):length(reference_sizes)]
    observed_sizes <- observed_sizes[(last_selected_scan_position + 1):length(observed_sizes)]

    iteration_count <- iteration_count + 1
  }

  return(data.frame(scan = assigned_observed, size = assigned_reference))
}


ladder_rsq_warning_helper <- function(
    fragments_trace,
    rsq_threshold) {
  rsq <- sapply(fragments_trace$local_southern_mod, function(x) suppressWarnings(summary(x$mod)$r.squared))
  if (any(rsq < rsq_threshold)) {
    size_ranges <- sapply(fragments_trace$local_southern_mod, function(x) x$mod$model$yi)
    size_ranges <- size_ranges[, which(rsq < rsq_threshold), drop = FALSE]
    size_ranges_vector <- vector("numeric", ncol(size_ranges))
    for (j in seq_along(size_ranges_vector)) {
      size_ranges_vector[j] <- paste0(size_ranges[1, j], "-", size_ranges[3, j])
    }
    warning(
      call. = FALSE,
      paste(
        "sample", fragments_trace$unique_id, "has badly fitting ladder for bp sizes:",
        paste0(size_ranges_vector, collapse = ", ")
      )
    )
  }
}


# bp sizing ---------------------------------------------------------------


local_southern <- function(x, y) {
  # do some quality control. There should be no missing values and vectors should be same length
  if (length(x) != length(y)) {
    stop(
      call. = FALSE,
      "local_southern error: ladder scan and size vectors different lengths"
    )
  } else if (any(is.na(x)) | any(is.na(y))) {
    stop(
      call. = FALSE,
      "local_southern error: missing values in ladder scan or size"
    )
  }


  # Sort the data points by x values
  sorted_indices <- order(x)
  x_sorted <- x[sorted_indices]
  y_sorted <- y[sorted_indices]

  # Function to calculate the fitting constants for each group of three neighboring points
  mod_list <- vector("list", length = length(x_sorted) - 2)

  for (i in 1:(length(x_sorted) - 2)) {
    xi <- x_sorted[i:(i + 2)]
    yi <- y_sorted[i:(i + 2)]
    mod_list[[i]] <- list(
      mod = lm(yi ~ xi),
      first = xi[1],
      last = xi[3]
    )
  }

  return(mod_list)
}

local_southern_predict <- function(local_southern_output, scans) {
  # total number of groups to brake the scans into:
  ladder_scan_pos <- sapply(local_southern_output, function(fit) fit$first)

  # Find the nearest ladder position for each scan position
  nearest_ladder_index <- sapply(scans, function(scan) which.min(abs(scan - ladder_scan_pos)))

  # Assign the scan positions to corresponding groups based on nearest ladder position
  scan_split <- split(scans, nearest_ladder_index)
  size_split <- vector("list", length = length(scan_split))
  for (i in seq_along(scan_split)) {
    if (i == 1 | i == length(scan_split)) {
      size_split[[i]] <- stats::predict(local_southern_output[[i]]$mod, data.frame(xi = scan_split[[i]]))
    } else {
      lower_prediction <- stats::predict(local_southern_output[[i - 1]]$mod, data.frame(xi = scan_split[[i]]))
      upper_prediction <- stats::predict(local_southern_output[[i]]$mod, data.frame(xi = scan_split[[i]]))
      size_split[[i]] <- (lower_prediction + upper_prediction) / 2
    }
  }

  size <- unlist(size_split)

  return(size)
}


# ladder fixing -----------------------------------------------------------


ladder_fix_helper <- function(fragments_trace,
                              replacement_ladder_df) {

  fragments_trace$ladder_df <- replacement_ladder_df
  ladder_df <- fragments_trace$ladder_df[which(!is.na(fragments_trace$ladder_df$size)), ]
  ladder_df <- ladder_df[which(!is.na(ladder_df$scan)), ]
  fragments_trace$local_southern_mod <- local_southern(ladder_df$scan, ladder_df$size)

  predicted_size <- local_southern_predict(local_southern_output = fragments_trace$local_southern_mod, scans = fragments_trace$scan)

  fragments_trace$trace_bp_df <- data.frame(
    unique_id = rep(fragments_trace$unique_id, length(fragments_trace$scan)),
    scan = fragments_trace$scan,
    size = predicted_size,
    signal = fragments_trace$raw_data,
    ladder_signal = fragments_trace$raw_ladder,
    off_scale = fragments_trace$scan %in% fragments_trace$off_scale_scans
  )

  # make a warning if one of the ladder modes is bad
  ladder_rsq_warning_helper(fragments_trace,
    rsq_threshold = 0.998
  )

  return(fragments_trace)
}


ladder_self_mod_predict <- function(fragments_trace,
                                    size_threshold,
                                    size_tolerance,
                                    rsq_threshold) {

  ladder_sizes <- fragments_trace$ladder_df[which(!is.na(fragments_trace$ladder_df$size)), "size"]
  ladder_peaks <- fragments_trace$ladder_df$scan

  mod_validations <- vector("list", length(fragments_trace$local_southern_mod))
  for (i in seq_along(fragments_trace$local_southern_mod)) {
    predictions <- stats::predict(fragments_trace$local_southern_mod[[i]]$mod, newdata = data.frame(xi = ladder_peaks))
    low_size_threshold <- fragments_trace$local_southern_mod[[i]]$mod$model$yi[1] - size_threshold
    high_size_threshold <- fragments_trace$local_southern_mod[[i]]$mod$model$yi[3] + size_threshold

    predictions_close <- predictions[which(predictions > low_size_threshold & predictions < high_size_threshold)]
    mod_sizes <- fragments_trace$local_southern_mod[[i]]$mod$model$yi

    ladder_hits <- sapply(predictions_close, function(x) {
      diff <- ladder_sizes - x
      ladder_hit <- ladder_sizes[which(diff > -size_tolerance & diff < size_tolerance)]
      if (length(ladder_hit) == 0) {
        return(NA_real_)
      } else {
        return(ladder_hit)
      }
    })


    ladder_hits <- ladder_hits[!ladder_hits %in% mod_sizes & !is.na(ladder_hits)]

    mod_validations[[i]]$predictions <- predictions
    mod_validations[[i]]$predictions_close <- predictions_close
    mod_validations[[i]]$mod_sizes <- mod_sizes
    mod_validations[[i]]$rsq <- summary(fragments_trace$local_southern_mod[[i]]$mod)$r.squared
    mod_validations[[i]]$ladder_hits <- ladder_hits
    mod_validations[[i]]$mod <- fragments_trace$local_southern_mod[[i]]$mod
  }

  # predicted to be a good model if it has some ladder hits and good rsq
  # use the good models to predict the rest of the peaks
  valid_models_tf <- sapply(mod_validations, function(x) ifelse(length(x$ladder_hits) > 0 & x$rsq > rsq_threshold, TRUE, FALSE))
  valid_models <- mod_validations[which(valid_models_tf)]
  predicted_sizes_list <- lapply(valid_models, function(x) {
    sizes <- vector("numeric", length(x$predictions))
    for (i in seq_along(x$predictions)) {
      if (x$predictions[i] %in% x$predictions_close) {
        sizes[i] <- x$predictions[i]
      } else {
        sizes[i] <- NA_real_
      }
    }

    return(sizes)
  })


  predicted_sizes_matrix <- do.call(cbind, predicted_sizes_list)
  predicted_sizes_avg <- numeric(nrow(predicted_sizes_matrix))
  for (i in seq_along(predicted_sizes_avg)) {
    predicted_sizes_avg[i] <- median(predicted_sizes_matrix[i, ],
      na.rm = TRUE
    )
  }

  confirmed_sizes <- unique(unlist(lapply(valid_models, function(x) x$mod_sizes)))
  unconfirmed_sizes <- ladder_sizes[which(!ladder_sizes %in% confirmed_sizes)]

  assigned_size <- sapply(predicted_sizes_avg, function(x) {
    diff <- ladder_sizes - x
    ladder_hit <- ladder_sizes[which(diff > -size_tolerance * 2 & diff < size_tolerance * 2)]
    if (length(ladder_hit) == 0) {
      return(NA_real_)
    } else {
      return(ladder_hit)
    }
  })

  assigned_df <- data.frame(
    scan = ladder_peaks[which(assigned_size %in% unconfirmed_sizes)],
    size = assigned_size[which(assigned_size %in% unconfirmed_sizes)]
  )

  # bind rows back with the sizes that were inferred to be correct

  ladder_df <- rbind(
    assigned_df,
    fragments_trace$ladder_df[which(fragments_trace$ladder_df$size %in% confirmed_sizes), ]
  )

  # now just rerun the bp sizing
  fixed_fragment_trace <- ladder_fix_helper(
    fragments_trace,
    ladder_df
  )

  return(fixed_fragment_trace)
}



# ladder ------------------------------------------------------------------


#' Ladder and bp sizing
#'
#' Find the ladder peaks in and use that to call bp size
#'
#' @param fragments_trace list from 'read_fsa' function
#' @param ladder_channel string: which channel in the fsa file contains the
#'        ladder signal
#' @param signal_channel string: which channel in the fsa file contains the data
#'        signal
#' @param ladder_sizes numeric vector: bp sizes of ladder used in fragment analysis.
#'        defaults to GeneScan™ 500 LIZ™
#' @param ladder_start_scan numeric: indicate the scan number to start looking for
#'        ladder peaks. Usually this can be automatically found (when set to NULL) since
#'        there's a big spike right at the start. However, if your ladder peaks
#'        are taller than the big spike, you will need to set this starting scan
#'        number manually.
#' @param minimum_peak_signal numeric: minimum height of peak from smoothed signal.
#' @param zero_floor logical: if set to TRUE, all negative values will be set to zero.
#'        This can help deal with cases where there are peaks in the negative direction
#'        that interfere with peak detection.
#' @param scan_subset numeric vector (length 2): filter the ladder and data signal
#'        between the selected scans (eg scan_subset = c(3000, 5000)).
#'        to pracma::savgol().
#' @param max_combinations numeric: what is the maximum number of ladder
#'        combinations that should be tested
#' @param ladder_selection_window numeric: in the ladder assigning algorithm,
#'        the we iterate through the scans in blocks and test their linear fit ( We can assume that the ladder is linear over a short distance)
#'        This value defines how large that block of peaks should be.
#' @param smoothing_window numeric: ladder signal smoothing window size passed
#' @param warning_rsq_threshold The value for which this function will warn you when parts of the ladder have R-squared values below the specified threshold.
#' @param show_progress_bar show progress bar
#'
#' @return This function modifies list of fragments_trace objects in place with the ladder assigned and base pair calculated.
#' @export
#'
#' @details
#' This function takes a list of fragments_trace files (the output from read_fsa) and identifies
#' the ladders in the ladder channel which is used to call the bp size. The output
#' is a list of fragments_traces. bp sizes are assigned using the local Southern
#' method. Basically, for each data point, linear models are made for the lower
#' and upper 3 size standard and the predicted sizes are averaged.
#'
#' Use [plot_data_channels()] to plot the raw data on the fsa file to identify which channel the ladder and data are in.
#'
#' The ladder peaks are assigned from largest to smallest. I would recommend excluding
#' size standard peaks less than 50 bp (eg size standard 35 bp).
#'
#' Each ladder should be manually inspected to make sure that is has been correctly
#' assigned.
#'
#' @seealso [plot_data_channels()] to plot the raw data in all channels. [plot_ladders()] to plot the assigned ladder
#' peaks onto the raw ladder signal. [fix_ladders_interactive()] to fix ladders with
#' incorrectly assigned peaks.
#'
#'
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' # Manually inspect the ladders
#' plot_ladders(fsa_list[1])
#'
find_ladders <- function(
    fragments_trace,
    ladder_channel = "DATA.105",
    signal_channel = "DATA.1",
    ladder_sizes = c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    ladder_start_scan = NULL,
    minimum_peak_signal = NULL,
    zero_floor = FALSE,
    scan_subset = NULL,
    ladder_selection_window = 5,
    max_combinations = 2500000,
    smoothing_window = 21,
    warning_rsq_threshold = 0.998,
    show_progress_bar = TRUE) {

  fit_ladder <- function(
      ladder,
      scans,
      sample_id) {
    if (is.null(ladder_start_scan)) {
      ladder_start_scan <- which.max(ladder) + 50
    }

    if (zero_floor) {
      ladder <- pmax(ladder, 0)
    }

    ladder_df <- data.frame(signal = ladder, scan = scans)
    ladder_df <- ladder_df[which(ladder_df$scan >= ladder_start_scan), ]
    ladder_df$detrended_signal <- detrend_signal(ladder_df$signal)
    ladder_df$smoothed_signal <- pracma::savgol(
      ladder_df$detrended_signal,
      smoothing_window
    )

    ladder_peaks <- find_ladder_peaks(
      ladder_df = ladder_df,
      n_reference_sizes = length(ladder_sizes),
      minimum_peak_signal = minimum_peak_signal,
      sample_id = sample_id
    )

    peaks_fit_df <- ladder_iteration(
      reference_sizes = rev(ladder_sizes), # start away from the spike going backwards
      observed_sizes = rev(ladder_peaks),
      choose = ladder_selection_window,
      max_combinations = max_combinations
    )

    peaks_not_fit <- ladder_peaks[which(!ladder_peaks %in% peaks_fit_df$scan)]

    peaks_not_fit_df <- data.frame(
      scan = peaks_not_fit,
      size = rep(NA_real_, length(peaks_not_fit))
    )


    combined_ladder_peaks <- rbind(peaks_fit_df, peaks_not_fit_df)
    combined_ladder_peaks <- combined_ladder_peaks[order(combined_ladder_peaks$scan), ]

    return(combined_ladder_peaks)
  }

  for (i in seq_along(fragments_trace)) {
    if (show_progress_bar) {
      pb <- utils::txtProgressBar(min = 0, max = length(fragments_trace), style = 3)
    }

    # populate the ladder and data channels with the supplied channel name

    fragments_trace[[i]]$raw_ladder <- fragments_trace[[i]]$fsa$Data[[ladder_channel]]
    fragments_trace[[i]]$raw_data <- fragments_trace[[i]]$fsa$Data[[signal_channel]]
    fragments_trace[[i]]$scan <- 0:(length(fragments_trace[[i]]$fsa$Data[[signal_channel]]) - 1)
    fragments_trace[[i]]$off_scale_scans <- fragments_trace[[i]]$fsa$Data$OfSc.1

    # make sure that the scan window is at least same length as length of size standards
    if (ladder_selection_window > length(ladder_sizes)) {
      ladder_selection_window <- length(ladder_sizes)
    }

    # allow user to subset to particular scans
    if (!is.null(scan_subset)) {
      fragments_trace[[i]]$raw_ladder <- fragments_trace[[i]]$raw_ladder[scan_subset[1]:scan_subset[2]]
      fragments_trace[[i]]$raw_data <- fragments_trace[[i]]$raw_data[scan_subset[1]:scan_subset[2]]
      fragments_trace[[i]]$scan <- fragments_trace[[i]]$scan[scan_subset[1]:scan_subset[2]]

      # set spike location since it's automatically set usually, and user may select scans to start after
      ladder_start_scan <- scan_subset[1]
    }

    # ladder
    ladder_df <- fit_ladder(
      ladder = fragments_trace[[i]]$raw_ladder,
      scans = fragments_trace[[i]]$scan,
      sample_id = fragments_trace[[i]]$unique_id
    )

    fragments_trace[[i]]$ladder_df <- ladder_df

    # predict bp size
    ladder_df <- ladder_df[which(!is.na(ladder_df$size)), ]
    ladder_df <- ladder_df[which(!is.na(ladder_df$scan)), ]
    fragments_trace[[i]]$local_southern_mod <- local_southern(ladder_df$scan, ladder_df$size)

    predicted_size <- local_southern_predict(local_southern_output = fragments_trace[[i]]$local_southern_mod, scans = fragments_trace[[i]]$scan)

    fragments_trace[[i]]$trace_bp_df <- data.frame(
      unique_id = rep(fragments_trace[[i]]$unique_id, length(fragments_trace[[i]]$scan)),
      scan = fragments_trace[[i]]$scan,
      size = predicted_size,
      signal = fragments_trace[[i]]$raw_data,
      ladder_signal = fragments_trace[[i]]$raw_ladder,
      off_scale = fragments_trace[[i]]$scan %in% fragments_trace[[i]]$off_scale_scans
    )

    # make a warning if one of the ladder modes is bad
    ladder_rsq_warning_helper(fragments_trace[[i]],
      rsq_threshold = warning_rsq_threshold
    )

    if (show_progress_bar) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  invisible()
}



#' Fix ladders manually
#'
#' Manually assign the ladder peaks for samples in a fragments_trace_list
#'
#' @param fragments_trace_list list of fragments_trace objects
#' @param ladder_df_list a list of dataframes, with the names being the unique id
#' and the value being a dataframe. The dataframe has two columns, size (indicating
#' the bp of the standard) and scan (the scan value of the ladder peak). It's
#' critical that the element name in the list is the unique id of the sample.
#'
#' @return This function modifies list of fragments_trace objects in place with the selected ladders fixed.
#' @export
#'
#' @details
#' This function returns a fragments_trace list the same length as was supplied.
#' It goes through each sample and either just returns the same fragments_trace
#' if the unique id doesn't match the samples that need the ladder fixed, or if
#' it is one to fix, it will use the supplied dataframe in the ladder_df_list
#' as the ladder. It then reruns the bp sizing methods on those samples.
#'
#' This is best used with [fix_ladders_interactive()] that can generate a `ladder_df_list`.
#'
#'
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
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
#'   "20230413_A01.fsa" = data.frame(
#'     size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#'     scan = c(1971, 2208, 2269, 2329, 2581, 2888, 3228, 3479, 3543, 3872, 4170, 4412, 4460)
#'   )
#' )
#'
#' fix_ladders_manual(
#'   fsa_list,
#'   example_list
#' )
#'
fix_ladders_manual <- function(fragments_trace_list,
                               ladder_df_list) {
  samples_to_fix <- names(ladder_df_list)
  for (i in seq_along(fragments_trace_list)) {
    if (fragments_trace_list[[i]]$unique_id %in% samples_to_fix) {
      message(paste("Fixing ladder for", fragments_trace_list[[i]]$unique_id))

      tmp_ladder_df <- ladder_df_list[[which(names(ladder_df_list) == fragments_trace_list[[i]]$unique_id)]]

      # do some quality control of the df user supplied
      if (!any(colnames(tmp_ladder_df) == "scan") | !any(colnames(tmp_ladder_df) == "size")) {
        stop(
          call. = FALSE,
          "Dataframe must contain columns 'size' and 'scan'"
        )
      }

      fragments_trace_list[[i]] <-  ladder_fix_helper(
        fragments_trace_list[[i]],
        replacement_ladder_df = tmp_ladder_df
        )

    }
  }

  invisible()
}






