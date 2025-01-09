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

  #allow user to set min signal

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



ladder_fit_cor <- function(fragments_trace){
  ladder_df <- fragments_trace$ladder_df[order(fragments_trace$ladder_df$size),]
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
    fragments_trace,
    rsq_threshold) {
  
  cor_list <- ladder_fit_cor(fragments_trace)
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
        "sample", fragments_trace$unique_id, "has badly fitting ladder for bp sizes:",
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
#' @param minimum_peak_signal numeric: minimum signal of peak from smoothed signal.
#' @param scan_subset numeric vector (length 2): filter the ladder and data signal
#'        between the selected scans (eg scan_subset = c(3000, 5000)).
#'        to pracma::savgol().
#' @param max_combinations numeric: what is the maximum number of ladder
#'        combinations that should be tested
#' @param ladder_selection_window numeric: in the ladder assigning algorithm,
#'        the we iterate through the scans in blocks and test their linear fit ( We can assume that the ladder is linear over a short distance)
#'        This value defines how large that block of peaks should be.
#' @param warning_rsq_threshold The value for which this function will warn you when parts of the ladder have R-squared values below the specified threshold.
#' @param show_progress_bar show progress bar
#'
#' @return This function modifies list of fragments_trace objects in place with the ladder assigned and base pair calculated.
#' @export
#'
#' @details
#' This function takes a list of fragments_trace files (the output from read_fsa) and identifies
#' the ladders in the ladder channel which is used to call the bp size. The output
#' is a list of fragments_traces. 
#' 
#' In this package, base pair (bp) sizes are assigned using a generalized additive model (GAM) with cubic regression splines. The model is fit to known ladder fragment sizes and their corresponding scan positions, capturing the relationship between scan number and bp size. Once trained, the model predicts bp sizes for all scans by interpolating between the known ladder points. This approach provides a flexible and accurate assignment of bp sizes, accommodating the slightly non-linear relationship.
#'
#' Use [plot_data_channels()] to plot the raw data on the fsa file to identify which channel the ladder and data are in.
#'
#' The ladder peaks are assigned from largest to smallest. I would recommend excluding
#' size standard peaks less than 50 bp (eg size standard 35 bp).
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
    scan_subset = NULL,
    ladder_selection_window = 5,
    max_combinations = 2500000,
    warning_rsq_threshold = 0.998,
    show_progress_bar = TRUE) {

  fit_ladder <- function(
      ladder,
      scans,
      sample_id) {
    if (is.null(ladder_start_scan)) {
      ladder_start_scan <- which.max(ladder) + 50
    }

    ladder_df <- data.frame(signal = ladder, scan = scans)
    ladder_df <- ladder_df[which(ladder_df$scan >= ladder_start_scan), ]
    ladder_df$detrended_signal <- detrend_signal(ladder_df$signal)
    ladder_df$smoothed_signal <- pracma::savgol(
      ladder_df$detrended_signal,
      21
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

  if (show_progress_bar) {
    pb <- utils::txtProgressBar(min = 0, max = length(fragments_trace), style = 3)
  }

  for (i in seq_along(fragments_trace)) {
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

    # ladder correlation stats
    # make a warning if one of the ladder modes is bad
    ladder_rsq_warning_helper(fragments_trace[[i]],
      rsq_threshold = warning_rsq_threshold
    )

    predicted_size <- predict_bp_size(
      ladder_df = ladder_df,
      scans = fragments_trace[[i]]$scan
    )

    fragments_trace[[i]]$trace_bp_df <- data.frame(
      unique_id = rep(fragments_trace[[i]]$unique_id, length(fragments_trace[[i]]$scan)),
      scan = fragments_trace[[i]]$scan,
      size = predicted_size,
      signal = fragments_trace[[i]]$raw_data,
      ladder_signal = fragments_trace[[i]]$raw_ladder,
      off_scale = fragments_trace[[i]]$scan %in% fragments_trace[[i]]$off_scale_scans
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
#' @param warning_rsq_threshold The value for which this function will warn you when parts of the ladder have R-squared values below the specified threshold.
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
#'  "20230413_A07.fsa" = data.frame(
#'    size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#'    scan = c(1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)
#'  )
#' )
#'
#' fix_ladders_manual(
#'   fsa_list,
#'   example_list
#' )
#'
fix_ladders_manual <- function(fragments_trace_list,
                               ladder_df_list,
                               warning_rsq_threshold = 0.998) {
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

      fragments_trace_list[[i]]$ladder_df <- tmp_ladder_df
      predicted_size <- predict_bp_size(
        fragments_trace_list[[i]]$ladder_df,
        fragments_trace_list[[i]]$scan
      )
    
      fragments_trace_list[[i]]$trace_bp_df <- data.frame(
        unique_id = rep(fragments_trace_list[[i]]$unique_id, length(fragments_trace_list[[i]]$scan)),
        scan = fragments_trace_list[[i]]$scan,
        size = predicted_size,
        signal = fragments_trace_list[[i]]$raw_data,
        ladder_signal = fragments_trace_list[[i]]$raw_ladder,
        off_scale = fragments_trace_list[[i]]$scan %in% fragments_trace_list[[i]]$off_scale_scans
      )
    
      # make a warning if one of the ladder modes is bad
      ladder_rsq_warning_helper(fragments_trace_list[[i]],
        rsq_threshold = warning_rsq_threshold
      )
    }
  }

  invisible()
}
