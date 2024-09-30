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
  peak_scan_period <- round(median(scans_diffs))


  return(peak_scan_period)
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
  if (is.na(fragments_repeat$get_allele_peak()$allele_size)) {
    df <- data.frame(
      "unique_id" = character(),
      "scan" = numeric(),
      "size" = numeric(),
      "signal" = numeric()
    )

    return(df)
  }

  fragment_window_positions <- which(fragments_repeat$trace_bp_df$size > fragments_repeat$get_allele_peak()$allele_size - fragment_window & fragments_repeat$trace_bp_df$size < fragments_repeat$get_allele_peak()$allele_size + fragment_window)
  window_df <- fragments_repeat$trace_bp_df[fragment_window_positions, ]
  main_peak_scan <- window_df[which(window_df$size == fragments_repeat$get_allele_peak()$allele_size), "scan"]

  peak_scan_period <- find_scan_period(window_df, main_peak_scan)

  pos_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_period,
    direction = 1,
    window = scan_peak_window
  )

  neg_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_period,
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
  if (is.na(fragments_repeat$get_allele_peak()$allele_size)) {
    df <- data.frame(
      "unique_id" = character(),
      "scan" = numeric(),
      "size" = numeric(),
      "signal" = numeric()
    )

    return(df)
  }

  fragment_window_positions <- which(fragments_repeat$trace_bp_df$size > fragments_repeat$get_allele_peak()$allele_size - fragment_window & fragments_repeat$trace_bp_df$size < fragments_repeat$get_allele_peak()$allele_size + fragment_window)
  window_df <- fragments_repeat$trace_bp_df[fragment_window_positions, ]
  main_peak_scan <- window_df[which(window_df$size == fragments_repeat$get_allele_peak()$allele_size), "scan"]


  # determine period

  peak_scan_period <- round(size_period / median(diff(window_df$size)))

  pos_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_period,
    direction = 1,
    window = scan_peak_window
  )

  neg_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
    main_peak_scan,
    peak_scan_period,
    direction = -1,
    window = scan_peak_window
  )

  peak_table <- fragments_repeat$trace_bp_df
  peak_table <- peak_table[which(peak_table$scan %in% c(neg_peaks, main_peak_scan, pos_peaks)), ]
  peak_table <- peak_table[which(peak_table$size > fragments_repeat$.__enclos_env__$private$min_bp_size & peak_table$size < fragments_repeat$.__enclos_env__$private$max_bp_size), ]


  return(peak_table)
}



# find_size_batch_correction_factor

find_batch_correction_factor <- function(
  fragments_list,
  trace_window_size = 50, smoothing_window = 301){
  # make df for all samples of plate id, batch sample id
  metadata_list <- lapply(fragments_list, function(x){
    df <- data.frame(
      unique_id = x$unique_id,
      batch_run_id = x$batch_run_id,
      batch_sample_id = x$batch_sample_id,
      batch_sample_modal_repeat = x$batch_sample_modal_repeat,
      allele_height = x$get_allele_peak()$allele_height
    )
    return(df)
  })

  metadata_df <- do.call(rbind, metadata_list)

  # filter for batch correction samples that have a trace (height > 100). perhaps a little arbitrary
  correction_sample_df <- metadata_df[which(!is.na(metadata_df$batch_run_id) & !is.na(metadata_df$batch_sample_id)), , drop = FALSE]
  correction_sample_df <- correction_sample_df[which(!is.na(correction_sample_df$allele_height) & correction_sample_df$allele_height > 100), ]

  if(length(unique(na.omit(correction_sample_df$batch_run_id))) <= 1){
    message("Batch correction was not carried out. There needs to be more than one 'batch_run_id' that have samples of adequate quality (allele_height > 100)")
    # need to now exit out of this function early and make sure batch correction is FALSE in the outer environment of the other function
    assign("batch_correction", FALSE, envir = parent.frame())
    return(NULL)
  }

  # Use match to align the order of correction_sample_fragments with correction_sample_df
  matched_indices <- match(correction_sample_df$unique_id, names(fragments_list))
  correction_sample_fragments <- fragments_list[matched_indices]
  correction_sample_df$smoothed_modal_size <- sapply(correction_sample_fragments, function(fragment){
      df <- fragment$trace_bp_df[which(fragment$trace_bp_df$size < fragment$get_allele_peak()$allele_size + trace_window_size & fragment$trace_bp_df$size > fragment$get_allele_peak()$allele_size - trace_window_size ), ]
      df$smoothed <- pracma::savgol(df$signal, smoothing_window)
      smoothed_modal_size <- df[which.max(df$smoothed), "size"]
      return(smoothed_modal_size)
    })
  # scale it for the model
  correction_sample_df$smoothed_modal_size_scaled <- correction_sample_df$smoothed_modal_size - median(correction_sample_df$smoothed_modal_size)
  if(length(unique(na.omit(correction_sample_df$batch_sample_id ))) <= 1){
    # Calculate batch effect directly from the smoothed modal sizes within the single batch
    #deal with case where there is replicates of the same sample by first splitting by batch
    correction_sample_df_split <- split(correction_sample_df, correction_sample_df$batch_run_id)
    batch_effect <- sapply(correction_sample_df_split, function(x) median(x$smoothed_modal_size_scaled))

    batch_effects_df <- data.frame(
      batch_run_id = names(correction_sample_df_split),
      batch_effect = batch_effect
    )
  } else{
    # used mixed model rather than fixed effects model for more flexible in handling unbalanced designs and non-overlapping batches.
    model <- lme4::lmer(smoothed_modal_size_scaled ~ batch_sample_id + (1|batch_run_id), data = correction_sample_df)
    batch_effects_df <- data.frame(
      batch_run_id = row.names(lme4::ranef(model)$batch_run_id),
      batch_effect = lme4::ranef(model)$batch_run_id[[1]]
    )
  }

  # do some checks to see if any batches have not been corrected
  fragments_batch_runs <- sapply(fragments_list, function(x) x$batch_run_id)
  unique_fragments_batch_runs <- unique(fragments_batch_runs)
  if(any(!unique_fragments_batch_runs %in% batch_effects_df$batch_run_id)){
    non_corrected_batches <- unique_fragments_batch_runs[which(!unique_fragments_batch_runs %in% batch_effects_df$batch_run_id)]
    if(length(non_corrected_batches) == 1 && is.na(non_corrected_batches)){
      warning(
        call. = FALSE,
        "Samples with 'batch_run_id' 'NA' had no batch correction carried out. If not intended, check metadata and quality of the 'batch_sample_id' samples."
      )
    } else{
      warning(
        call. = FALSE,
        paste0("The following 'batch_run_id' were not corrected: ", non_corrected_batches,
        ". If not intended, check metadata and quality of the 'batch_sample_id' samples.")
      )
    }
  }
  message("Correcting batch effects")
  # save correction factor for each class object but only if it was actually in the mod
  for (i in seq_along(fragments_list)) {
    if(fragments_list[[i]]$batch_run_id %in% batch_effects_df$batch_run_id){
      # Made the plate id explicitly match the list name for cases when the plate name is a number. It could cause subsetting issues
      fragments_list[[i]]$.__enclos_env__$private$batch_correction_factor <- batch_effects_df[which(batch_effects_df$batch_run_id == fragments_list[[i]]$batch_run_id), "batch_effect"]
    } 
  }

  invisible()
}

# repeat length correction -------------------------------------------------------

model_repeat_length <- function(
    fragments_list,
    repeat_size) {
  calling_close_neighbouring_repeats <- function(controls_fragments) {
    # use np_repeat to accurately call the repeat length of the neighboring peaks
    # extract a dataframe of the called repeats that can then be used to make a model
    controls_fragments_df_list <- lapply(controls_fragments, function(x) {
      df_length <- nrow(x$repeat_table_df)
      # identify peaks close to modal peak and at least 20% as high
      main_peak_delta <- x$repeat_table_df$size - x$get_allele_peak()$allele_size
      height_prop <- x$repeat_table_df$height / x$get_allele_peak()$allele_height
      peak_cluster <- vector("logical", length = nrow(x$repeat_table_df))
      for (i in seq_along(main_peak_delta)) {
        if (abs(main_peak_delta[[i]]) < 30 & height_prop[[i]] > 0.2) {
          peak_cluster[[i]] <- TRUE
        } else {
          peak_cluster[[i]] <- FALSE
        }
      }
      cluster_df <- x$repeat_table_df[peak_cluster, ]
      cluster_df_length <- nrow(cluster_df)
      # use np_repeat method to accurately call the neighboring repeats
      data.frame(
        unique_id = rep(x$unique_id, cluster_df_length),
        size = cluster_df$size,
        validated_repeats = np_repeat(
          size = cluster_df$size,
          main_peak_size = x$get_allele_peak()$allele_size,
          main_peak_repeat = x$batch_sample_modal_repeat,
          repeat_size = repeat_size
        ),
        height = cluster_df$height,
        batch_run_id = rep(x$batch_run_id, cluster_df_length)
      )
    })

    controls_repeats_df <- do.call(rbind, controls_fragments_df_list)
  }

  ## first pull out a data.frame for all samples with a column that indicates if it's a positive control or not
  extracted <- lapply(fragments_list, function(x) {
    data.frame(
      unique_id = x$unique_id,
      allele_size = x$get_allele_peak()$allele_size,
      batch_run_id = x$batch_run_id,
      batch_sample_id = x$batch_sample_id,
      batch_sample_modal_repeat = x$batch_sample_modal_repeat
    )
  })
  extracted_df <- do.call(rbind, extracted)

  # Check to see if there are controls, if there are none, give error
  if (!any(!is.na(extracted_df$batch_sample_modal_repeat))) {
    stop("No samples with batch_sample_modal_repeat were detected. Ensure that the metadata has been added to the samples with 'add_metadata()'.",
      call. = FALSE
    )
  }
  # pull out the controls
  controls_df <- extracted_df[which(!is.na(extracted_df$batch_sample_modal_repeat)), , drop = FALSE]
  controls_fragments <- fragments_list[which(names(fragments_list) %in% controls_df$unique_id)]
  controls_repeats_df <- calling_close_neighbouring_repeats(controls_fragments)

  # Check to see if there are controls for each plate, if there are no controls for a plate, give error
  all_batch_run_ids <- lapply(fragments_list, function(x) x$batch_run_id)
  control_batch_run_ids <- unique(controls_repeats_df$batch_run_id)
  if (length(unique(control_batch_run_ids)) != length(unique(all_batch_run_ids))) {
    plates_missing_controls <- paste0(all_batch_run_ids[which(!all_batch_run_ids %in% control_batch_run_ids)], collapse = ", ")
    stop(paste("Plate(s)", plates_missing_controls, "have no repeat-length control samples"),
      call. = FALSE
    )
  }

  # identify size stds with shared id for more quality control
  # can compare to each other to make sure that the same peak in the distribution has been selected as the modal   
  standard_sample_ids <- sapply(controls_fragments, function(x) x$batch_sample_id)
  warning_message_list <- list()
  if (any(!is.na(standard_sample_ids))) {

    # find the smoothed modal size
    controls_df$smoothed_modal_size <- sapply(controls_fragments, function(fragment){
      df <- fragment$trace_bp_df[which(fragment$trace_bp_df$size < fragment$get_allele_peak()$allele_size + 50 & fragment$trace_bp_df$size > fragment$get_allele_peak()$allele_size - 50 ), ]
      df$smoothed <- pracma::savgol(df$signal, 301)
      smoothed_modal_size <- df[which.max(df$smoothed), "size"]
      return(smoothed_modal_size)
    })

    controls_df$modal_offset <- controls_df$smoothed_modal_size - controls_df$allele_size
    controls_by_sample_ids <- split(controls_df, controls_df$batch_sample_id)
    controls_by_sample_ids_looks_off <- sapply(controls_by_sample_ids, function(batch_sample_df){
      sample_diff_whole_repeat_apart <- sapply(batch_sample_df$modal_offset , function(x) any(abs(x - batch_sample_df$modal_offset ) > repeat_size*0.8) )
      return(any(sample_diff_whole_repeat_apart))
    })
    if(any(controls_by_sample_ids_looks_off)){
      off_sample_id <- names(controls_by_sample_ids_looks_off)[which(controls_by_sample_ids_looks_off)]
      # store the warning message to generate a single warning message later in the function
      warning_message_list[[1]] <- paste0(
        "It looks like the following 'batch_sample_id' need inspecting because the modal peak may have changed across 'batch_run_id': ",
        paste(off_sample_id, collapse = ", "), 
        ". It's possible that the modal peak has shifted to a different spot in the distribution in at least of one the runs. ",
        "Use plot_batch_correction_samples() to visualize and identify these samples, and update metadata with the correct repeat length of the modal peak for the appropriate sample."
      )
    }
  }

  message(paste0("Repeat correction model: ", length(unique(controls_repeats_df$unique_id)), " samples used to build model"))

  # Can now make a model based on the bp size and the known repeat size
  if (length(unique(controls_repeats_df$batch_run_id)) == 1) {
    # when there's only one plate just set up none lm
    correction_mods <- stats::lm(validated_repeats ~ size, data = controls_repeats_df)
  } else {
    # when there are multiple samples a linear model can be made using the modal peak and the known repeat length of the modal peak
    correction_mods <- lm(validated_repeats ~ size * batch_run_id, data = controls_repeats_df)
  }
  repeat_bp_size <- round(1 / correction_mods$coefficients[2], 2)
  message(paste0("Repeat correction model: ", repeat_bp_size, " bp increase per repeat"))


  # check to see if any samples look off
  controls_repeats_df$predicted_repeat <- stats::predict.lm(correction_mods, controls_repeats_df)
  controls_repeats_df$residuals <- correction_mods$residuals

  if (any(abs(controls_repeats_df$residuals) > 0.3)) {
    samples_all_controls <- unique(controls_repeats_df$unique_id)
    samples_high_diff <- unique(controls_repeats_df[which(abs(controls_repeats_df$residuals) > 0.5), "unique_id"])

    warning_message_2 <- paste0(
      "The following samples may be off (based on model residuals) and need investigation: ",
      paste(samples_high_diff, collapse = ", "),
      ". It's possible that at least one of these samples has the incorrect repeat length indicated in the metadata."
    )
    if(length(warning_message_list) == 1){
      warning(call. = FALSE, paste0(warning_message_list[[1]], "\n\n", warning_message_2))
    } else{
      warning(call. = FALSE, warning_message_2)
    }
  } else if(length(warning_message_list) == 1){
    warning(call. = FALSE, warning_message_list[[1]])
  }

  for (i in seq_along(fragments_list)) {
    fragments_list[[i]]$.__enclos_env__$private$repeat_correction_mod <- correction_mods
  }

  invisible()
}


# call_repeats ------------------------------------------------------------

#' Call Repeats for Fragments
#'
#' This function calls the repeat lengths for a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#' @param assay_size_without_repeat An integer specifying the assay size without repeat for repeat calling. This is the length of the sequence flanking the repeat in the PCR product.
#' @param repeat_size An integer specifying the repeat size for repeat calling. Default is 3.
#' @param force_whole_repeat_units A logical value specifying if the peaks should be forced to be whole repeat units apart. Usually the peaks are slightly under the whole repeat unit if left unchanged.
#' @param correction A character vector of either "batch" to carry out a batch correction from common samples across runs (known repeat length not required), or "repeat" to use samples with validated modal repeat lengths to correct the repeat length. Requires metadata to be added (see [add_metadata()]) with both "batch" and "repeat" requiring \code{"batch_run_id"}, "batch" requiring (\code{"batch_sample_id"}) and "repeat" requiring \code{"batch_sample_modal_repeat"} (but also benefits from having \code{"batch_sample_id"}).
#' @param batch_correction A logical specifying if the size should be adjusted across fragment analysis runs. Requires metadata to be added to specify samples (\code{"batch_sample_id"}) common across runs (\code{"batch_run_id"})(see [add_metadata()]).
#' @param repeat_calling_algorithm A character specifying the repeat calling algorithm. Options: \code{"none"}, \code{"fft"}, or \code{"size_period"} (see details section for more information on these).
#' @param repeat_calling_algorithm_size_window_around_allele A numeric value for how big of a window around the tallest peak should be used to find the peak periodicity. Used for both \code{"fft"} and \code{"size_period"}. For \code{"fft"}, you want to make sure that this window is limited to where there are clear peaks. For \code{"size_period"}, it will not make a big difference.
#' @param repeat_calling_algorithm_peak_assignment_scan_window A numeric value for the scan window when assigning the peak. This is used for both \code{"fft"} and \code{"size_period"}. When the scan period is determined, the algorithm jumps to the predicted scan for the next peak. This value opens a window of the neighboring scans to pick the tallest in.
#' @param repeat_calling_algorithm_size_period A numeric value \code{"size_period"} algorithm to set the peak periodicity by bp size. This is the key variable to change for \code{"size_period"}. In fragment analysis, the peaks are usually slightly below the actual repeat unit size.
#'
#' @return This function modifies list of fragments objects in place with repeats added.
#'
#' @details
#' This function has a lot of different options features for determining the repeat length of your samples. This includes i) an option to force the peaks to be whole repeat units apart, ii) corrections to correct batch effects or accurately call repeat length by comparing to samples of known length, and iii) algorithms or re-calling the peaks to remove any contaminating peaks or shoulder-peaks.
#'
#' --- force_whole_repeat_units ---
#' 
#' The `force_whole_repeat_units` option aims to correct for the systematic underestimation in fragment sizes that occurs in capillary electrophoresis. It is independent to the algorithms described above and can be used in conjunction. It modifies repeat lengths in a way that helps align peaks with the underlying repeat pattern, making the repeat lengths whole units (rather than ~0.9 repeats). The calculated repeat lengths start from the main peak's repeat length and increases in increments of the specified `repeat_size` in either direction. This option basically enables you to get exactly the same result as expansion_index values calculated from data from Genemapper.
#' 
#' --- correction ---
#'
#' There are two main correction approaches that are somewhat related: either 'batch' or 'repeat'. Batch correction is relatively simple and just requires you to link samples across batches to correct batch-batch variation in repeat sizes. However, even though the repeat size that is return will be precise, it will not be accurate and underestimates the real repeat length. By contrast, repeat correction can be used to accurately call repeat lengths (which also corrects the batch effects). However, the repeat correction will only be as good as your sample used to call the repeat length so this is a challenging and advanced feature. You need to use a sample that reliably returns the same peak as the modal peak, or you need to be willing to understand the shape of the distribution and manually validate the repeat length of each batch_sample_id for each run. 
#' 
#' Batch correction uses common sample(s) across fragment analysis runs to correct systematic batch effects that occur with repeat-containing amplicons in capillary electrophoresis. There are slight fluctuations of size across runs for amplicons containing repeats that result in systematic differences around 1-3 base pairs. So, if samples are to be analyzed for different runs, the absolute bp size is not comparable unless this batch effect is corrected. This is only relevant when the absolute size of a amplicons are compared for grouping metrics as described above (otherwise instability metrics are all relative and it doesn’t matter that there’s systematic batch effects across runs) or when plotting traces from different runs. This correction can be achieved by running a couple of samples in every fragment analysis run, or having a single run that takes a couple of samples from every run together, thereby linking them. These samples are then indicated in the metadata with batch_run_id (to group samples by fragment analysis run) and batch_sample_id (to enable linking samples across batches) (see [add_metadata()]).
#' 
#' Samples with known and validated repeat size can be used to accurately call the repeat length (and therefore also correct batch effects). Similar to batch correction, batch_run_id (to group samples by fragment analysis run) and batch_sample_id (to enable linking samples across batches) are used, but importantly batch_sample_modal_repeat is also set (see [add_metadata()]). The batch_sample_modal_repeat is the validated repeat length of the modal repeat of the sample. This validated repeat length is then used to call the repeat length of the modal repeat for each sample (by each batch_run_id). Importantly, this correction requires you to know with confidence the repeat length of the modal peak of the sample. Therefore it's important that the sample used for repeat correction has a clear and prominent modal peak. If the repeat length is very long, it's common for the modal peak of a sample to change so if you use this feature you're going to have to understand the shape of the distribution of your sample and double check that the correct peak has been called as the modal peak after you have used [find_alleles()]. If a different peak is selected as the modal peak than usual, you need to go back to the metadata and adjust the repeat size of the size standard (For example, your size standard sample has been validated to have 120 repeats. You run [find_alelles()] and look at the distribution of peaks and notice that the peak one repeat unit higher is the modal peak this time. Therefore, you're going to need to set the batch_sample_modal_repeat as 121 in the metadata just for that batch_run_id. In the other runs you would keep the batch_sample_modal_repeat as 120.).
#' 
#' --- repeat_calling_algorithm ---
#' 
#' The `fft` or `size_period` algorithms both re-call the peaks based on empirically determined (`fft`) or specified (`size_period`) periodicity of the peaks. The main application of these algorithms is to solve the issue of contaminating peaks in the expected regular pattern of peaks. The `fft` approach applies a fourier transform to the peak signal to determine the underlying periodicity of the signal. `size_period` is similar and simpler, where instead of automatically figuring out the periodicity, we as users specify the periodicity (since we usually know the size distance between repeat units). We can use the periodicity to jump between peaks.
#' @seealso [find_alleles()], [add_metadata()]
#'
#' @export
#' 
#' @importFrom lme4 lmer
#' @importFrom lme4 ranef
#'
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[c(90:94)], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' fragments_list <- find_fragments(
#'   fsa_list,
#'   min_bp_size = 300
#' )
#'
#' find_alleles(fragments_list)
#' 
#' add_metadata(fragments_list)
#'
#' # Simple conversion from bp size to repeat size
#' call_repeats(
#'   fragments_list,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(fragments_list[1], xlim = c(120, 170))
#'
#' # Use force_whole_repeat_units algorithm to make sure called
#' # repeats are the exact number of bp apart
#'
#' call_repeats(
#'   fragments_list,
#'   force_whole_repeat_units = TRUE,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(fragments_list[1], xlim = c(120, 170))
#'
#' 
#' # apply batch correction
#' call_repeats(
#'   fragments_list,
#'   correction = "batch",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#' 
#' plot_traces(fragments_list[1], xlim = c(120, 170))
#' 
#' # apply repeat correction
#' call_repeats(
#'   fragments_list,
#'   correction = "repeat",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#' 
#' plot_traces(fragments_list[1], xlim = c(120, 170))
#'
#' # use different algorithms to call the repeats to ensure only periodic peaks are called
#'
#' # fft to automatically find peak period
#' call_repeats(
#'   fragments_list,
#'   repeat_calling_algorithm = "fft",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(fragments_list[1], xlim = c(120, 170))
#'
#' # size_period to manually supply the peak period
#' call_repeats(
#'   fragments_list,
#'   repeat_calling_algorithm = "size_period",
#'   repeat_calling_algorithm_size_period = 2.75,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(fragments_list[1], xlim = c(120, 170))
#' 
call_repeats <- function(
    fragments_list,
    assay_size_without_repeat = 87,
    repeat_size = 3,
    force_whole_repeat_units = FALSE,
    correction = "none",
    repeat_calling_algorithm = "none",
    repeat_calling_algorithm_size_window_around_allele = repeat_size * 5,
    repeat_calling_algorithm_peak_assignment_scan_window = 3,
    repeat_calling_algorithm_size_period = repeat_size * 0.93) {
 
  ### in this function, we are doing three key things
      #### 1) use repeat_calling_algorithm to find repeats and generate a new repeat table dataframe
      #### 2) apply batch correction or repeat correction
      #### 3) call repeats with or without force whole repeat units
  
  # check to make sure all the required inputs for the function have been given
  if (fragments_list[[1]]$.__enclos_env__$private$find_main_peaks_used == FALSE) {
    stop(paste0(fragments_list[[1]]$unique_id, " requires main alleles to be identified before repeats can be called. Find alleles using 'find_main_peaks()' within the class, or use 'find_alleles()' to find the main peaks across a list of 'fragments_repeats' objects"),
      call. = FALSE
    )
  }
  
  # first use repeat_calling_algorithm to find repeats and generate a new repeat table dataframe
  fragments_list <- lapply(fragments_list, function(fragment){

      # only continue from here if main peaks were successfully found, otherwise, don't return repeat data (ie it can be an empty df)
      if (is.na(fragment$get_allele_peak()$allele_size) | is.na(fragment$get_allele_peak()$allele_height)) {
        fragment$.__enclos_env__$private$repeats_not_called_reason <- "No main peaks"
        # populate with empty dataframe to help the rest of the pipeline
        fragment$repeat_table_df <- data.frame(
          unique_id = character(),
          size = numeric(),
          height = numeric(),
          calculated_repeats = numeric(),
          off_scale = logical()
        )

        # exit lapply early 
        return(fragment)
      } 
      # repeat calling algorithm
      if (repeat_calling_algorithm == "none") {

        repeat_table_df <- data.frame(
          unique_id = fragment$peak_table_df$unique_id,
          size = fragment$peak_table_df$size, 
          height = fragment$peak_table_df$height,
          calculated_repeats = (fragment$peak_table_df$size- assay_size_without_repeat) / repeat_size,
          off_scale = ifelse(any(colnames(fragment$peak_table_df) == "off_scale"),
          fragment$peak_table_df$off_scale,
            rep(FALSE, nrow(fragment$peak_table_df))
          )
        )
      } else if (repeat_calling_algorithm == "fft") {
        # check to see that fragments repeats has trace data since that is required.
        if (is.null(fragment$trace_bp_df)) {
          stop("fft algorithm requires trace data. Use fsa samples rather than peak table is inputs into the pipeline.",
            call. = FALSE
          )
        }

        # this function generates a new peak table
        fft_peak_df <- fft_repeat_caller(fragment,
          fragment_window = repeat_calling_algorithm_size_window_around_allele,
          scan_peak_window = repeat_calling_algorithm_peak_assignment_scan_window
        )

        repeat_table_df <- data.frame(
          unique_id = fft_peak_df$unique_id,
          size = fft_peak_df$size, 
          height = fft_peak_df$signal,
          calculated_repeats = (fft_peak_df$size - assay_size_without_repeat) / repeat_size,
          off_scale = fft_peak_df$off_scale
        )
      } else if (repeat_calling_algorithm == "size_period") {
        # check to see that fragments repeats has trace data since that is required.
        if (is.null(fragment$trace_bp_df)) {
          stop("size_period algorithm requires trace data. Use fsa samples rather than peak table is inputs into the pipeline.",
            call. = FALSE
          )
        }
        size_period_df <- size_period_repeat_caller(fragment,
          size_period = repeat_calling_algorithm_size_period,
          fragment_window = repeat_calling_algorithm_size_window_around_allele,
          scan_peak_window = repeat_calling_algorithm_peak_assignment_scan_window
        )
        repeat_table_df <- data.frame(
          unique_id = size_period_df$unique_id,
          size = size_period_df$size, #use the original size so that the correction in the repeat table is only applied to repeats rather than size
          height = size_period_df$signal,
          calculated_repeats = (size_period_df$size - assay_size_without_repeat) / repeat_size,
          off_scale = size_period_df$off_scale
        )
      } else {
        stop(
          call. = FALSE,
          "Invalid repeat calling algorithm selected"
        )
      }
      fragment$repeat_table_df <- repeat_table_df
      return(fragment)
  })

  # now we can do #2 and find correction factors
  if (correction == "repeat") {
    model_repeat_length(
      fragments_list = fragments_list,
      repeat_size = repeat_size
    )
  } else if(correction == "batch"){
    find_batch_correction_factor(fragments_list)
  } else if(correction != "none"){
    stop(call. = FALSE, "Invalid correction type. Select either 'repeat' or 'batch'")
  }

  # call repeats for each sample
  fragments_list <- lapply(
    fragments_list,
    function(fragment) {
      repeat_table_df <- fragment$repeat_table_df
      # only continue from here if there actually is data
      if(nrow(repeat_table_df) == 0){
        fragment$repeat_table_df$repeats <- numeric()
        return(fragment) # return early
      }

      if(correction == "batch"){
        # re-calculate calculated_repeats repeats but now including batch correction
        repeat_table_df$calculated_repeats <- (repeat_table_df$size - assay_size_without_repeat - fragment$.__enclos_env__$private$batch_correction_factor) / repeat_size
      } else if(correction == "repeat"){
        # Predicted modal repeat size and calculate a repeat correction factor
        repeat_table_df$batch_run_id <- rep(fragment$batch_run_id, nrow(repeat_table_df))
        modal_row_df <- repeat_table_df[which(repeat_table_df$size == fragment$get_allele_peak()$allele_size), ]
        predicted_modal_repeat <- stats::predict.lm(
          fragment$.__enclos_env__$private$repeat_correction_mod, 
          modal_row_df
        )
        fragment$.__enclos_env__$private$repeat_correction_factor <- predicted_modal_repeat - modal_row_df$calculated_repeats

        # apply correction factor to all calculated repeats
        repeat_table_df$calculated_repeats  <- repeat_table_df$calculated_repeats + fragment$.__enclos_env__$private$repeat_correction_factor         
      }

      # Finally call repeats with or without forcing whole repeat units
      if (force_whole_repeat_units) {
        repeat_table_df$repeats <- np_repeat(
          size = repeat_table_df$size,
          main_peak_size = fragment$get_allele_peak()$allele_size,
          main_peak_repeat = repeat_table_df$calculated_repeats[which(repeat_table_df$size == fragment$get_allele_peak()$allele_size)],
          repeat_size = repeat_size
        )
      } else{
        repeat_table_df$repeats <- repeat_table_df$calculated_repeats
      }

      # Finally save main peak repeat length and repeats data
      fragment$repeat_table_df <- repeat_table_df
      allele_subset <- repeat_table_df$repeats[which(repeat_table_df$size == fragment$get_allele_peak()$allele_size)]
      fragment$set_allele_peak(unit = "repeats", value = allele_subset)
      
      # save useful info that is used elsewhere
      fragment$.__enclos_env__$private$repeat_size <- repeat_size
      fragment$.__enclos_env__$private$assay_size_without_repeat <- assay_size_without_repeat

      return(fragment)
    }
  )


  # need to go over samples and apply repeat to all traces if it exists
  if(correction == "repeat"){
    # need to figure out correction factor for samples that repeat lengths were not called because no alleles
    repeat_correction_list <- lapply(fragments_list, function(x){
      data.frame(batch_run_id = x$batch_run_id, repeat_correction_factor = x$.__enclos_env__$private$repeat_correction_factor)
    })
    repeat_correction_df <- do.call(rbind, repeat_correction_list)
    repeat_correction_factor_by_batch <- lapply(
      split(repeat_correction_df, repeat_correction_df$batch_run_id), 
      function(x) median(x$repeat_correction_factor, na.rm = TRUE)
    )
    for (i in seq_along(fragments_list)) {
      if(is.na(fragments_list[[i]]$.__enclos_env__$private$repeat_correction_factor)){
        fragments_list[[i]]$.__enclos_env__$private$repeat_correction_factor <- repeat_correction_factor_by_batch[[fragments_list[[i]]$batch_run_id]]
      }
      fragments_list[[i]]$trace_bp_df$calculated_repeats <- (fragments_list[[i]]$trace_bp_df$size - assay_size_without_repeat) / repeat_size
      fragments_list[[i]]$trace_bp_df$calculated_repeats <- fragments_list[[i]]$trace_bp_df$calculated_repeats + fragments_list[[i]]$.__enclos_env__$private$repeat_correction_factor 
    }
  } else if(correction == "batch"){
    for (i in seq_along(fragments_list)) {
      fragments_list[[i]]$trace_bp_df$calculated_repeats <- (fragments_list[[i]]$trace_bp_df$size - assay_size_without_repeat - fragments_list[[i]]$.__enclos_env__$private$batch_correction_factor) / repeat_size
    }
  } else{
    for (i in seq_along(fragments_list)) {
      fragments_list[[i]]$trace_bp_df$calculated_repeats <- (fragments_list[[i]]$trace_bp_df$size - assay_size_without_repeat) / repeat_size
    }
  }

  # loop over samples to give appropriate warnings about certain events
  repeats_not_called_reason <- sapply(fragments_list, function(x) x$.__enclos_env__$private$repeats_not_called_reason)
  if(any(repeats_not_called_reason %in% "No main peaks")){
    warning(
      paste0(
        "Repeats were not called in the following samples (no allele in sample): ", 
        paste0(names(repeats_not_called_reason)[which(repeats_not_called_reason == "No main peaks")], collapse = ", ")
    ),
      call. = FALSE
    )
  }    

  invisible()
}
