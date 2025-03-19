# instability index ---------------------------------------------------------
instability_index <- function(repeats,
                              signals,
                              index_peak_signal,
                              index_peak_repeat,
                              peak_threshold,
                              abs_sum = FALSE) {
  # apply signal threshold
  peak_over_threshold <- which(signals / index_peak_signal > peak_threshold)
  repeats <- repeats[peak_over_threshold]
  signals <- signals[peak_over_threshold]

  # normalized peak signal
  signals_normalized <- signals / sum(signals)

  # distance to index peak
  repeat_delta <- repeats - index_peak_repeat
  if (abs_sum == FALSE) {
    sum(signals_normalized * repeat_delta)
  } else if (abs_sum == TRUE) {
    sum(abs(signals_normalized * repeat_delta))
  }
}

# function for finding quantiles -----------------------------------------------

find_percentiles <- function(repeats,
                             signals,
                             index_peak_repeat,
                             type, # "percentile" or "repeat"
                             range,
                             col_prefix) {
  # if there are double peaks select the tallest of the peaks, otherwise approx interpolation doesn't work
  # also if there are no main peak called, filter out (for example samples used to call repeats but irrelevant for metrics)
  df_names <- paste(col_prefix, range, sep = "_")

  # Deal with case when there are no expansion peaks by returning 0s
  if (sum(repeats > index_peak_repeat) <= 1) {
    percentile_df <- as.data.frame(setNames(as.list(rep(0, length(range))), df_names))
  } else {
    unique_repeat_df <- aggregate(signals ~ repeats, FUN = max)
    cumsum_pct <- cumsum(unique_repeat_df$signals) / sum(unique_repeat_df$signals)
    repeat_delta <- unique_repeat_df$repeats - index_peak_repeat

    values <- vector("numeric", length(range))

    if (type == "percentile") {
      for (i in seq_along(range)) {
        values[[i]] <- approx(cumsum_pct,
          repeat_delta,
          xout = range[[i]],
          yleft = min(repeat_delta)
        )$y
      }
    } else if (type == "repeat") {
      for (i in seq_along(range)) {
        values[[i]] <- approx(repeat_delta,
          cumsum_pct,
          xout = range[[i]],
          yleft = min(cumsum_pct)
        )$y
      }
    }

    percentile_df <- as.data.frame(setNames(as.list(values), df_names))
  }

  return(percentile_df)
}


# skewness ------------------------------------------------------------------

fishers_skewness <- function(x, y) {
  mean_val <- sum(x * y)
  sd_val <- sqrt(sum(y * (x - mean_val)^2))

  skewness <- sum(y * (x - mean_val)^3) / sd_val^3

  return(skewness)
}


# kurtosis -----------------------------------------------------------------

fishers_kurtosis <- function(x, y) {
  mean_val <- sum(x * y)
  sd_val <- sqrt(sum(y * (x - mean_val)^2))

  kurtosis <- (sum(y * (x - mean_val)^4) / sd_val^4) - 3
  return(kurtosis)
}

# subsetting repeat table ---------------------------------------------------

repeat_table_subset <- function(repeat_table_df,
                                allele_signal,
                                index_repeat,
                                peak_threshold,
                                window_around_index_peak) {
  # Filter to include only the peaks above the certain threshold
  # signal threshold is set on the modal peak rather than the index peak
  repeat_table_df$peak_percent <- repeat_table_df$signal / allele_signal
  signal_filtered_df <- repeat_table_df[which(repeat_table_df$peak_percent > peak_threshold), ]

  # Ensure window_around_index_peak is exactly length 2
  if (length(window_around_index_peak) != 2) {
    stop("window_around_index_peak must be a vector of length 2")
  }

  # Filter to include only peaks of a certain size
  lower_lim <- ifelse(is.na(window_around_index_peak[1]),
    min(signal_filtered_df$repeats),
    index_repeat - abs(window_around_index_peak[1])
  )
  upper_lim <- ifelse(is.na(window_around_index_peak[1]),
    max(signal_filtered_df$repeats),
    index_repeat + abs(window_around_index_peak[2])
  )
  size_filtered_df <- signal_filtered_df[which(signal_filtered_df$repeats >= lower_lim & signal_filtered_df$repeats <= upper_lim), ]

  return(size_filtered_df)
}

# Calculate metrics -------------------------------------------------------

#' Calculate Repeat Instability Metrics
#'
#' This function computes instability metrics from a list of fragments data objects.
#'
#' @param fragments_list A list of "fragments" objects representing fragment data.
#' @param peak_threshold The threshold for peak signals to be considered in the calculations, relative to the modal peak signal of the expanded allele.
#' @param window_around_index_peak A numeric vector (length = 2) defining the range around the index peak. First number specifies repeats before the index peak, second after. For example, \code{c(-5, 40)} around an index peak of 100 would analyze repeats 95 to 140. The sign of the numbers does not matter (The absolute value is found).
#' @param percentile_range A numeric vector of percentiles to compute (e.g., c(0.5, 0.75, 0.9, 0.95)).
#' @param repeat_range A numeric vector specifying ranges of repeats for the inverse quantile computation.
#'
#' @return A data.frame with calculated instability metrics for each sample.
#' @details
#' Each of the columns in the supplied dataframe are explained below:
#'
#' ## General Information
#' - `unique_id`: A unique identifier for the sample (usually the fsa file name).
#'
#' ## Quality Control
#' - `QC_comments`: Quality control comments.
#' - `QC_modal_peak_signal`: Quality control status based on the modal peak signal (Low < 500, very low < 100).
#' - `QC_peak_number`: Quality control status based on the number of peaks (Low < 20, very low < 10).
#' - `QC_off_scale`: Quality control comments for off-scale peaks. Potential peaks that are off-scale are given. However, a caveat is that this could be from any of the channels (ie it could be from the ladder channel but is the same scan as the given repeat).
#'
#' ## General sample metrics
#' - `modal_peak_repeat`: The repeat size of the modal peak.
#' - `modal_peak_signal`: The signal of the modal peak.
#' - `index_peak_repeat`: The repeat size of the index peak (the repeat value closest to the modal peak of the index sample).
#' - `index_peak_signal`: The signal of the index peak.
#' - `index_weighted_mean_repeat`: The weighted mean repeat size (weighted on the signal of the peaks) of the index sample.
#' - `n_peaks_total`: The total number of peaks in the repeat table.
#' - `n_peaks_analysis_subset`: The number of peaks in the analysis subset.
#' - `n_peaks_analysis_subset_expansions`: The number of expansion peaks in the analysis subset.
#' - `min_repeat`: The minimum repeat size in the analysis subset.
#' - `max_repeat`: The maximum repeat size in the analysis subset.
#' - `mean_repeat`: The mean repeat size in the analysis subset.
#' - `weighted_mean_repeat`: The weighted mean repeat size (weight on peak signal) in the analysis subset.
#' - `median_repeat`: The median repeat size in the analysis subset.
#' - `max_signal`: The maximum peak signal in the analysis subset.
#' - `max_delta_neg`: The maximum negative delta to the index peak.
#' - `max_delta_pos`: The maximum positive delta to the index peak.
#' - `skewness`: The skewness of the repeat size distribution.
#' - `kurtosis`: The kurtosis of the repeat size distribution.
#'
#' ## Repeat instability metrics
#' - `modal_repeat_change`: The difference between the modal repeat and the index repeat.
#' - `average_repeat_change`: The weighted mean of the sample (weighted by peak signal) subtracted by the weighted mean repeat of the index sample(s).
#' - `instability_index_change`: The instability index of the sample subtracted by the instability index of the index sample(s). This will be very similar to the average_repeat_change, with the key difference of instability_index_change being that it is an internally calculated metric for each sample, and therefore the random slight fluctuations of bp size (or systematic if across plates for example) will be removed. However, it requires the index peak to be correctly set for each sample, and if set incorrectly, can produce large arbitrary differences.  
#' - `instability_index`: The instability index based on peak signal and distance to the index peak. (See Lee et al., 2010, \doi{10.1186/1752-0509-4-29}).
#' - `instability_index_abs`: The absolute instability index. The absolute value is taken for the "Change from the main allele".
#' - `expansion_index`: The instability index for expansion peaks only.
#' - `contraction_index`: The instability index for contraction peaks only.
#' - `expansion_ratio`: The ratio of expansion peaks' signals to the main peak signal. Also known as "peak proportional sum" (See Genetic Modifiers of Huntington’s Disease (GeM-HD) Consortium, 2019, \doi{10.1016/j.cell.2019.06.036}).
#' - `contraction_ratio`: The ratio of contraction peaks' signals to the main peak signal.
#' - `expansion_percentile_*`: The repeat size at specified percentiles of the cumulative distribution of expansion peaks.
#' - `expansion_percentile_for_repeat_*`: The percentile rank of specified repeat sizes in the distribution of expansion peaks.
#'
#' @export
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#' # import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
#' test_fragments <- trace_main(fsa_list, grouped = TRUE, metadata_data.frame = metadata)
#'
#' test_metrics_grouped <- calculate_instability_metrics(
#'   fragments_list = test_fragments,
#'   peak_threshold = 0.05,
#'   window_around_index_peak = c(-40, 40)
#' )
calculate_instability_metrics <- function(
    fragments_list,
    peak_threshold = 0.05,
    window_around_index_peak = c(NA, NA),
    percentile_range = c(0.5, 0.75, 0.9, 0.95),
    repeat_range = c(2, 5, 10, 20)) {
  # calculate metrics
  metrics_list <- lapply(fragments_list, function(fragments_repeats) {

    # return early under different situations and record a reason why
    if (nrow(fragments_repeats$repeat_table_df) == 0) {
      fragments_repeats$.__enclos_env__$private$metrics_qc_message <- "Skip reason: sample has no data"
      return(NULL)
    } else if (is.na(fragments_repeats$get_allele_peak()$allele_repeat)) {
      fragments_repeats$.__enclos_env__$private$metrics_qc_message <- "Skip reason: no allele found in sample"
      return(NULL)
    } else if (fragments_repeats$.__enclos_env__$private$assigned_index_peak_grouped == TRUE &&  is.na(fragments_repeats$get_index_peak()$index_repeat)){
      # because of the warning above we know that there is data in this sample, but the issue came from the index grouping
      fragments_repeats$.__enclos_env__$private$metrics_qc_message <- "Skip reason: Invalid index peak in sample grouping. Issue likely with `metrics_baseline_control` sample(s) that pairs with this sample."
      return(NULL)
    }

    # no issues so set this as blank in case calculate_instability_metrics was run with an issue previously
    fragments_repeats$.__enclos_env__$private$metrics_qc_message <- NA_character_



    # filter dataset to user supplied thresholds
    size_filtered_df <- repeat_table_subset(
      repeat_table_df = fragments_repeats$repeat_table_df,
      allele_signal = fragments_repeats$get_allele_peak()$allele_signal,
      index_repeat = fragments_repeats$get_index_peak()$index_repeat,
      peak_threshold = peak_threshold,
      window_around_index_peak = window_around_index_peak
    )

    # filter and calculate index samples if they exist
    if(!is.null(fragments_repeats$.__enclos_env__$private$index_samples) && length(fragments_repeats$.__enclos_env__$private$index_samples) > 0){

      # filter for index samples with data
      not_na_allele <- sapply(fragments_repeats$.__enclos_env__$private$index_samples, function(x) !is.na(x$allele_repeat))

      ## add height and signal sum filter here!
      ## also check to see if that completely removes them and makes those NA


      index_sample_list_filtered <- fragments_repeats$.__enclos_env__$private$index_samples[not_na_allele]

      if(length(index_sample_list_filtered) > 0){
        index_sample_list_filtered <- lapply(index_sample_list_filtered, function(x){
          list(
            x$allele_repeat,
            x$allele_signal,
            repeat_table_subset(
              repeat_table_df = x$repeat_table_df,
              allele_signal = x$allele_signal,
              index_repeat = x$allele_repeat,
              peak_threshold = peak_threshold,
              window_around_index_peak = window_around_index_peak
            )
          )
        })
  
  
        control_weighted_mean_repeat <- sapply(index_sample_list_filtered, function(x){
          weighted.mean(x[[3]]$repeats, x[[3]]$signal)
        })
        index_weighted_mean_repeat <- median(control_weighted_mean_repeat, na.rm = TRUE)
  
        control_instability_index <- sapply(index_sample_list_filtered, function(x){
          instability_index(
            # can use the modal as the index peak since these are the index samples
            repeats = x[[3]]$repeats,
            signals = x[[3]]$signal,
            index_peak_signal = x[[2]],
            index_peak_repeat = x[[1]],
            peak_threshold = peak_threshold,
            abs_sum = FALSE
          )
        })
        index_instability_index <- median(control_instability_index, na.rm = TRUE)
      } else{
        index_weighted_mean_repeat <- NA
        index_instability_index <- NA
      }
    } else{
      index_weighted_mean_repeat <- NA
      index_instability_index <- NA
    }

    # first subset to make some dataframe that are just for contractions or expansions
    size_filtered_df$repeat_delta_index_peak <- size_filtered_df$repeats - fragments_repeats$get_index_peak()$index_repeat
    expansion_filtered <- size_filtered_df[which(size_filtered_df$repeat_delta_index_peak >= 0), ]
    contraction_filtered <- size_filtered_df[which(size_filtered_df$repeat_delta_index_peak <= 0), ]

    # QCs
    QC_modal_peak_signal <- if (fragments_repeats$get_allele_peak()$allele_signal > 500) {
      NA_character_
    } else if (fragments_repeats$get_allele_peak()$allele_signal > 100) {
      "Low"
    } else {
      "Extremely low"
    }

    QC_peak_number <- if (nrow(fragments_repeats$repeat_table_df) > 20) {
      NA_character_
    } else if (nrow(fragments_repeats$repeat_table_df) > 10) {
      "Low"
    } else {
      "Extremely low"
    }

    QC_off_scale <- if (any(fragments_repeats$repeat_table_df$off_scale)) {
      paste(
        "The following repeats were determined off scale (check ladder too, could be scans in any channel):",
        paste(round(fragments_repeats$repeat_table_df[which(fragments_repeats$repeat_table_df$off_scale), "repeats"]), collapse = ", ")
      )
    } else {
      NA_character_
    }

    # make a wide dataframe
    metrics <- data.frame(
      unique_id = fragments_repeats$unique_id,
      QC_comments = NA_character_,
      QC_modal_peak_signal = QC_modal_peak_signal,
      QC_peak_number = QC_peak_number,
      QC_off_scale = QC_off_scale,
      modal_peak_repeat = fragments_repeats$get_allele_peak()$allele_repeat,
      modal_peak_signal = fragments_repeats$get_allele_peak()$allele_signal,
      index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
      index_peak_signal = fragments_repeats$get_index_peak()$index_signal,
      index_weighted_mean_repeat = index_weighted_mean_repeat,
      n_peaks_total = nrow(fragments_repeats$repeat_table_df),
      n_peaks_analysis_subset = nrow(size_filtered_df),
      n_peaks_analysis_subset_expansions = nrow(expansion_filtered),
      min_repeat = min(size_filtered_df$repeats),
      max_repeat = max(size_filtered_df$repeats),
      mean_repeat = mean(size_filtered_df$repeats),
      weighted_mean_repeat = weighted.mean(size_filtered_df$repeats, size_filtered_df$signal),
      median_repeat = median(size_filtered_df$repeats),
      max_signal = max(size_filtered_df$signal),
      max_delta_neg = min(size_filtered_df$repeat_delta_index_peak),
      max_delta_pos = max(size_filtered_df$repeat_delta_index_peak),
      skewness = fishers_skewness(size_filtered_df$repeats, size_filtered_df$signal),
      kurtosis = fishers_kurtosis(size_filtered_df$repeats, size_filtered_df$signal),
      modal_repeat_change = fragments_repeats$get_allele_peak()$allele_repeat - fragments_repeats$get_index_peak()$index_repeat,
      average_repeat_change = weighted.mean(size_filtered_df$repeats, size_filtered_df$signal) - index_weighted_mean_repeat,
      instability_index_change = instability_index(
        repeats = size_filtered_df$repeats,
        signals = size_filtered_df$signal,
        index_peak_signal = fragments_repeats$get_allele_peak()$allele_signal,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ) - index_instability_index,
      instability_index = instability_index(
        repeats = size_filtered_df$repeats,
        signals = size_filtered_df$signal,
        index_peak_signal = fragments_repeats$get_allele_peak()$allele_signal,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ),
      instability_index_abs = instability_index(
        repeats = size_filtered_df$repeats,
        signals = size_filtered_df$signal,
        index_peak_signal = fragments_repeats$get_allele_peak()$allele_signal,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = TRUE
      ),
      expansion_index = instability_index(
        repeats = expansion_filtered$repeats,
        signals = expansion_filtered$signal,
        index_peak_signal = fragments_repeats$get_allele_peak()$allele_signal,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ),
      contraction_index = instability_index(
        repeats = contraction_filtered$repeats,
        signals = contraction_filtered$signal,
        index_peak_signal = fragments_repeats$get_allele_peak()$allele_signal,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ),
      expansion_ratio = sum(expansion_filtered$peak_percent), 
      contraction_ratio = sum(contraction_filtered$peak_percent)
    )

    expansion_percentile <- find_percentiles(
      expansion_filtered$repeats,
      expansion_filtered$signal,
      fragments_repeats$get_index_peak()$index_repeat,
      type = "percentile",
      range = percentile_range,
      col_prefix = "expansion_percentile"
    )

    expansion_repeat <- find_percentiles(
      expansion_filtered$repeats,
      expansion_filtered$signal,
      fragments_repeats$get_index_peak()$index_repeat,
      type = "repeat",
      range = repeat_range,
      col_prefix = "expansion_percentile_for_repeat"
    )

    metrics <- cbind(metrics, expansion_percentile)
    metrics <- cbind(metrics, expansion_repeat)
        
    return(metrics)
  })

  metrics <- do.call(rbind, metrics_list)

  # add back in any samples that were removed earlier or failed to calculate metrics (they are returned as NULL and therefore not in the dataframe)
  missing_samples <- names(fragments_list)[!names(fragments_list) %in% metrics$unique_id]
  if (length(missing_samples) > 0) {
    metrics[nrow(metrics) + seq_along(missing_samples), "unique_id"] <- missing_samples
    rownames(metrics) <- metrics$unique_id

    # add in the reason for skip
    metrics$QC_comments <- sapply(fragments_list, function(x) x$.__enclos_env__$private$metrics_qc_message)[metrics$unique_id]

  }

  return(metrics)
}

