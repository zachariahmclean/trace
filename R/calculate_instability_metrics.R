# instability index ---------------------------------------------------------
instability_index <- function(repeats,
                              heights,
                              index_peak_height,
                              index_peak_repeat,
                              peak_threshold,
                              abs_sum = FALSE) {
  # apply height threshold
  peak_over_threshold <- which(heights / index_peak_height > peak_threshold)
  repeats <- repeats[peak_over_threshold]
  heights <- heights[peak_over_threshold]

  # normalized peak height
  heights_normalized <- heights / sum(heights)

  # distance to index peak
  repeat_delta <- repeats - index_peak_repeat
  if (abs_sum == FALSE) {
    sum(heights_normalized * repeat_delta)
  } else if (abs_sum == TRUE) {
    sum(abs(heights_normalized * repeat_delta))
  }
}

# function for finding quantiles -----------------------------------------------

find_percentiles <- function(repeats,
                             heights,
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
    unique_repeat_df <- aggregate(heights ~ repeats, FUN = max)
    cumsum_pct <- cumsum(unique_repeat_df$heights) / sum(unique_repeat_df$heights)
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
                                allele_height,
                                index_repeat,
                                peak_threshold,
                                window_around_index_peak) {
  # Filter to include only the peaks above the certain threshold
  # height threshold is set on the modal peak rather than the index peak
  repeat_table_df$peak_percent <- repeat_table_df$height / allele_height
  height_filtered_df <- repeat_table_df[which(repeat_table_df$peak_percent > peak_threshold), ]

  # Ensure window_around_index_peak is exactly length 2
  if (length(window_around_index_peak) != 2) {
    stop("window_around_index_peak must be a vector of length 2")
  }

  # Filter to include only peaks of a certain size
  lower_lim <- ifelse(is.na(window_around_index_peak[1]),
    min(height_filtered_df$repeats),
    index_repeat - abs(window_around_index_peak[1])
  )
  upper_lim <- ifelse(is.na(window_around_index_peak[1]),
    max(height_filtered_df$repeats),
    index_repeat + abs(window_around_index_peak[2])
  )
  size_filtered_df <- height_filtered_df[which(height_filtered_df$repeats >= lower_lim & height_filtered_df$repeats <= upper_lim), ]

  return(size_filtered_df)
}

# Calculate metrics -------------------------------------------------------

#' Calculate Repeat Instability Metrics
#'
#' This function computes instability metrics from a list of fragments_repeats data objects.
#'
#' @param fragments_list A list of "fragments_repeats" objects representing fragment data.
#' @param peak_threshold The threshold for peak heights to be considered in the calculations, relative to the modal peak height of the expanded allele.
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
#' - `QC_modal_peak_height`: Quality control status based on the modal peak height (Low < 500, very low < 100).
#' - `QC_peak_number`: Quality control status based on the number of peaks (Low < 20, very low < 10).
#' - `QC_off_scale`: Quality control comments for off-scale peaks. Potential peaks that are off-scale are given. However, a caveat is that this could be from any of the channels (ie it could be from the ladder channel but is the same scan as the given repeat).
#'
#' ## General sample metrics
#' - `modal_peak_repeat`: The repeat size of the modal peak.
#' - `modal_peak_height`: The height of the modal peak.
#' - `index_peak_repeat`: The repeat size of the index peak (the repeat value closest to the modal peak of the index sample).
#' - `index_peak_height`: The height of the index peak.
#' - `index_weighted_mean_repeat`: The weighted mean repeat size (weighted on the height of the peaks) of the index sample.
#' - `n_peaks_total`: The total number of peaks in the repeat table.
#' - `n_peaks_analysis_subset`: The number of peaks in the analysis subset.
#' - `n_peaks_analysis_subset_expansions`: The number of expansion peaks in the analysis subset.
#' - `min_repeat`: The minimum repeat size in the analysis subset.
#' - `max_repeat`: The maximum repeat size in the analysis subset.
#' - `mean_repeat`: The mean repeat size in the analysis subset.
#' - `weighted_mean_repeat`: The weighted mean repeat size (weight on peak height) in the analysis subset.
#' - `median_repeat`: The median repeat size in the analysis subset.
#' - `max_height`: The maximum peak height in the analysis subset.
#' - `max_delta_neg`: The maximum negative delta to the index peak.
#' - `max_delta_pos`: The maximum positive delta to the index peak.
#' - `skewness`: The skewness of the repeat size distribution.
#' - `kurtosis`: The kurtosis of the repeat size distribution.
#'
#' ## Repeat instability metrics
#' - `modal_repeat_delta`: The delta between the modal peak repeat and the index peak repeat.
#' - `average_repeat_gain`: The average repeat change: The weighted mean of the sample (weighted by peak height) subtracted by the weighted mean repeat of the index sample.
#' - `instability_index`: The instability index based on peak height and distance to the index peak. (See Lee et al., 2010, https://doi.org/10.1186/1752-0509-4-29).
#' - `instability_index_abs`: The absolute instability index. The absolute value is taken for the "Change from the main allele".
#' - `expansion_index`: The instability index for expansion peaks only.
#' - `contraction_index`: The instability index for contraction peaks only.
#' - `expansion_ratio`: The ratio of expansion peaks' heights to the main peak height. Also known as "peak proportional sum" (https://doi.org/10.1016/j.cell.2019.06.036).
#' - `contraction_ratio`: The ratio of contraction peaks' heights to the main peak height.
#' - `expansion_percentile_*`: The repeat size at specified percentiles of the cumulative distribution of expansion peaks.
#' - `expansion_percentile_for_repeat_*`: The percentile rank of specified repeat sizes in the distribution of expansion peaks.

#'
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#' metadata <- trace::metadata
#'
#' test_fragments <- peak_table_to_fragments(gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B",
#'   min_size_bp = 400
#' )
#'
#' add_metadata(
#'   fragments_list = test_fragments,
#'   metadata_data.frame = metadata
#' )
#'
#' find_alleles(
#'   fragments_list = test_fragments,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1
#' )
#'
#'
#' call_repeats(
#'   fragments_list = test_fragments,
#'   repeat_calling_algorithm = "none",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' assign_index_peaks(
#'   fragments_list = test_fragments,
#'   grouped = TRUE
#' )
#'
#'
#' # grouped metrics
#' # uses t=0 samples as indicated in metadata
#' test_metrics_grouped <- calculate_instability_metrics(
#'   fragments_list = test_fragments,
#'   peak_threshold = 0.05,
#'   window_around_index_peak = c(-40, 40),
#'   percentile_range = c(0.5, 0.75, 0.9, 0.95),
#'   repeat_range = c(2, 5, 10, 20)
#' )
calculate_instability_metrics <- function(
    fragments_list,
    peak_threshold = 0.05,
    window_around_index_peak = c(NA, NA),
    percentile_range = c(0.5, 0.75, 0.9, 0.95),
    repeat_range = c(2, 5, 10, 20)) {
  # calculate metrics
  metrics_list <- lapply(fragments_list, function(fragments_repeats) {
    # check to make sure all the required steps for the function have been done
    if(fragments_repeats$.__enclos_env__$private$find_main_peaks_used == FALSE){
      stop(paste0(fragments_repeats$unique_id, " requires called alleles to calculate repeat instability metrics. Use 'find_alleles()'."),
          call. = FALSE
      )
    } 
    if(fragments_repeats$.__enclos_env__$private$assigned_index_peak_used == FALSE){
      stop(paste0(fragments_repeats$unique_id, " requires an index peak to calculate repeat instability metrics. Use 'assign_index_peaks' to set the index peaks."),
          call. = FALSE
      )
    } 

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
      allele_height = fragments_repeats$get_allele_peak()$allele_height,
      index_repeat = fragments_repeats$get_index_peak()$index_repeat,
      peak_threshold = peak_threshold,
      window_around_index_peak = window_around_index_peak
    )

    if(!is.null(fragments_repeats$.__enclos_env__$private$index_samples) && length(fragments_repeats$.__enclos_env__$private$index_samples) > 0){
      control_weighted_mean_repeat <- sapply(fragments_repeats$.__enclos_env__$private$index_samples, function(x){
        control_filtered_df <- repeat_table_subset(
          repeat_table_df = x[[2]],
          allele_height = x[[2]][which(x[[2]]$repeats == x[[1]]), "height"],
          index_repeat = x[[1]],
          peak_threshold = peak_threshold,
          window_around_index_peak = window_around_index_peak
        )

        weighted.mean(control_filtered_df$repeats, control_filtered_df$height)
      })

      index_weighted_mean_repeat <- median(control_weighted_mean_repeat, na.rm = TRUE)
    } else{
      index_weighted_mean_repeat <- NA
    }

    # first subset to make some dataframe that are just for contractions or expansions
    size_filtered_df$repeat_delta_index_peak <- size_filtered_df$repeats - fragments_repeats$get_index_peak()$index_repeat
    expansion_filtered <- size_filtered_df[which(size_filtered_df$repeat_delta_index_peak >= 0), ]
    contraction_filtered <- size_filtered_df[which(size_filtered_df$repeat_delta_index_peak <= 0), ]

    # QCs
    QC_modal_peak_height <- if (fragments_repeats$get_allele_peak()$allele_height > 500) {
      NA_character_
    } else if (fragments_repeats$get_allele_peak()$allele_height > 100) {
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
      QC_modal_peak_height = QC_modal_peak_height,
      QC_peak_number = QC_peak_number,
      QC_off_scale = QC_off_scale,
      modal_peak_repeat = fragments_repeats$get_allele_peak()$allele_repeat,
      modal_peak_height = fragments_repeats$get_allele_peak()$allele_height,
      index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
      index_peak_height = fragments_repeats$get_index_peak()$index_height,
      index_weighted_mean_repeat = index_weighted_mean_repeat,
      n_peaks_total = nrow(fragments_repeats$repeat_table_df),
      n_peaks_analysis_subset = nrow(size_filtered_df),
      n_peaks_analysis_subset_expansions = nrow(expansion_filtered),
      min_repeat = min(size_filtered_df$repeats),
      max_repeat = max(size_filtered_df$repeats),
      mean_repeat = mean(size_filtered_df$repeats),
      weighted_mean_repeat = weighted.mean(size_filtered_df$repeats, size_filtered_df$height),
      median_repeat = median(size_filtered_df$repeats),
      max_height = max(size_filtered_df$height),
      max_delta_neg = min(size_filtered_df$repeat_delta_index_peak),
      max_delta_pos = max(size_filtered_df$repeat_delta_index_peak),
      skewness = fishers_skewness(size_filtered_df$repeats, size_filtered_df$height),
      kurtosis = fishers_kurtosis(size_filtered_df$repeats, size_filtered_df$height),
      modal_repeat_delta = fragments_repeats$get_allele_peak()$allele_repeat - fragments_repeats$get_index_peak()$index_repeat,
      average_repeat_gain = weighted.mean(size_filtered_df$repeats, size_filtered_df$height) - index_weighted_mean_repeat,
      instability_index = instability_index(
        repeats = size_filtered_df$repeats,
        heights = size_filtered_df$height,
        index_peak_height = fragments_repeats$get_allele_peak()$allele_height,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ),
      instability_index_abs = instability_index(
        repeats = size_filtered_df$repeats,
        heights = size_filtered_df$height,
        index_peak_height = fragments_repeats$get_allele_peak()$allele_height,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = TRUE
      ),
      expansion_index = instability_index(
        repeats = expansion_filtered$repeats,
        heights = expansion_filtered$height,
        index_peak_height = fragments_repeats$get_allele_peak()$allele_height,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ),
      contraction_index = instability_index(
        repeats = contraction_filtered$repeats,
        heights = contraction_filtered$height,
        index_peak_height = fragments_repeats$get_allele_peak()$allele_height,
        index_peak_repeat = fragments_repeats$get_index_peak()$index_repeat,
        peak_threshold = peak_threshold,
        abs_sum = FALSE
      ),
      expansion_ratio = sum(expansion_filtered$peak_percent) - 1, # remove the main peak by subtracting 1
      contraction_ratio = sum(contraction_filtered$peak_percent) - 1
    )

    expansion_percentile <- find_percentiles(
      expansion_filtered$repeats,
      expansion_filtered$height,
      fragments_repeats$get_index_peak()$index_repeat,
      type = "percentile",
      range = percentile_range,
      col_prefix = "expansion_percentile"
    )

    expansion_repeat <- find_percentiles(
      expansion_filtered$repeats,
      expansion_filtered$height,
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

