
#' Find Alleles
#'
#' This function identifies main allele within each fragment object.
#'
#' @param fragments_list A list of fragment objects containing peak data.
#' @param peak_region_size_gap_threshold Gap threshold for identifying peak regions. The peak_region_size_gap_threshold is a parameter used to determine the maximum allowed gap between peak sizes within a peak region. Adjusting this parameter affects the size range of peaks that can be grouped together in a region. A smaller value makes it more stringent, while a larger value groups peaks with greater size differences, leading to broader peak regions that may encompass wider size ranges.
#' @param peak_region_height_threshold_multiplier Multiplier for the peak height threshold. The peak_region_height_threshold_multiplier parameter allows adjusting the threshold for identifying peak regions based on peak heights. Increasing this multiplier value will result in higher thresholds, making it more stringent to consider peaks as part of a peak region. Conversely, reducing the multiplier value will make the criteria less strict, potentially leading to more peaks being grouped into peak regions. It's important to note that this parameter's optimal value depends on the characteristics of the data and the specific analysis goals. Choosing an appropriate value for this parameter can help in accurately identifying meaningful peak regions in the data.
#'
#' @return A list of fragments with identified main alleles.
#'
#' @details This function finds the main alleles for each fragment in the list by identifying clusters of peaks ("peak regions") with the highest signal intensities. This is based on the idea that PCR amplicons of repeats have clusters of peaks (from somatic mosaicism and PCR artifacts) that help differentiate the main allele of interest from capillary electrophoresis noise/contamination.
#'
#' The tallest of peaks will be selected as the allele. This means that if your sample has multiple alleles, you'll need to make sure that your data is subsetted to exclude that second allele. For example in a human cell line with a shorter repeat, we recommend using `min_bp_size` in [find_fragments()] to make sure that the smaller allele is excluded. Otherwise, the smaller taller allele will be set as the allele.
#' 
#' The parameters `peak_region_height_threshold_multiplier` and `peak_region_size_gap_threshold` will only need to be adjusted if peaks are not being found for some reason. They influence the criteria for identifying peak regions, and finding the right balance between them is crucial.
#' @export
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' fragments_list <- find_fragments(fsa_list,
#'   min_bp_size = 300
#' )
#'
#'
#' find_alleles(
#'   fragments_list,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1
#' )
find_alleles <- function(
    fragments_list,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1) {
  # internal helper functions
  find_peak_regions <- function(height, size) {
    peak_regions <- rep(NA_real_, length(height))
    mean_height <- mean(height) * peak_region_height_threshold_multiplier
    # loop over each fragment and check to see if it's within the thresholds
    for (i in seq_along(height)) {
      if (height[i] < mean_height || i == 1 || i == length(height)) {
        peak_regions[i] <- NA_real_
      } else if (height[i - 1] < mean_height && height[i + 1] < mean_height) {
        peak_regions[i] <- NA_real_
      } else {
        # check to see if peaks before it are within the size threshold
        current_size <- size[i]
        valid_lower_peaks <- which(size < current_size & size > current_size - peak_region_size_gap_threshold & height > mean_height)
        unique_regions <- unique(na.omit(peak_regions))
        if (length(valid_lower_peaks) > 0) {
          if (length(unique_regions) > 0) {
            peak_regions[i] <- unique_regions[length(unique_regions)]
          } else {
            peak_regions[i] <- 1
          }
        } else {
          if (length(unique_regions) > 0) {
            peak_regions[i] <- unique_regions[length(unique_regions)] + 1
          } else {
            peak_regions[i] <- 1
          }
        }
      }
    }

    return(peak_regions)
  }

  main_peaks <- lapply(fragments_list, function(fragment) {
    # the main idea here is that PCR generates clusters of peaks around the main alleles.
    # find the cluster of peaks and pick the tallest within each cluster
    # then of those clusters, pick out the tallest of them all


    # first select if working off repeat size or bp size
    fragment_height <- if (is.null(fragment$repeat_table_df)) fragment$peak_table_df$height else fragment$repeat_table_df$height
    fragment_sizes <- if (is.null(fragment$repeat_table_df)) fragment$peak_table_df$size else fragment$repeat_table_df$repeats

    # Find peak regions
    peak_regions <- find_peak_regions(fragment_height, fragment_sizes)

    # find all possible peaks
    all_peaks <- pracma::findpeaks(fragment_height, peakpat = "[+]{1,}[0]*[-]{1,}")

    # Find unique peak regions
    unique_regions <- unique(na.omit(peak_regions))
    top_regional_peaks_positions <- numeric(length(unique_regions))

    # Find the tallest peak within each peak region
    for (i in seq_along(unique_regions)) {
      region_positions <- which(peak_regions == i)

      if (any(region_positions %in% all_peaks[, 2])) {
        # Find the position of the tallest peak within the maxima positions
        peak_region_subset <- all_peaks[which(all_peaks[, 2] %in% region_positions), , drop = FALSE]
        top_regional_peaks_positions[i] <- peak_region_subset[which.max(peak_region_subset[, 1]), 2]
      } else {
        # just pick the tallest if somehow the peak region doesn't have a peak called
        # not sure if this will happen, just dealing with a possible case
        top_regional_peaks_positions[i] <- region_positions[which.max(fragment_height[region_positions])][1]
      }
    }

    # Now we need to pick the tallest of the candidates
    top_regional_peaks_positions <-
      top_regional_peaks_positions[order(fragment_height[top_regional_peaks_positions], decreasing = TRUE)][1]
    

    if (length(top_regional_peaks_positions) == 0) {
      warning(paste0(fragment$unique_id, ": No main alleles identified"))
    }

    # change this so that it populates either repeat, size
    if (is.null(fragment$repeat_table_df)) {
      fragment$set_allele_peak(unit = "size", value = fragment_sizes[top_regional_peaks_positions])
    } else {
      fragment$set_allele_peak(unit = "repeats", value = fragment_sizes[top_regional_peaks_positions])
    }

    # peak_regions
    fragment$.__enclos_env__$private$peak_regions <- peak_regions

    return(fragment)
  })

  invisible()
}

