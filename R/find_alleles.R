
#' Find Alleles in Fragments List
#'
#' This function identifies main alleles within each fragment in a list of fragments.
#'
#' @param fragments_list A list of fragment objects containing peak data.
#' @param number_of_peaks_to_return Number of main peaks to be returned for each fragment. Must either be 1 or 2 if a normal sized allele is also to be found.
#' @param peak_region_size_gap_threshold Gap threshold for identifying peak regions. The peak_region_size_gap_threshold is a parameter used to determine the maximum allowed gap between peak sizes within a peak region. Adjusting this parameter affects the size range of peaks that can be grouped together in a region. A smaller value makes it more stringent, while a larger value groups peaks with greater size differences, leading to broader peak regions that may encompass wider size ranges.
#' @param peak_region_height_threshold_multiplier Multiplier for the peak height threshold. The peak_region_height_threshold_multiplier parameter allows adjusting the threshold for identifying peak regions based on peak heights. Increasing this multiplier value will result in higher thresholds, making it more stringent to consider peaks as part of a peak region. Conversely, reducing the multiplier value will make the criteria less strict, potentially leading to more peaks being grouped into peak regions. It's important to note that this parameter's optimal value depends on the characteristics of the data and the specific analysis goals. Choosing an appropriate value for this parameter can help in accurately identifying meaningful peak regions in the data.
#'
#' @return A list of fragments with identified main alleles.
#'
#' @details This function finds the main alleles for each fragment in the list by identifying clusters of peaks ("peak regions")
#' with the highest signal intensities. This is based on the idea that PCR amplicons of repeats have broad peaks and PCR artifacts that help identifying the alleles.
#'  The number of peaks to be returned, and the parameters for identifying peak regions can be customized.
#'  It's important to note that both peak_region_height_threshold_multiplier and peak_region_size_gap_threshold influence the criteria for identifying peak regions, and finding the right balance between them is crucial.
#' @export
#'
#' @examples
#' file_list <- trace::cell_line_fsa_list
#'
#' test_ladders <- find_ladders(file_list[1])
#'
#' fragments_list <- find_fragments(test_ladders,
#'   min_bp_size = 300
#' )
#'
#'
#' test_alleles <- find_alleles(
#'   fragments_list = fragments_list,
#'   number_of_peaks_to_return = 1,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1
#' )
find_alleles <- function(
    fragments_list,
    number_of_peaks_to_return = 1,
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

  # candidate alleles ------------------------------------------------------------
  find_candidate_peaks <- function(height,
                                   size,
                                   peak_regions) {
    # find all possible peaks
    all_peaks <- pracma::findpeaks(height, peakpat = "[+]{1,}[0]*[-]{1,}")

    # Find unique peak regions
    unique_regions <- unique(na.omit(peak_regions))
    top_regional_peaks <- numeric(length(unique_regions))

    # Find the tallest peak within each peak region
    for (i in seq_along(unique_regions)) {
      region_positions <- which(peak_regions == i)

      if (any(region_positions %in% all_peaks[, 2])) {
        # Find the position of the tallest peak within the maxima positions
        peak_region_subset <- all_peaks[which(all_peaks[, 2] %in% region_positions), , drop = FALSE]
        top_regional_peaks[i] <- peak_region_subset[which.max(peak_region_subset[, 1]), 2]
      } else {
        # just pick the tallest if somehow the peak region doesn't have a peak called
        # not sure if this will happen, just dealing with a possible case
        top_regional_peaks[i] <- region_positions[which.max(height[region_positions])][1]
      }
    }

    # Now we need to pick the alleles
    # do a first pass and if only one significant peak region found when we expect two, see if there are two significant maxima in the region
    # this is for human patient data and to identify alleles close in size and homozygous alleles

    if (length(top_regional_peaks) == 1 &
      number_of_peaks_to_return == 2) {
      region_positions <- which(peak_regions == 1)
      region_maxima <- all_peaks[which(all_peaks[, 2] %in% region_positions), 2]
      significant_maxima <-
        region_maxima[which(height[region_maxima] > max(height[region_positions]) * 0.5)]
      # chose the two tallest maxima if more than one peak has now been found
      if (length(significant_maxima) > 1) {
        sig_maxima_height <- height[significant_maxima]
        second_tallest_height <-
          sig_maxima_height[order(sig_maxima_height, decreasing = TRUE)][2]
        top_regional_peaks <-
          significant_maxima[which(height[significant_maxima] >= second_tallest_height)][1:2]
      } # deal with case where the peak after the wt peak is pretty high, perhaps indicating heterozygous +1
      else if (height[top_regional_peaks[1] + 1] / height[top_regional_peaks[1]] > 0.5) {
        top_regional_peaks <- c(top_regional_peaks[1], top_regional_peaks[1] + 1)
      } # homozygous
      else {
        top_regional_peaks <- c(top_regional_peaks[1], top_regional_peaks[1])
      }
    } else {
      # otherwise pick out the desired number of peaks based on height

      top_regional_peaks <-
        top_regional_peaks[order(height[top_regional_peaks], decreasing = TRUE)][1:number_of_peaks_to_return]
    }


    return(top_regional_peaks)
  }


  main_peaks <- lapply(fragments_list, function(fragment) {
    # the main idea here is that PCR generates clusters of peaks around the main alleles.
    # find the cluster of peaks and pick the tallest within

    # do some checks
    if (!number_of_peaks_to_return %in% c(1, 2)) {
      stop("number_of_peaks_to_return must be 1 or 2",
        call. = FALSE
      )
    }

    # first select if working off repeat size or bp size
    fragment_height <- if (is.null(fragment$repeat_table_df)) fragment$peak_table_df$height else fragment$repeat_table_df$height
    fragment_sizes <- if (is.null(fragment$repeat_table_df)) fragment$peak_table_df$size else fragment$repeat_table_df$repeats

    # Find peak regions
    peak_regions <- find_peak_regions(fragment_height, fragment_sizes)

    top_regional_peaks <- find_candidate_peaks(
      fragment_height,
      fragment_sizes,
      peak_regions = peak_regions
    )

    if (length(top_regional_peaks) == 0) {
      warning(paste0(fragment$unique_id, ": No main alleles identified"))
    }

    # change this so that it populates either repeat, size

    if (number_of_peaks_to_return == 2) {
      # shorter repeat allele
      if (is.null(fragment$repeat_table_df)) {
        fragment$set_allele(allele = 2, unit = "size", value = fragment_sizes[top_regional_peaks[1]])
      } else {
        fragment$set_allele(allele = 2, unit = "repeats", value = fragment_sizes[top_regional_peaks[1]])
      }

      # longer repeat allele
      if (is.null(fragment$repeat_table_df)) {
        fragment$set_allele(allele = 1, unit = "size", value = fragment_sizes[top_regional_peaks[2]])
      } else {
        fragment$set_allele(allele = 1, unit = "repeats", value = fragment_sizes[top_regional_peaks[2]])
      }
    } else if (number_of_peaks_to_return == 1) {
      # longer repeat allele
      if (is.null(fragment$repeat_table_df)) {
        fragment$set_allele(allele = 1, unit = "size", value = fragment_sizes[top_regional_peaks[1]])
      } else {
        fragment$set_allele(allele = 1, unit = "repeats", value = fragment_sizes[top_regional_peaks[1]])
      }
    }

    # peak_regions
    fragment$.__enclos_env__$private$peak_regions <- peak_regions

    return(fragment)
  })

  return(main_peaks)
}



