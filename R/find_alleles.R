
#' Find Alleles
#'
#' This function identifies main allele within each fragment object.
#'
#' @param fragments_list A list of fragment objects containing peak data.
#' @param config_file The file path to a YAML file containing a full list of parameters. This provides a central place to adjust parameters for the pipeline. Use the following command to make a copy of the YAML file: `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
#' @param ... additional parameters from any of the functions in the pipeline detailed below may be passed to this function. This overwrites values in the `config_file`. These parameters include:
#'   \itemize{
#'     \item `number_of_alleles` Number of alleles to be returned for each fragment. Must either be 1 or 2. Being able to identify two alleles is for cases when you are analyze different human samples with a normal and expanded alleles and you can't do the preferred option of simply ignoring the normal allele in [find_fragments()] (eg setting the min_bp_size above the normal allele bp size). Default: `1`.
#'     \item `peak_region_size_gap_threshold` Gap threshold for identifying peak regions. The peak_region_size_gap_threshold is a parameter used to determine the maximum allowed gap between peak sizes within a peak region. Adjusting this parameter affects the size range of peaks that can be grouped together in a region. A smaller value makes it more stringent, while a larger value groups peaks with greater size differences, leading to broader peak regions that may encompass wider size ranges. Default: `6`.
#'     \item `peak_region_signal_threshold_multiplier` Multiplier for the peak signal threshold. The peak_region_signal_threshold_multiplier parameter allows adjusting the threshold for identifying peak regions based on peak signals. Increasing this multiplier value will result in higher thresholds, making it more stringent to consider peaks as part of a peak region. Conversely, reducing the multiplier value will make the criteria less strict, potentially leading to more peaks being grouped into peak regions. It's important to note that this parameter's optimal value depends on the characteristics of the data and the specific analysis goals. Choosing an appropriate value for this parameter can help in accurately identifying meaningful peak regions in the data. Default: `1`.
#'    }
#' @return This function modifies list of fragments objects in place with alleles added.
#'
#' @details This function finds the main alleles for each fragment in the list by identifying clusters of peaks ("peak regions") with the highest signal intensities. This is based on the idea that PCR amplicons of repeats have clusters of peaks (from somatic mosaicism and PCR artifacts) that help differentiate the main allele of interest from capillary electrophoresis noise/contamination.
#'
#' If number_of_alleles = 1, the tallest of peaks will be selected as the allele. This means that if your sample has multiple alleles, you have two options i) make sure that your data is subsetted to only include the allele of interest (using `min_bp_size` in [find_fragments()] to make sure that the smaller allele is excluded), or ii) setting number_of_alleles = 2, which will pick the two tallest peaks in their respective peak regions and set the main allele as the larger repeat size, and allele_2 as the shorter repeat size. We recommend the subsetting approach since that is far simpler and less likely to fail, and the second option only if you're doing an experiment analysis a large number of human samples where both the normal and expanded allele repeat lengths vary, which makes it very difficult to find a common bp size that excludes the normal allele.
#' 
#' The parameters `peak_region_signal_threshold_multiplier` and `peak_region_size_gap_threshold` will only need to be adjusted in rare cases if peaks are not being found for some reason. They influence the criteria for identifying peak regions. peak_region_signal_threshold_multiplier is multiplied to the mean height of all the peaks to create a hight threshold for inclusion into the peak region, so most of the time it's already a very low value and probably only needs to be changed if you have very few peaks. peak_region_size_gap_threshold is the distance between the peaks, either bp size, or repeats if repeats have already been called. 
#' @export
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' find_fragments(fsa_list,
#'   min_bp_size = 300
#' )
#'
#'
#' find_alleles(
#'   fsa_list,
#'   peak_region_signal_threshold_multiplier = 1
#' )
find_alleles <- function(
    fragments_list,
    config_file = NULL,
     ...) {
  # internal helper functions
  find_peak_regions <- function(signal, size) {
    
    peak_regions <- rep(NA_real_, length(signal))
    mean_signal <- mean(signal) * config$peak_region_signal_threshold_multiplier
    # loop over each fragment and check to see if it's within the thresholds
    for (i in seq_along(signal)) {
      if (signal[i] < mean_signal || i == 1 || i == length(signal)) {
        peak_regions[i] <- NA_real_
      } else if (signal[i - 1] < mean_signal && signal[i + 1] < mean_signal) {
        peak_regions[i] <- NA_real_
      } else {
        # check to see if peaks before it are within the size threshold
        current_size <- size[i]
        valid_lower_peaks <- which(size < current_size & size > current_size - config$peak_region_size_gap_threshold & signal > mean_signal)
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


  # NEED TO VALIDATE INPUTS
 
  config <- load_config(config_file, ...)

  # prepare output file
  output <- trace_output$new("find_alleles")

  lapply(fragments_list, function(fragment) {
    # the main idea here is that PCR generates clusters of peaks around the main alleles.
    # find the cluster of peaks and pick the tallest within each cluster
    # then of those clusters, pick out the tallest of them all


    # first select if working off repeat size or bp size, and return warning if going off repeats
    fragment_signal <- if (!is.null(fragment$peak_table_df)) fragment$peak_table_df$signal else fragment$repeat_table_df$signal
    fragment_sizes <- if (!is.null(fragment$peak_table_df)) fragment$peak_table_df$size else fragment$repeat_table_df$repeats
    if (is.null(fragment$peak_table_df) & config$peak_region_size_gap_threshold == 6){
      if(!any( sapply(output$warning_message, function(x) grep("Alleles were called on repeat size", x)) )){
        output$set_status(
          "warning", 
          "Alleles were called on repeat size. The default peak_region_size_gap_threshold is set expecting bp size, so may need to be decreased (eg, 6 / 3 repeats = 2 for the value in repeat units). This is probably only relevant if selecting two alleles."
        )
      }      
    }

    # Find peak regions
    peak_regions <- find_peak_regions(fragment_signal, fragment_sizes)

    # find all possible peaks
    all_peaks <- pracma::findpeaks(fragment_signal, peakpat = "[+]{1,}[0]*[-]{1,}")

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
        top_regional_peaks_positions[i] <- region_positions[which.max(fragment_signal[region_positions])][1]
      }
    }

    if(config$number_of_alleles == 1){
      top_regional_peaks_positions <-
        top_regional_peaks_positions[order(fragment_signal[top_regional_peaks_positions], decreasing = TRUE)][1]
      
      # change this so that it populates either repeat, size
      if (is.null(fragment$repeat_table_df)) {
        fragment$set_allele_peak(allele = 1, unit = "size", value = fragment_sizes[top_regional_peaks_positions])
      } else {
        fragment$set_allele_peak(allele = 1, unit = "repeats", value = fragment_sizes[top_regional_peaks_positions])
      }
      # set the second allele to NA just in case find_alleles was reran with 1 allele and 2 was used previously
      fragment$set_allele_peak(allele = 2, unit = "repeats", value = NA) 
    } 

    # below is all the code to deal with cases when the user wants to return two peaks
    # it significantly complicates the selection of peaks, but basically the two tallest peak regions will be selected
    # the main largest bp size of the options will then be set as the main allele
    # do a first pass and if only one significant peak region found when we expect two, see if there are two significant maxima in the region
    # this is for human patient data and to identify alleles close in size and homozygous alleles

    if(config$number_of_alleles == 2){
      if(length(top_regional_peaks_positions) == 1){
        region_positions <- which(peak_regions == 1)
        region_maxima <- all_peaks[which(all_peaks[, 2] %in% region_positions), 2]
        significant_maxima <-
          region_maxima[which(fragment_signal[region_maxima] > max(fragment_signal[region_positions]) * 0.5)]
        # chose the two tallest maxima if more than one peak has now been found
        if (length(significant_maxima) > 1) {
          sig_maxima_signal <- fragment_signal[significant_maxima]
          second_tallest_signal <- sig_maxima_signal[order(sig_maxima_signal, decreasing = TRUE)][2]
          top_regional_peaks_positions <- significant_maxima[which(fragment_signal[significant_maxima] >= second_tallest_signal)][1:2]
        } # deal with case where the peak after the wt peak is pretty high, perhaps indicating heterozygous +1
        else if (fragment_signal[top_regional_peaks_positions[1] + 1] / fragment_signal[top_regional_peaks_positions[1]] > 0.5) {
          top_regional_peaks_positions <- c(top_regional_peaks_positions[1], top_regional_peaks_positions[1] + 1)
        } # homozygous
        else {
          top_regional_peaks_positions <- c(top_regional_peaks_positions[1], top_regional_peaks_positions[1])
        }
      } else{
        top_regional_peaks_positions <- top_regional_peaks_positions[order(fragment_signal[top_regional_peaks_positions], decreasing = TRUE)][1:2]
        top_regional_peaks_positions <- sort(top_regional_peaks_positions)
      }
      
      # change this so that it populates either repeat, size
      if (is.null(fragment$repeat_table_df)) {
        fragment$set_allele_peak(allele = 1, unit = "size", value = fragment_sizes[top_regional_peaks_positions[2]])
        fragment$set_allele_peak(allele = 2, unit = "size", value = fragment_sizes[top_regional_peaks_positions[1]]) # smaller of the size
      } else {
        fragment$set_allele_peak(allele = 1, unit = "repeats", value = fragment_sizes[top_regional_peaks_positions[2]])
        fragment$set_allele_peak(allele = 2, unit = "repeats", value = fragment_sizes[top_regional_peaks_positions[1]]) # smaller of the size
      }
    } else if(config$number_of_alleles != 1){
      stop(call. = FALSE, "invalid 'number_of_alleles', must either be 1 or 2")
    }

    return(fragment)
  })

  # give warning for samples with no alleles called
  no_alleles <- names(fragments_list)[sapply(fragments_list, function(x) is.na(x$get_allele_peak()$allele_signal))]
  if(length(no_alleles) > 0){
    output$set_status(
      "warning", 
      paste(
        "Alleles were not called for the following samples:",
        paste0(no_alleles, collapse = ", ")
      )
    )
  }
  return(output)
}

