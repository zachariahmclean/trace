

#' Extract traces
#'
#' Extract the raw trace from a list of fragments objects
#'
#' @param fragments_trace_list a list of fragments objects
#'
#' @return returns a dataframe of the raw trace data. Each row representing a single scan.
#' @export
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' extracted_traces <- extract_trace_table(fsa_list)
#'
extract_trace_table <- function(fragments_trace_list) {
  # turn the output into a dataframe
  plate_list <- lapply(fragments_trace_list, function(x) {
    x$trace_bp_df
  })

  plate_combined_df <- do.call(rbind, plate_list)

  return(plate_combined_df)
}



#' Extract ladder summary
#'
#' Extract a table summarizing the ladder models
#'
#' @param fragments_trace_list a list of fragments trace objects
#' @param sort A logical statement for if the samples should be ordered by average ladder R-squared.
#'
#' @return a dataframe of ladder quality information
#' @export
#'
#' @details
#' The ladder peaks are assigned using a custom algorithm that maximizes the fit of detected ladder peaks and given base-pair sizes. The base pair is assigned using the local Southern method. Basically, for each data point, linear models are made for the lower and upper 3 size standard and the predicted sizes are averaged.
#'
#' This function summarizes the R-squared values of these individual linear models.
#'
#'
#' @examples
#'
#'   fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#'
#'   find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#'   extract_ladder_summary(fsa_list, sort = TRUE)
extract_ladder_summary <- function(
    fragments_trace_list,
    sort = FALSE){

  summary_list <- lapply(fragments_trace_list, function(x){
    rsq <- sapply(x$local_southern_mod, function(mod_list) suppressWarnings(summary(mod_list$mod)$r.squared))

    data.frame(
      unique_id = x$unique_id,
      avg_rsq = mean(rsq),
      min_rsq = min(rsq)
    )
  })

  summary_df <- do.call(rbind, summary_list)

  if(sort == TRUE){
    summary_df <- summary_df[order(summary_df$avg_rsq), ]
  }

  return(summary_df)

}


# Extract alleles -------------------------------------------------------

#' Extract Modal Peaks
#'
#' Extracts modal peak information from each sample in a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#'
#' @return A data.frame containing modal peak information for each sample.
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#'
#' test_fragments <- peak_table_to_fragments(gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B",
#'   min_size_bp = 400
#' )
#'
#' find_alleles(
#'   fragments_list = test_fragments,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1
#' )
#'
#' extract_alleles(test_fragments)
#'
extract_alleles <- function(fragments_list) {
  extracted <- lapply(fragments_list, function(x) {
    data.frame(
      unique_id = x$unique_id,
      size = x$get_allele_peak()$allele_size,
      repeats = x$get_allele_peak()$allele_repeat,
      height = x$get_allele_peak()$allele_height
    )
  })
  extracted_df <- do.call(rbind, extracted)

  return(extracted_df)
}

# Extract fragments -------------------------------------------------------

#' Extract All Fragments
#'
#' Extracts peak data from each sample in a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#'
#' @return A data.frame containing peak data for each sample.
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
#'   fragments_list = test_fragments
#' )
#'
#' call_repeats(
#'   fragments_list = test_fragments,
#'   repeat_calling_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' extract_alleles(test_fragments)
#'
extract_fragments <- function(fragments_list) {
  suppressWarnings(
    extracted <- lapply(fragments_list, function(x) {
      if (is.null(x$peak_table_df) & is.null(x$repeat_table_df)) {
        return(NULL)
      } else if (!is.null(x$peak_table_df) & is.null(x$repeat_table_df)) {
        df_length <- nrow(x$peak_table_df)
        data.frame(
          unique_id = rep(x$unique_id, df_length),
          main_peak_size = rep(x$get_allele_peak()$allele_size, df_length),
          main_peak_height = rep(x$get_allele_peak()$allele_height, df_length),
          height = x$peak_table_df$height,
          size = x$peak_table_df$size,
          peak_region = x$.__enclos_env__$private$peak_regions
        )
      } else if (!is.null(x$repeat_table_df)) {
        df_length <- nrow(x$repeat_table_df)
        data.frame(
          unique_id = rep(x$unique_id, df_length),
          main_peak_repeat = rep(x$get_allele_peak()$allele_repeat, df_length),
          main_peak_height = rep(x$get_allele_peak()$allele_height, df_length),
          height = x$repeat_table_df$height,
          repeats = x$repeat_table_df$repeats,
          peak_region = x$.__enclos_env__$private$peak_regions
        )
      }
    })
  )
  extracted_df <- do.call(rbind, extracted)


  return(extracted_df)
}
