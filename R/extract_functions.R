

#' Extract traces
#'
#' Extract the raw trace from a list of fragments objects
#'
#' @param fragments_list a list of fragments objects
#'
#' @return A dataframe of the raw trace data. Each row representing a single scan.
#' @export
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' extracted_traces <- extract_trace_table(fsa_list)
#'
extract_trace_table <- function(fragments_list) {
  # turn the output into a dataframe
  plate_list <- lapply(fragments_list, function(x) {
    x$trace_bp_df
  })

  plate_combined_df <- do.call(rbind, plate_list)

  return(plate_combined_df)
}



#' Extract ladder summary
#'
#' Extract a table summarizing the ladder models
#'
#' @param fragments_list a list of fragments trace objects
#' @param sort A logical statement for if the samples should be ordered by average ladder R-squared.
#'
#' @return a dataframe of ladder quality information
#' @export
#'
#' @details
#' The ladder peaks are assigned using a custom algorithm that maximizes the fit of detected ladder peaks and given base-pair sizes. This function summarizes the R-squared values of these individual correlations.
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
    fragments_list,
    sort = FALSE){

  # test to make sure that fragments trace objects
  if(any(sapply(fragments_list, function(x) class(x)[1] != "fragments"))){
    stop(call. = FALSE, "Wrong objects supplied. Please supply a list of 'fragments' objects")
  }
  
  summary_list <- lapply(fragments_list, function(fragment){
    cor_list <- ladder_fit_cor(fragment)
    rsq <- sapply(cor_list, function(x) x$rsq)

    data.frame(
      unique_id = fragment$unique_id,
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
#' @param fragments_list A list of fragments objects containing fragment data.
#'
#' @return A dataframe containing modal peak information for each sample
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#'
#' test_fragments <- genemapper_table_to_fragments(gm_raw,
#'   dye_channel = "B",
#'   min_size_bp = 400
#' )
#'
#' find_alleles(
#'   fragments_list = test_fragments,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_signal_threshold_multiplier = 1
#' )
#'
#' extract_alleles(test_fragments)
#'
extract_alleles <- function(fragments_list) {
  extracted <- lapply(fragments_list, function(x) {
    df <- as.data.frame(x$get_allele_peak())
    df$unique_id <- x$unique_id
    df[,c(ncol(df),1:(ncol(df)-1))]
  })
  extracted_df <- do.call(rbind, extracted)

  return(extracted_df)
}

# Extract fragments -------------------------------------------------------

#' Extract All Fragments
#'
#' Extracts peak data from each sample in a list of fragments.
#'
#' @param fragments_list A list of fragments objects containing fragment data.
#'
#' @return A dataframe containing peak data for each sample
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#' metadata <- trace::metadata
#'
#' test_fragments <- genemapper_table_to_fragments(gm_raw,
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
          main_peak_signal = rep(x$get_allele_peak()$allele_signal, df_length),
          signal = x$peak_table_df$signal,
          size = x$peak_table_df$size
        )
      } else if (!is.null(x$repeat_table_df)) {
        df_length <- nrow(x$repeat_table_df)
        data.frame(
          unique_id = rep(x$unique_id, df_length),
          main_peak_repeat = rep(x$get_allele_peak()$allele_repeat, df_length),
          main_peak_signal = rep(x$get_allele_peak()$allele_signal, df_length),
          signal = x$repeat_table_df$signal,
          repeats = x$repeat_table_df$repeats
        )
      }
    })
  )
  extracted_df <- do.call(rbind, extracted)


  return(extracted_df)
}



#' Extract repeat correction summary
#'
#' Extracts a table summarizing the model used to correct repeat length
#'
#' @param fragments_list A list of fragments class objects obtained from the [call_repeats()] function when the `correction = "repeat"` parameter is used.
#' @export
#' @return A data.frame
#' @details
#' For each of the samples used for repeat correction, this table pulls out the modal repeat length called by the model (`allele_repeat`), how far that sample is on average from the linear model in repeat units by finding the average residuals (`avg_residual`), and the absolute value of the `avg_residual` (`abs_avg_residual`)
#' 
#' @examples
#'
#'
#' fsa_list <- lapply(cell_line_fsa_list[16:19], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' find_fragments(fsa_list, min_bp_size = 300)
#'
#' test_alleles <- find_alleles(
#'   fsa_list 
#' )
#' 
#' add_metadata(
#'   fsa_list,
#'   metadata
#' )
#'
#'
#' call_repeats(
#'   fragments_list = fsa_list,
#'   correction = "repeat"
#' )
#'
#' # finally extract repeat correction summary
#' extract_repeat_correction_summary(fsa_list)
#'
#'
extract_repeat_correction_summary <- function(
  fragments_list
){
  # first do some validation to check if it's valid that they are trying to use this function
  if(is.null(fragments_list[[1]]$.__enclos_env__$private$repeat_correction_mod)){
    stop(call. = FALSE, "No repeat correction model detected in the first sample. You must have used correction = 'repeat' in call_repeats() to use this function.")
  }
  first_model_df <- fragments_list[[1]]$.__enclos_env__$private$repeat_correction_mod
  identical_model_test <- logical(length(fragments_list))
  for (i in seq_along(fragments_list)) {
    identical_model_test[i] <- identical(first_model_df, fragments_list[[i]]$.__enclos_env__$private$repeat_correction_mod)
  }
  if (!all(identical_model_test)) {
    stop("The supplied fragments list must come from the same 'call_repeats' function output", call. = FALSE)
  }

  #generate table of average residual and allele
  controls_repeats_df <- fragments_list[[1]]$.__enclos_env__$private$repeat_correction_mod$model
  controls_repeats_df$unique_id <- sub("\\.[0-9]+$", "", row.names(controls_repeats_df))
  controls_repeats_df$residuals <- fragments_list[[1]]$.__enclos_env__$private$repeat_correction_mod$residuals
  
  controls_repeats_df_split <- split(controls_repeats_df, controls_repeats_df$unique_id)
  controls_repeats_df_split_summarized <- lapply(controls_repeats_df_split, function(x){
    data.frame(
      unique_id = unique(x$unique_id),
      avg_residual = mean(x$residuals)
    )
  })
  
  controls_repeats_summarized <- do.call(rbind, controls_repeats_df_split_summarized)
  controls_repeats_summarized$batch_run_id <- sapply(fragments_list[controls_repeats_summarized$unique_id], function(x) x$batch_run_id)
  controls_repeats_summarized$batch_sample_id <- sapply(fragments_list[controls_repeats_summarized$unique_id], function(x) x$batch_sample_id)
  controls_repeats_summarized$batch_sample_modal_repeat <- sapply(fragments_list[controls_repeats_summarized$unique_id], function(x) x$batch_sample_modal_repeat)
  controls_repeats_summarized$allele_repeat <- sapply(fragments_list[controls_repeats_summarized$unique_id], function(x) x$get_allele_peak()$allele_repeat)
  controls_repeats_summarized$abs_avg_residual <- abs(controls_repeats_summarized$avg_residual)
  # reorder cols
  controls_repeats_summarized <- controls_repeats_summarized[ ,names(controls_repeats_summarized)[c(1,3:6,2,7)]]

  return(controls_repeats_summarized)
}






