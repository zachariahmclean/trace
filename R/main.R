#TOODO
 # take stop calls out of modules and put into the main function to handle weather to continue or stop
 # each function returns either the object or an error
   # checkerror e = tryCatch(5 + "s" , error = function(e) e) 
  # any(class(e) == "simpleError")
 
#' Main function for sample processing
#'
#' The main function for the trace package that handles processing of samples through the pipeline ready for the calculation of repeat instability metrics.
#'
#' @param fragments_list A list of fragments objects containing fragment data, generated with either [read_fsa()], [peak_table_to_fragments()], or [repeat_table_to_repeats()].
#' @param ... additional parameters from any of the functions in the pipeline detailed below may be passed to this function. This overwrites values in the `config_file`.
#' @param input_type Type of data input, either "fsa", "fragments", or "repeats", which are detailed below.
#' @param metadata_data.frame metadata passed to [add_metadata()] for grouping samples for metrics calculations or batch correction.
#' @param config_file A YAML file containing a full list of parameters that can be adjusted for the pipeline if many need to be changed. Use the following command to make a copy of the YAML file: `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
#'
#' @return A list of fragments objects ready for calculation of instability metrics using [calculate_instability_metrics()]
#'
#' @details
#' This function goes through the full pipeline applying the library of functions within this package. To adjust parameters, you can either pass them directing to this function, or you can edit a configuration file and supply it to "config_file". 
#' 
#' fsa pipeline: [add_metadata()] (only if metadata_data.frame supplied), [find_ladders()], [fix_ladders_manual()] (only if ladder_df_list is supplied), [find_fragments()], [find_alleles()], [call_repeats()], [assign_index_peaks()].
#' 
#' fragments pipeline: [add_metadata()] (only if metadata_data.frame supplied), [find_alleles()], [call_repeats()], [assign_index_peaks()].
#' 
#' repeats pipeline: [add_metadata()] (only if metadata_data.frame supplied), [find_alleles()], [assign_index_peaks()].
#' 
#' 
#' @export
#' 
#' @examples
#' 
#' # import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
#' fragments_list <- trace_main(cell_line_fsa_list, grouped = TRUE, metadata_data.frame = metadata)
#' metrics <- calculate_instability_metrics(
#'   fragments_list = fragments_list,
#'   peak_threshold = 0.05,
#'   window_around_index_peak = c(-40, 40),
#'   percentile_range = c(0.5, 0.75, 0.9, 0.95),
#'   repeat_range = c(2, 5, 10, 20)
#' )
#'
 trace_main <- function(
  fragments_list, 
  ...,
  input_type = "fsa",
  metadata_data.frame = NULL,
  config_file = NULL
){
   
  # Import config file if not supplied by user
  config <- load_config(config_file, ...)


   # I think input_type is duplicative of the method used to import fragments_list
   
   
  sample_processed <- switch(input_type,
    fsa = trace_fsa(
      fragments_list,
      config = config,
      metadata_data.frame = metadata_data.frame,
      ladder_df_list = ladder_df_list
    ),
    fragments = trace_fragments(
      fragments_list,
      config = config,
      metadata_data.frame = metadata_data.frame
    ),
    repeats = trace_repeats(
      fragments_list,
      config = config,
      metadata_data.frame = metadata_data.frame
    )
  )
 
  return(sample_processed)
 }


# Helper function to load the configuration file
load_config <- function(config_file, ...) {
  if (is.null(config_file)) {
    config_file <- system.file("extdata/trace_config.yaml", package = "trace")
  }

  config <- yaml::read_yaml(config_file)

  user_args <- list(...)

  if(length(user_args) > 0){
    for (i in seq_along(config)) {
      for (param in names(user_args)) {
        if(param %in% names(config[[i]])){
          config[[i]][[param]] <- user_args[[param]]
        }
      }   
    }
  }

  # need to add index_override_dataframe if supplied since it cant be in config file
  if(!is.null(user_args$index_override_dataframe)){
    config$assign_index_peaks$index_override_dataframe <- user_args$index_override_dataframe
  } else{
    config$assign_index_peaks$index_override_dataframe <- NULL
  }

  # add ladder fixing if requried
  if(!is.null(user_args$ladder_df_list)){
    config$fix_ladders_manual$ladder_df_list <- user_args$ladder_df_list
  } else{
    config$fix_ladders_manual$ladder_df_list <- NULL
  }

  return(config)
}



## fsa pipeline
trace_fsa <-  function(x,
  config,
  metadata_data.frame,
  index_override_dataframe,
  ladder_df_list
) {

  if(!is.null(metadata_data.frame)){
    add_metadata(
      x,
      metadata_data.frame = metadata_data.frame,
      unique_id = config$add_metadata$unique_id,
      metrics_group_id = config$add_metadata$metrics_group_id,
      metrics_baseline_control = config$add_metadata$metrics_baseline_control,
      batch_run_id = config$add_metadata$batch_run_id,
      batch_sample_id = config$add_metadata$batch_sample_id,
      batch_sample_modal_repeat = config$add_metadata$batch_sample_modal_repeat
    )
  }

  find_ladders(
    x,
    ladder_channel = config$find_ladders$ladder_channel,
    signal_channel = config$find_ladders$signal_channel,
    ladder_sizes = config$find_ladders$ladder_sizes,
    ladder_start_scan = config$find_ladders$ladder_start_scan,
    minimum_peak_signal = config$find_ladders$minimum_peak_signal,
    scan_subset = config$find_ladders$scan_subset,
    ladder_selection_window = config$find_ladders$ladder_selection_window,
    max_combinations = config$find_ladders$max_combinations,
    warning_rsq_threshold = config$find_ladders$warning_rsq_threshold,
    show_progress_bar = config$find_ladders$show_progress_bar
  )

  if(!is.null(config$fix_ladders_manual$ladder_df_list)){
    fix_ladders_manual(x, config$fix_ladders_manual$ladder_df_list, config$find_ladders$warning_rsq_threshold)
  }

  x <- find_fragments(
    x,
    smoothing_window = config$find_fragments$smoothing_window,
    minimum_peak_signal = config$find_fragments$minimum_peak_signal,
    min_bp_size = config$find_fragments$min_bp_size,
    max_bp_size = config$find_fragments$max_bp_size,
    peakpat = config$find_fragments$peakpat
  )

 find_alleles(
    x,
    number_of_alleles = config$find_alleles$number_of_alleles,
    peak_region_size_gap_threshold = config$find_alleles$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$find_alleles$peak_region_signal_threshold_multiplier
  )

  call_repeats(
    x,
    assay_size_without_repeat = config$call_repeats$assay_size_without_repeat,
    repeat_size = config$call_repeats$repeat_size,
    correction = config$call_repeats$correction,
    force_whole_repeat_units = config$call_repeats$force_whole_repeat_units,
    force_repeat_pattern = config$call_repeats$force_repeat_pattern,
    force_repeat_pattern_size_period = config$call_repeats$force_repeat_pattern_size_period,
    force_repeat_pattern_size_window = config$call_repeats$force_repeat_pattern_size_window
  )

  assign_index_peaks(
    x,
    grouped = config$assign_index_peaks$grouped,
    index_override_dataframe = config$assign_index_peaks$index_override_dataframe
  )

  return(x)
}

## fragments pipeline
trace_fragments <-  function(x,
  config,
  metadata_data.frame) {
  if(!is.null(metadata_data.frame)){
    add_metadata(
      x,
      metadata_data.frame = metadata_data.frame,
      unique_id = config$add_metadata$unique_id,
      metrics_group_id = config$add_metadata$metrics_group_id,
      metrics_baseline_control = config$add_metadata$metrics_baseline_control,
      batch_run_id = config$add_metadata$batch_run_id,
      batch_sample_id = config$add_metadata$batch_sample_id,
      batch_sample_modal_repeat = config$add_metadata$batch_sample_modal_repeat
    )
  }

  find_alleles(
    x,
    number_of_alleles = config$find_alleles$number_of_alleles,
    peak_region_size_gap_threshold = config$find_alleles$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$find_alleles$peak_region_signal_threshold_multiplier
  )

  call_repeats(
    x,
    assay_size_without_repeat = config$call_repeats$assay_size_without_repeat,
    repeat_size = config$call_repeats$repeat_size,
    correction = config$call_repeats$correction,
    force_whole_repeat_units = config$call_repeats$force_whole_repeat_units,
    force_repeat_pattern = config$call_repeats$force_repeat_pattern,
    force_repeat_pattern_size_period = config$call_repeats$force_repeat_pattern_size_period,
    force_repeat_pattern_size_window = config$call_repeats$force_repeat_pattern_size_window
  )

  assign_index_peaks(
    x,
    grouped = config$assign_index_peaks$grouped,
    index_override_dataframe = config$assign_index_peaks$index_override_dataframe
  )

  return(x)
}

## repeats pipeline
trace_repeats <- function(x,
  config,
  metadata_data.frame ) {
  if(!is.null(metadata_data.frame)){
    add_metadata(
      x,
      metadata_data.frame = metadata_data.frame,
      unique_id = config$add_metadata$unique_id,
      metrics_group_id = config$add_metadata$metrics_group_id,
      metrics_baseline_control = config$add_metadata$metrics_baseline_control,
      batch_run_id = config$add_metadata$batch_run_id,
      batch_sample_id = config$add_metadata$batch_sample_id,
      batch_sample_modal_repeat = config$add_metadata$batch_sample_modal_repeat
    )
  }

  message("Finding alleles")
  
  # there's an issue that the default peak_region_size_gap_threshold in the config file is for fragments
  # if the user hasn't uploaded their own value, change it to 2 (for two repeats)
  if(is.null(config_file)){
    config$find_alleles$peak_region_size_gap_threshold <- 2
    message("overriding peak_region_size_gap_threshold to 2")
  }

  find_alleles(
    x,
    number_of_alleles = config$find_alleles$number_of_alleles,
    peak_region_size_gap_threshold = config$find_alleles$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$find_alleles$peak_region_signal_threshold_multiplier
  )

  message("Assigning index peaks")
  assign_index_peaks(
    x,
    grouped = config$assign_index_peaks$grouped,
    index_override_dataframe = config$assign_index_peaks$index_override_dataframe
  )

  return(x)
}
