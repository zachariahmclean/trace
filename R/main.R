#TOODO
 # take stop calls out of modules and put into the main function to handle weather to continue or stop
 # each function returns either the object or an error
   # checkerror e = tryCatch(5 + "s" , error = function(e) e) 
  # any(class(e) == "simpleError")
 
#' Main function for sample processing
#'
#' The main function for the trace package that handles processing of samples through the pipeline ready for the calculation of repeat instability metrics.
#'
#' @param fragments_list A list of fragments objects containing fragment data, generated with either [read_fsa()], [genemapper_table_to_fragments()], or [repeat_table_to_fragments()].
#' @param input_type Type of data input, either "fsa", "fragments", or "repeats", which are detailed below.
#' @param metadata_data.frame metadata passed to [add_metadata()] for grouping samples for metrics calculations or batch correction.
#' @param index_override_dataframe
#' @param ladder_df_list
#' @param config_file A YAML file containing a full list of parameters that can be adjusted for the pipeline if many need to be changed. Use the following command to make a copy of the YAML file: `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
#' @param ... additional parameters from any of the functions in the pipeline detailed below may be passed to this function. This overwrites values in the `config_file`.
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
  input_type = "fsa",
  metadata_data.frame = NULL,
  index_override_dataframe = NULL,
  ladder_df_list = NULL,
  config_file = NULL,
  ...
){
   
  # Import config file if not supplied by user
  config <- load_config(config_file, ...)


  if(!is.null(metadata_data.frame)){
    add_metadata(
      fragments_list,
      metadata_data.frame = metadata_data.frame,
      unique_id = config$unique_id,
      metrics_group_id = config$metrics_group_id,
      metrics_baseline_control = config$metrics_baseline_control,
      batch_run_id = config$batch_run_id,
      batch_sample_id = config$batch_sample_id,
      batch_sample_modal_repeat = config$batch_sample_modal_repeat
    )
  }

  sample_processed <- switch(input_type,
    fsa = trace_fsa(
      fragments_list,
      config = config,
      index_override_dataframe = index_override_dataframe,
      ladder_df_list = ladder_df_list
    ),
    fragments = trace_fragments(
      fragments_list,
      config = config,
      index_override_dataframe = index_override_dataframe
    ),
    repeats = trace_repeats(
      fragments_list,
      config = config,
      index_override_dataframe = index_override_dataframe
    )
  )
 
  return(sample_processed)
 }


# Helper function to load the configuration file
load_config <- function(config_file, ...) {
  if (is.null(config_file)) {
    config_file <- system.file("extdata/trace_config.yaml", package = "trace")
  }

  # read in config and flatten
  config_nested <- yaml::read_yaml(config_file)
  config <- list()
  for (i in seq_along(config_nested)) {
    config <- c(config, config_nested[[i]])
  }

  # override config file with user supplied arguments
  user_args <- list(...)

  if(length(user_args) > 0){
	  for(arg in names(user_args)){
		config[[arg]] <- user_args[[arg]]
	  }
  }

  return(config)
}



## fsa pipeline
trace_fsa <-  function(x,
  config,
  index_override_dataframe,
  ladder_df_list
) {

  message("Finding ladders")

  find_ladders(
    x,
    ladder_channel = config$ladder_channel,
    signal_channel = config$signal_channel,
    ladder_sizes = config$ladder_sizes,
    ladder_start_scan = config$ladder_start_scan,
    minimum_ladder_signal = config$minimum_ladder_signal,
    scan_subset = config$scan_subset,
    ladder_selection_window = config$ladder_selection_window,
    max_combinations = config$max_combinations,
    warning_rsq_threshold = config$warning_rsq_threshold,
    show_progress_bar = config$show_progress_bar
  )

  if(!is.null(ladder_df_list)){
    fix_ladders_manual(x, ladder_df_list, config$warning_rsq_threshold)
  }

  message("Finding fragments")

  find_fragments(
    x,
    smoothing_window = config$smoothing_window,
    minimum_peak_signal = config$minimum_peak_signal,
    min_bp_size = config$min_bp_size,
    max_bp_size = config$max_bp_size,
    peakpat = config$peakpat
  )

  trace_fragments(x,
    config = config,
    index_override_dataframe = index_override_dataframe
  )

  return(x)
}

## fragments pipeline
trace_fragments <-  function(x,
  config,
  index_override_dataframe) {
  
  message("Finding alleles")

  find_alleles(
    x,
    number_of_alleles = config$number_of_alleles,
    peak_region_size_gap_threshold = config$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$peak_region_signal_threshold_multiplier
  )

  message("Calling repeats")

  call_repeats(
    x,
    assay_size_without_repeat = config$assay_size_without_repeat,
    repeat_size = config$repeat_size,
    correction = config$correction,
    force_whole_repeat_units = config$force_whole_repeat_units,
    force_repeat_pattern = config$force_repeat_pattern,
    force_repeat_pattern_size_period = config$force_repeat_pattern_size_period,
    force_repeat_pattern_size_window = config$force_repeat_pattern_size_window
  )

  message("Assigning index peaks")

  assign_index_peaks(
    x,
    grouped = config$grouped,
    index_override_dataframe = config$index_override_dataframe
  )

  return(x)
}

## repeats pipeline
trace_repeats <- function(x,
  config,
  index_override_dataframe ) {

  message("Finding alleles")
  
  # there's an issue that the default peak_region_size_gap_threshold in the config file is for fragments
  # if the user hasn't uploaded their own value, change it to 2 (for two repeats)
  if(is.null(config_file)){
    config$peak_region_size_gap_threshold <- 2
    message("overriding peak_region_size_gap_threshold to 2")
  }

  find_alleles(
    x,
    number_of_alleles = config$number_of_alleles,
    peak_region_size_gap_threshold = config$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$peak_region_signal_threshold_multiplier
  )

  message("Assigning index peaks")
  assign_index_peaks(
    x,
    grouped = config$grouped,
    index_override_dataframe = config$index_override_dataframe
  )

  return(x)
}
