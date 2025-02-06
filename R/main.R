trace_main <- function(
  fragments_list, 
  metadata_data.frame = NULL,
  ladder_channel = "DATA.105",
  signal_channel = "DATA.1",
  ladder_sizes = c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
  ladder_start_scan = NULL,
  min_bp_size = 100,
  number_of_alleles = 1,
  assay_size_without_repeat = 87,
  repeat_size = 3,
  correction = "none",
  force_whole_repeat_units = FALSE,
  force_repeat_pattern = FALSE,
  force_repeat_pattern_size_period = repeat_size * 0.93,
  force_repeat_pattern_size_window = 0.5,
  grouped = FALSE,
  index_override_dataframe = NULL,
  config_file = NULL
){

  # Determine the class of the fragments_list
  sample_import <- determine_class(fragments_list)

  # Import config file if not supplied by user
  config <- load_config(config_file)

  # # determine the exact class of the fragments_list and assign the appropriate S4 class
  # unique_class <- unique(sapply(fragments_list, function(x) class(x)[1]))
  # if(length(unique_class) > 1 | !unique_class %in% c("fragments_trace", "fragments_repeats")){
  #   stop(call. = FALSE,
  #   "fragments_list must be a list of objects generated using either read_fsa(), peak_table_to_fragments(), or repeat_table_to_repeats().")
  # } else if(unique_class == "fragments_trace"){
  #   sample_import <- new("traceFsa")
  #   sample_import@fragments_list <- fragments_list
  # } else{
  #   contains_peak_table_df <- sapply(fragments_list, function(x) !is.null(x$peak_table_df))
  #   contains_repeat_table_df <- sapply(fragments_list, function(x) !is.null(x$repeat_table_df))

  #   if(all(contains_repeat_table_df & !contains_peak_table_df)){
  #     sample_import <- new("traceRepeats")
  #     sample_import@fragments_list <- fragments_list
  #   } else if(all(contains_peak_table_df)){
  #     sample_import <- new("traceFragments")
  #     # this is problematic because when if somebody runs their fragments sample, then it adds the repeats, modifies in place, so clone here
  #     sample_import@fragments_list <- lapply(fragments_list, function(x) x$clone())
  #   } else{
  #     stop(call. = FALSE,
  #       "Import error. Not all samples have peak_table_df or repeat_table_df.")
  #   }
  # }

  # # import config file if not supplied by user
  # if(is.null(config_file)){
  #   config_file <- system.file("extdata/trace_config.yaml", package = "trace")
  # }
  # config <- yaml::read_yaml(config_file)

  sample_processed <- traceMain(
    sample_import,
    config = config,
    metadata_data.frame = metadata_data.frame,
    ladder_channel = ladder_channel,
    signal_channel = signal_channel,
    ladder_sizes = ladder_sizes,
    ladder_start_scan = ladder_start_scan,
    min_bp_size = min_bp_size,
    number_of_alleles = number_of_alleles,
    assay_size_without_repeat = assay_size_without_repeat,
    repeat_size = repeat_size,
    correction = correction,
    force_whole_repeat_units = force_whole_repeat_units,
    force_repeat_pattern = force_repeat_pattern,
    force_repeat_pattern_size_period = force_repeat_pattern_size_period,
    force_repeat_pattern_size_window = force_repeat_pattern_size_window,
    grouped = grouped,
    index_override_dataframe = index_override_dataframe,
    config_file = config_file)

  return(sample_processed@fragments_list)
  
}

# Helper function to determine the class of the fragments_list
determine_class <- function(fragments_list) {
  unique_class <- unique(sapply(fragments_list, function(x) class(x)[1]))
  
  if (length(unique_class) > 1 || !unique_class %in% c("fragments_trace", "fragments_repeats")) {
    stop("fragments_list must be a list of objects generated using either read_fsa(), peak_table_to_fragments(), or repeat_table_to_repeats().")
  }
  
  if (unique_class == "fragments_trace") {
    sample_import <- new("traceFsa")
  } else {
    contains_peak_table_df <- sapply(fragments_list, function(x) !is.null(x$peak_table_df))
    contains_repeat_table_df <- sapply(fragments_list, function(x) !is.null(x$repeat_table_df))
    
    if (all(contains_repeat_table_df & !contains_peak_table_df)) {
      sample_import <- new("traceRepeats")
    } else if (all(contains_peak_table_df)) {
      sample_import <- new("traceFragments")
      sample_import@fragments_list <- lapply(fragments_list, function(x) x$clone())
    } else {
      stop("Import error. Not all samples have peak_table_df or repeat_table_df.")
    }
  }
  
  sample_import@fragments_list <- fragments_list
  return(sample_import)
}

# Helper function to load the configuration file
load_config <- function(config_file) {
  if (is.null(config_file)) {
    config_file <- system.file("extdata/trace_config.yaml", package = "trace")
  }
  return(yaml::read_yaml(config_file))
}


# set up s4 classes
traceMainClass <- setClass("traceMainClass", 
  slots = c(
    fragments_list = "list"
  )
)

traceFsa <- setClass("traceFsa", contains = "traceMainClass")
traceFragments <- setClass("traceFragments", contains = "traceMainClass")
traceRepeats <- setClass("traceRepeats", contains = "traceMainClass")

# main function
setGeneric("traceMain", function(x, ...) standardGeneric("traceMain"))

## fsa pipeline
setMethod("traceMain", "traceFsa", function(x,
  config,
  metadata_data.frame,
  ladder_channel,
  signal_channel,
  ladder_sizes,
  ladder_start_scan,
  min_bp_size,
  number_of_alleles,
  assay_size_without_repeat,
  repeat_size ,
  correction ,
  force_whole_repeat_units ,
  force_repeat_pattern ,
  force_repeat_pattern_size_period ,
  force_repeat_pattern_size_window ,
  grouped ,
  index_override_dataframe ,
  config_file 
) {

  if(!is.null(metadata_data.frame)){
    add_metadata(
      x@fragments_list,
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
    x@fragments_list,
    ladder_channel = ladder_channel,
    signal_channel = signal_channel,
    ladder_sizes = ladder_sizes,
    ladder_start_scan =ladder_start_scan,
    minimum_peak_signal = config$find_ladders$minimum_peak_signal,
    scan_subset = config$find_ladders$scan_subset,
    ladder_selection_window = config$find_ladders$ladder_selection_window,
    max_combinations = config$find_ladders$max_combinations,
    warning_rsq_threshold = config$find_ladders$warning_rsq_threshold,
    show_progress_bar = config$find_ladders$show_progress_bar
  )

  x@fragments_list <- find_fragments(
    x@fragments_list,
    smoothing_window = config$find_fragments$smoothing_window,
    minimum_peak_signal = config$find_fragments$minimum_peak_signal,
    min_bp_size = min_bp_size,
    max_bp_size = config$find_fragments$max_bp_size,
    peakpat = config$find_fragments$peakpat
  )

 find_alleles(
    x@fragments_list,
    number_of_alleles = number_of_alleles,
    peak_region_size_gap_threshold = config$find_alleles$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$find_alleles$peak_region_size_gap_threshold
  )

  call_repeats(
    x@fragments_list,
    assay_size_without_repeat = assay_size_without_repeat,
    repeat_size = repeat_size,
    correction = correction,
    force_whole_repeat_units = force_whole_repeat_units,
    force_repeat_pattern = force_repeat_pattern,
    force_repeat_pattern_size_period = force_repeat_pattern_size_period,
    force_repeat_pattern_size_window = force_repeat_pattern_size_window
  )

  assign_index_peaks(
    x@fragments_list,
    grouped = grouped,
    index_override_dataframe = index_override_dataframe
  )

  return(x)
})

## fragments pipeline
setMethod("traceMain", "traceFragments", function(x,
  config,
  metadata_data.frame,
  ladder_channel,
  signal_channel,
  ladder_sizes,
  ladder_start_scan,
  min_bp_size,
  number_of_alleles,
  assay_size_without_repeat,
  repeat_size ,
  correction ,
  force_whole_repeat_units ,
  force_repeat_pattern ,
  force_repeat_pattern_size_period ,
  force_repeat_pattern_size_window ,
  grouped ,
  index_override_dataframe ,
  config_file ) {
  if(!is.null(metadata_data.frame)){
    add_metadata(
      x@fragments_list,
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
    x@fragments_list,
    number_of_alleles = number_of_alleles,
    peak_region_size_gap_threshold = config$find_alleles$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$find_alleles$peak_region_signal_threshold_multiplier
  )

  call_repeats(
    x@fragments_list,
    assay_size_without_repeat = assay_size_without_repeat,
    repeat_size = repeat_size,
    correction = correction,
    force_whole_repeat_units = force_whole_repeat_units,
    force_repeat_pattern = force_repeat_pattern,
    force_repeat_pattern_size_period = force_repeat_pattern_size_period,
    force_repeat_pattern_size_window = force_repeat_pattern_size_window
  )

  assign_index_peaks(
    x@fragments_list,
    grouped = grouped,
    index_override_dataframe = index_override_dataframe
  )

  return(x)
})

## repeats pipeline
setMethod("traceMain", "traceRepeats", function(x,
  config,
  metadata_data.frame,
  ladder_channel,
  signal_channel,
  ladder_sizes,
  ladder_start_scan,
  min_bp_size,
  number_of_alleles,
  assay_size_without_repeat,
  repeat_size ,
  correction ,
  force_whole_repeat_units ,
  force_repeat_pattern ,
  force_repeat_pattern_size_period ,
  force_repeat_pattern_size_window ,
  grouped ,
  index_override_dataframe ,
  config_file ) {
  if(!is.null(metadata_data.frame)){
    add_metadata(
      x@fragments_list,
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
    x@fragments_list,
    number_of_alleles = number_of_alleles,
    peak_region_size_gap_threshold = config$find_alleles$peak_region_size_gap_threshold,
    peak_region_signal_threshold_multiplier = config$find_alleles$peak_region_size_gap_threshold
  )

  message("Assigning index peaks")
  assign_index_peaks(
    x@fragments_list,
    grouped = grouped,
    index_override_dataframe = index_override_dataframe
  )

  return(x)
})
