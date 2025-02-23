#TOODO
 # take stop calls out of modules and put into the main function to handle weather to continue or stop
 # each function returns either the object or an error
   # checkerror e = tryCatch(5 + "s" , error = function(e) e) 
  # any(class(e) == "simpleError")
 
#' Main function for sample processing
#'
#' The main function for the trace package that handles processing of samples through the pipeline ready for the calculation of repeat instability metrics.
#'
#' @param fragments_list A list of fragments objects containing fragment data, generated with either [read_fsa()], [size_table_to_fragments()], [genemapper_table_to_fragments()], or [repeat_table_to_fragments()].
#' @param metadata_data.frame metadata passed to [add_metadata()] for grouping samples for metrics calculations or batch correction.
#' @param index_override_dataframe index_override_dataframe A data.frame to manually set index peaks. See [assign_index_peaks()].
#' @param ladder_df_list a list of dataframes, with the names being the unique id
#'                       and the value being a dataframe. The dataframe has two columns, size (indicating
#'                       the bp of the standard) and scan (the scan value of the ladder peak). It's
#'                       critical that the element name in the list is the unique id of the sample. 
#'                       Either manually figure out what scan the ladder peaks should be and generate the list, or use [fix_ladders_interactive()] to interactively generate the ladder_df_list.
#' @param config_file The file path to a YAML file containing a full list of parameters. This provides a central place to adjust parameters for the pipeline. Use the following command to make a copy of the YAML file: `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
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
#' @importFrom  yaml read_yaml
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
  metadata_data.frame = NULL,
  index_override_dataframe = NULL,
  ladder_df_list = NULL,
  config_file = NULL,
  ...
){
   
  # Import config file if not supplied by user
  config <- load_config(config_file, ...)
   
  # set input type
  input_type <- sapply(fragments_list, function(x) x$input_method)
  input_type <- unique(input_type)
  if(length(input_type) > 1){
    stop(call. = FALSE, "'fragments_list' must be imported using the same method. Use either generated with either read_fsa(), size_table_to_fragments(), genemapper_table_to_fragments(), or [repeat_table_to_fragments().")
  }

  if(!is.null(metadata_data.frame)){

    add_metadata_status <- add_metadata(
      fragments_list,
      metadata_data.frame = metadata_data.frame
    )
    print(add_metadata_status)
  }

  sample_processed <- switch(input_type,
    fsa = trace_fsa(
      fragments_list,
      config = config,
      index_override_dataframe = index_override_dataframe,
      ladder_df_list = ladder_df_list
    ),
    size = trace_fragments(
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

  # check if config file has already been established previously as is 'trace_config' class and return early. 
  # This is for use in main function
  if("trace_config" %in% class(config_file)){
    config <- config_file
  } else{
    # read in default config if not supplied
    if (is.null(config_file)) {
      config_file <- system.file("extdata/trace_config.yaml", package = "trace")
    }

    # read in config and flatten
    config_nested <- yaml::read_yaml(config_file)
    config <- list()
    for (i in seq_along(config_nested)) {
      config <- c(config, config_nested[[i]])
    }
  }
  
  # override config file with user supplied arguments
  user_args <- list(...)

  if(length(user_args) > 0){
	  for(arg in names(user_args)){
		config[[arg]] <- user_args[[arg]]
	  }
  }

  # set config as S3 class
  class(config) <- "trace_config"

  return(config)
}



## fsa pipeline
trace_fsa <-  function(x,
  config,
  index_override_dataframe,
  ladder_df_list
) {
  message("Finding ladders")

  find_ladders_status <- find_ladders(x, config)
  print(find_ladders_status)

  if(!is.null(ladder_df_list)){
    fix_ladders_manual_status <- fix_ladders_manual(x, ladder_df_list, config$warning_rsq_threshold)
    print(fix_ladders_manual_status)
  }

  message("Finding fragments")

  find_fragments_status <- find_fragments(x, config)
  print(find_fragments_status)

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

  find_alleles_status <- find_alleles(x,config)
  print(find_alleles_status)

  message("Calling repeats")

  call_repeats(x, config)

  message("Assigning index peaks")

  assign_index_peaks(x, config, index_override_dataframe = config$index_override_dataframe)

  return(x)
}

## repeats pipeline
trace_repeats <- function(x,
  config,
  index_override_dataframe ) {

  message("Finding alleles")
  
  # # there's an issue that the default peak_region_size_gap_threshold in the config file is for fragments
  # # if the user hasn't uploaded their own value, change it to 2 (for two repeats)
  # if(is.null(config_file)){
  #   config$peak_region_size_gap_threshold <- 2
  #   message("overriding peak_region_size_gap_threshold to 2")
  # }

  find_alleles_status <- find_alleles(x, config)
  print(find_alleles_status)

  message("Assigning index peaks")

  assign_index_peaks(x, config, index_override_dataframe = config$index_override_dataframe)

  return(x)
}
