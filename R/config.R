#' Load Trace Config File
#'
#' Load configuration file that can be passed to individual trace modules.
#'
#' @param config_file The file path to a YAML file containing a full list of parameters. If NULL, it will load a default YAML.
#'
#' @return A trace_config object.
#'
#' @details
#' Use this function to load in a config file for passing to individual trace modules. This is only required if creating a custom pipeline rather than the main function. This provides a central place to adjust parameters for the pipeline. Either edit the YAML directly, or modify the object (which is basically just a list).
#' 
#' Use the following command to make a copy of the YAML file: `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
#' 
#' @importFrom  yaml read_yaml
#' @export
#' 
#' @examples
#' 
#' config <- load_config()
#'
load_config <- function(config_file = NULL) {
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

  # correct type of NA for certain parameters
  na_params <- c("ladder_start_scan", "minimum_ladder_signal",
                  "min_scan", "max_scan")
  for (param in na_params) {
    if(config[[param]] == "NA"){
      config[[param]] <- NA_real_
    }
  }

  # set config as S3 class
  class(config) <- "trace_config"

  return(config)
}


 update_config <- function(config, ...){
     
  # override config file with user supplied arguments
  user_args <- list(...)

  for(arg in names(user_args)){
  config[[arg]] <- user_args[[arg]]
  }
   
  validate_inputs(config)

  return(config)
 }

validate_inputs <- function(config){

  expected <- list(
    ladder_channel = list(type = "character", length = "single", allow_na = FALSE, allow_inf = FALSE),
    signal_channel = list(type = "character", length = "single", allow_na = FALSE, allow_inf = FALSE),
    ladder_sizes = list(type = "numeric", length = "multiple", allow_na = FALSE, allow_inf = FALSE),
    ladder_start_scan = list(type = "numeric", length = "single", allow_na = TRUE, allow_inf = FALSE),
    minimum_ladder_signal = list(type = "numeric", length = "single", allow_na = TRUE, allow_inf = TRUE),
    min_scan = list(type = "numeric", length = "single", allow_na = TRUE, allow_inf = TRUE),
    max_scan = list(type = "numeric", length = "single", allow_na = TRUE, allow_inf = TRUE),
    ladder_selection_window = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    max_combinations = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    warning_rsq_threshold = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    show_progress_bar = list(type = "logical", length = "single", allow_na = FALSE, allow_inf = FALSE),
    smoothing_window = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    minimum_peak_signal = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = TRUE),
    min_bp_size = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = TRUE),
    max_bp_size = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = TRUE),
    peak_scan_ramp = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    number_of_alleles = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    peak_region_size_gap_threshold = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    peak_region_signal_threshold_multiplier = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    assay_size_without_repeat = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    repeat_size = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    correction = list(type = "character", length = "single", allow_na = FALSE, allow_inf = FALSE),
    force_whole_repeat_units = list(type = "logical", length = "single", allow_na = FALSE, allow_inf = FALSE),
    force_repeat_pattern = list(type = "logical", length = "single", allow_na = FALSE, allow_inf = FALSE),
    force_repeat_pattern_size_period = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    force_repeat_pattern_size_window = list(type = "numeric", length = "single", allow_na = FALSE, allow_inf = FALSE),
    grouped = list(type = "logical", length = "single", allow_na = FALSE, allow_inf = FALSE)
  )

  # Check for unexpected parameters
  for (param in names(config)) {
    if (!param %in% names(expected)) {
      stop(paste("Unexpected parameter:", param))
    }
  }

  # Validate each parameter
  for (param in names(expected)) {
    
    if (!param %in% names(config)) {
      stop(paste("Missing parameter:", param))
    }

    # Check for NA values
    if (any(is.na(config[[param]])) && !expected[[param]]$allow_na) {
      stop(paste("Parameter", param, "contains NA values but they are not allowed"))
    }

    # Check for Inf values
    if (any(is.infinite(config[[param]])) && !expected[[param]]$allow_inf) {
      stop(paste("Parameter", param, "contains Inf values but they are not allowed"))
    }

    # Check type (if not NA or Inf)
    if (!any(is.na(config[[param]])) && !any(is.infinite(config[[param]]))) {
      # Handle integer as a valid subtype of numeric
      if (expected[[param]]$type == "numeric" && inherits(config[[param]], "integer") ) {
        # Allow integer to pass as numeric
        next
      }

      # Check type
      if (!inherits(config[[param]], expected[[param]]$type)) {
        stop(paste("Parameter", param, "should be of type", expected[[param]]$type, "but is of type", typeof(config[[param]]) ))
      }
    }

    # Check length
    if (expected[[param]]$length == "single") {
      if (length(config[[param]]) != 1) {
        stop(paste("Parameter", param, "should be a single value"))
      }
    } else if (expected[[param]]$length == "multiple") {
      if (length(config[[param]]) < 2) {
        stop(paste("Parameter", param, "should be a vector of multiple values"))
      }
    }
  }


  invisible()

}