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

  return(config)
 }