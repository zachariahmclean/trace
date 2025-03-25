#' Main function for sample processing
#'
#' The main function for the trace package that handles processing of samples through the pipeline ready for the calculation of repeat instability metrics.
#'
#' @param fragments_list A list of fragments objects containing fragment data, generated with either [read_fsa()], [size_table_to_fragments()], [genemapper_table_to_fragments()], or [repeat_table_to_fragments()].
#' @param metadata_data.frame metadata passed to [add_metadata()] for grouping samples for metrics calculations or batch correction.
#' @param index_override_dataframe A data.frame to manually set index peaks. Column 1: unique sample IDs, Column 2: desired index peaks (the order of the columns is important since the information is pulled by column position rather than column name). Closest peak in each sample is selected so the number needs to just be approximate. Default: `NULL`. See [assign_index_peaks()].
#' @param ladder_df_list A list of dataframes, with the names being the unique id and the value being a dataframe. The dataframe has two columns, size (indicating the bp of the standard) and scan (the scan value of the ladder peak). It's critical that the element name in the list is the unique id of the sample. Either manually figure out what scan the ladder peaks should be and generate the list, or use [fix_ladders_interactive()] to interactively generate the ladder_df_list.
#' @param config_file The file path to a YAML file containing a full list of parameters. This provides a central place to adjust parameters for the pipeline. Use the following command to make a copy of the YAML file: `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
#' @param ... additional parameters from any of the functions in the pipeline detailed below may be passed to this function. This overwrites values in the `config_file`. These parameters are grouped by functionality:
#'
#' **Ladder-related parameters:**
#'   \itemize{
#'     \item `ladder_channel`: string, which channel in the fsa file contains the ladder signal. Default: `"DATA.105"`.
#'     \item `signal_channel`: string, which channel in the fsa file contains the data signal. Default: `"DATA.1"`.
#'     \item `ladder_sizes`: numeric vector, bp sizes of ladder used in fragment analysis. Default: `c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)`.
#'     \item `ladder_start_scan`: single numeric indicating the scan number to start looking for ladder peaks. Usually this can be automatically found (when set to NA). Default: `NA`.
#'     \item `minimum_ladder_signal`: single numeric for minimum signal of peak from the raw signal. Default: `NA`.
#'     \item `min_scan`: single numeric indicating the lower scan limit to filter out scans below. Default: `NA`.
#'     \item `max_scan`: single numeric indicating the upper scan limit to filter out scans above. Default: `NA`.
#'     \item `ladder_selection_window`: single numeric for the ladder assigning algorithm. We iterate through the scans in blocks and test their linear fit (We can assume that the ladder is linear over a short distance). This value defines how large that block of peaks should be. Default: `5`.
#'     \item `max_combinations`: single numeric indicating what is the maximum number of ladder combinations that should be tested. Default: `2500000`.
#'     \item `warning_rsq_threshold`: single numeric for the value for which this function will warn you when parts of the ladder have R-squared values below the specified threshold. Default: `0.998`.
#'     \item `show_progress_bar`: single logical for showing progress bar. Default: `TRUE`.
#'   }
#'
#' **Peak-finding parameters:**
#'   \itemize{
#'     \item `smoothing_window`: Single numeric value for signal smoothing window size passed to pracma::savgol(). Default: `21`.
#'     \item `minimum_peak_signal`: Single numeric value for the minimum peak signal for a valid peak. To have no minimum signal set as "-Inf". Default: `20`
#'     \item `min_bp_size`: Single numeric value for minimum bp size of peaks to consider. Default: `200`.
#'     \item `max_bp_size`: Single numeric value for maximum bp size of peaks to consider. Default: `1000`.
#'     \item `peak_scan_ramp`: Single numeric value to indicate how many scans (increasing in signal) should be either side of the peak maxima. Default: `5`.
#'   }
#'
#' **Allele-calling parameters:**
#'   \itemize{
#'     \item `number_of_alleles`: Number of alleles to be returned for each fragment. Must either be 1 or 2. Default: `1`.
#'     \item `peak_region_size_gap_threshold`: Gap threshold for identifying peak regions. Default: `6`.
#'     \item `peak_region_signal_threshold_multiplier`: Multiplier for the peak signal threshold. Default: `1`.
#'   }
#'
#' **Repeat-calling parameters:**
#'   \itemize{
#'     \item `assay_size_without_repeat`: An integer specifying the assay size without repeat for repeat calling. Default: `87`.
#'     \item `repeat_size`: An integer specifying the repeat size for repeat calling. Default: `3`.
#'     \item `correction`: A character vector of either "batch" to carry out a batch correction from common samples across runs (known repeat length not required), or "repeat" to use samples with validated modal repeat lengths to correct the repeat length. Default: `"none"`.
#'     \item `force_whole_repeat_units`: A logical value specifying if the peaks should be forced to be whole repeat units apart. Default: `FALSE`.
#'     \item `force_repeat_pattern`: A logical value specifying if the peaks should be re-called to fit the specific repeat unit pattern. Default: `FALSE`.
#'     \item `force_repeat_pattern_size_period`: A numeric value to set the peak periodicity bp size. Default: `2.79`.
#'     \item `force_repeat_pattern_size_window`: A numeric value for the size window when assigning the peak. Default: `0.5`.
#'   }
#'
#' **Index peak assignment parameters:**
#'   \itemize{
#'     \item `grouped`: Logical value indicating whether samples should be grouped to share a common index peak. Default: `FALSE`.
#'   }
#'
#' @return A list of fragments objects ready for calculation of instability metrics using [calculate_instability_metrics()]
#'
#' @details
#' This function processes samples through the full pipeline, applying the library of functions within this package. Parameters can be adjusted either by passing them directly to this function or by editing a configuration file and supplying it to `config_file`. Below is a breakdown of the pipeline stages and their key functionalities:
#'
#' **Ladder-related parameters:**
#' The ladder is used to calibrate the base pair (bp) sizes of fragments. The ladder peaks are identified in the ladder channel using \code{\link{find_ladders}}, and a generalized additive model (GAM) with cubic regression splines is used to fit the relationship between scan numbers and bp sizes. This allows for accurate interpolation of bp sizes for all scans. Manual inspection of ladder assignments is recommended to ensure correctness. If ladders are broken, first try and adjust parameters. In a last resort, ladders may be manually set by using the `ladder_df_list` parameter. To generate this list, use the helper shiny app \code{\link{fix_ladders_interactive}}.
#'
#' **Peak-finding parameters:**
#' Fragment peaks are identified in the continuous trace data using \code{\link{find_fragments}}. The signal is smoothed using a Savitzky-Golay filter, and peaks are detected based on their signal intensity and bp size. Parameters such as `min_bp_size` and `max_bp_size` allow filtering peaks outside the desired range, while `peak_scan_ramp` controls the number of scans around the peak maxima.
#'
#' **Allele-calling parameters:**
#' The main alleles within each fragment are identified by clustering peaks into "peak regions" using \code{\link{find_alleles}}. The tallest peak in each region is selected as the main allele. If `number_of_alleles` is set to 2, the two tallest peaks in their respective regions are selected, with the larger repeat size assigned as the main allele.
#'
#' **Repeat-calling parameters:**
#' Repeat lengths are calculated using \code{\link{call_repeats}}, based on the assay size without the repeat and the specified repeat size. Options include batch correction to account for run-to-run variability and repeat correction using samples with known modal repeat lengths. The `force_whole_repeat_units` option ensures repeat lengths are whole numbers, while `force_repeat_pattern` re-calls peaks to fit the expected repeat pattern.
#'
#' **Index peak assignment parameters:**
#' The index peak is the reference repeat length used for instability metrics calculations. It is typically the inherited repeat length or the modal repeat length at a starting time point. Samples can be grouped to share a common index peak using \code{\link{assign_index_peaks}}, or the index peak can be manually overridden using `index_override_dataframe`.
#'
#' **Metadata:**
#' Metadata is used to group samples for metrics calculations, batch correction, or repeat length validation. Use \code{\link{add_metadata}} to add metadata to the fragments list. Required (even if the column is left blank) metadata fields include:
#'   - `batch_run_id`: Groups samples by fragment analysis run.
#'   - `batch_sample_id`: Links samples across batches for batch or repeat correction.
#'   - `batch_sample_modal_repeat`: Specifies the validated repeat length for samples used in repeat correction.
#'   - `metrics_group_id`: Groups samples for shared index peak assignment.
#'   - `metrics_baseline_control`: Identifies samples used as baseline controls for index peak assignment.
#'
#' @export
#' 
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#' 
#' # import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
#' fragments_list <- trace(fsa_list, grouped = TRUE, metadata_data.frame = metadata)
#' 
#' metrics <- calculate_instability_metrics(
#'   fragments_list,
#'   peak_threshold = 0.05,
#'   window_around_index_peak = c(-40, 40),
#'   percentile_range = c(0.5, 0.75, 0.9, 0.95),
#'   repeat_range = c(2, 5, 10, 20)
#' )
 trace <- function(
  fragments_list, 
  metadata_data.frame = NULL,
  index_override_dataframe = NULL,
  ladder_df_list = NULL,
  config_file = NULL,
  ...
){
   
  # copy each object to make sure that they are not modified in place
  # otherwise sometimes repeated running of pipeline results in unexpected errors
  fragments_list <- lapply(fragments_list, function(x) x$clone())
   
  # Import config file if not supplied by user
  config <- load_config(config_file)
  config <- update_config(config, ...)
   
  # set input type
  input_type <- sapply(fragments_list, function(x) x$input_method)
  input_type <- unique(input_type)
  if(length(input_type) > 1){
    stop(call. = FALSE, "'fragments_list' must be imported using the same method generated with either read_fsa(), size_table_to_fragments(), genemapper_table_to_fragments(), or repeat_table_to_fragments().")
  } 

  if(!is.null(metadata_data.frame)){

    add_metadata_status <- add_metadata(
      fragments_list,
      metadata_data.frame = metadata_data.frame
    )
    add_metadata_status$print_status()
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

## fsa pipeline
trace_fsa <-  function(x,
  config,
  index_override_dataframe,
  ladder_df_list
) {
  message("Finding ladders")
  find_ladders_status <- find_ladders(x, config)
  find_ladders_status$print_status()

  if(!is.null(ladder_df_list)){
    fix_ladders_manual_status <- fix_ladders_manual(x, ladder_df_list, config$warning_rsq_threshold)
    fix_ladders_manual_status$print_status()
  }

  message("Finding fragments")
  find_fragments_status <- find_fragments(x, config)
  find_fragments_status$print_status()

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
  find_alleles_status$print_status()

  message("Calling repeats")
  call_repeats_status <- call_repeats(x, config)
  call_repeats_status$print_status()

  message("Assigning index peaks")
  assign_index_peaks_status <- assign_index_peaks(x, config, index_override_dataframe = index_override_dataframe)
  assign_index_peaks_status$print_status()

  return(x)
}

## repeats pipeline
trace_repeats <- function(x,
  config,
  index_override_dataframe ) {

  message("Finding alleles")
  
  # there's an issue that the default peak_region_size_gap_threshold in the config file is for fragments
  # if the user hasn't uploaded their own value, give warning
  if(config$peak_region_size_gap_threshold == 6){
    config$peak_region_size_gap_threshold <- 2
    message("WARNING! peak_region_size_gap_threshold is at the default value of 6 which is likely not appropriate when starting with the data here that already has the repeats called. If analyzing triplet repeats, a more appropriate value may be 'peak_region_size_gap_threshold = 2'.")
  }

  find_alleles_status <- find_alleles(x, config)
  find_alleles_status$print_status()

  message("Assigning index peaks")
  assign_index_peaks_status <- assign_index_peaks(x, config, index_override_dataframe = index_override_dataframe)
  assign_index_peaks_status$print_status()

  return(x)
}
