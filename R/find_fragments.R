# find fragments ------------------------------------------------------------




#' Find fragment peaks
#'
#' Find fragment peaks in continuous trace data and convert to fragments
#' class.
#'
#' @param fragments_list A list of fragments objects containing fragment data.
#' @param config A trace_config object generated using [load_config()].
#' @param ... additional parameters from any of the functions in the pipeline detailed below may be passed to this function. This overwrites values in the `config`. These parameters include:
#'   \itemize{
#'    \item `smoothing_window` numeric, signal smoothing window size passed to pracma::savgol(). Default: `21`.
#'    \item `minimum_peak_signal` numeric, minimum signal of the raw trace. To have no minimum signal set as "-Inf". Default: `20`.
#'    \item `min_bp_size` numeric, minimum bp size of peaks to consider. Default: `100`.
#'    \item `max_bp_size` numeric, maximum bp size of peaks to consider. Default: `1000`.
#'    \item `peak_scan_ramp` Single numeric value to indicate how many scans (increasing in signal) should be either side of the peak maxima. Default: `5`.
#'  }
#'
#' @return a list of fragments objects.
#' @export
#'
#' @importFrom pracma findpeaks
#' @importFrom pracma savgol
#'
#' @details
#' 
#' This takes in a list of fragments objects and returns a list of new fragments objects.
#' 
#' This function is basically a wrapper around [pracma::findpeaks()]. If your amplicon is large, there may be fewer scans that make up individual peak. So for example you may want to set peak_scan_ramp as a smaller value.
#'
#' If too many and inappropriate peaks are being called, this may also be solved with the different repeat calling algorithms in [call_repeats()].
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#' config <- load_config()
#'
#' find_ladders(fsa_list, config)
#'
#' find_fragments(fsa_list,
#'   config,
#'   min_bp_size = 300
#' )
#'
#'
#' # Manually inspect the ladders
#' plot_traces(fsa_list,
#'   show_peaks = TRUE, n_facet_col = 1,
#'   xlim = c(400, 550), ylim = c(0, 1200)
#' )
find_fragments <- function(
    fragments_list,
    config,
    ...) {
  find_fragment_peaks <- function(trace_bp_df) {

    if(config$smoothing_window %% 2 != 1){
      stop("smoothing_window must be an odd integer value")
    }

    smoothed_signal <- pracma::savgol(
      trace_bp_df$signal,
      config$smoothing_window
    )

    # call all peaks regardless of height
    peaks <- pracma::findpeaks(smoothed_signal,
      minpeakheight = -Inf,
      peakpat = sprintf('[+]{%d,}[0]*[-]{%d,}', config$peak_scan_ramp, config$peak_scan_ramp)
    )

    n_scans <- length(trace_bp_df$signal)
    window_width <- 3

    # go through raw signal and make sure that the identified scan in the smoothed signal is still the highest
    # it will also deal with cases where the scans have the same signal (which.max will chose first)
    peak_position <- numeric(nrow(peaks))
    for (i in seq_along(peak_position)) {
      if (peaks[i, 2] + window_width > 1 & peaks[i, 2] + window_width < n_scans) { # make sure that the subsetting would be in bounds when taking window into account
        max_peak <- which.max(trace_bp_df$signal[(peaks[i, 2] - window_width):(peaks[i, 2] + window_width)])

        peak_position[i] <- peaks[i, 2] - window_width - 1 + max_peak
      } else {
        peak_position[i] <- peaks[i, 2]
      }
    }

    df <- trace_bp_df[peak_position, c("scan", "size", "signal", "off_scale")]
    colnames(df) <- c("scan", "size", "signal", "off_scale")

    # filter for minimum_peak_signal. Do it here rather than in findpeaks so that it is filtered on the raw signal value
    df <- df[which(df$signal > config$minimum_peak_signal), , drop = FALSE]

    # remove shoulder peaks
    df2 <- deshoulder(df, shoulder_window = 1.5)

    return(df2)
  }

  # prepare output file
  output <- trace_output$new("find_fragments")

  # load config
  config <- tryCatch(
    update_config(config, ...),
    error = function(e) e
  )
  if("error" %in% class(config)){
    output$set_status(
      "error", 
      config$message
    )
    return(output)
  }


  fragments_list <- lapply(fragments_list, function(x) {
    # find peak table
    df <- tryCatch(
      find_fragment_peaks(x$trace_bp_df),
      error = function(e) e
    )
    if("error" %in% class(df)){
      output$set_status(
        "error", 
        paste0("There was an error finding fragments for ", x$unique_id, ":\n", df$message)
      )
      return(output)
    }

    df$unique_id <- rep(x$unique_id, nrow(df))
    df <- df[which(df$size > config$min_bp_size & df$size < config$max_bp_size), ,drop = FALSE]
    x$peak_table_df <- df
    x$.__enclos_env__$private$min_bp_size <- config$min_bp_size
    x$.__enclos_env__$private$max_bp_size <- config$max_bp_size
    
    return(x)
  })

  return(output)
}

