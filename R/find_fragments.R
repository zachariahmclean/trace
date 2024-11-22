# find fragments ------------------------------------------------------------




#' Find fragment peaks
#'
#' Find fragment peaks in continuous trace data and convert to fragments_repeats
#' class.
#'
#' @param fragments_trace_list A list of fragments_trace objects containing fragment data.
#' @param smoothing_window numeric: signal smoothing window size passed to pracma::savgol()
#' @param minimum_peak_signal numeric: minimum signal of peak from smoothed trace
#' @param min_bp_size numeric: minimum bp size of peaks to consider
#' @param max_bp_size numeric: maximum bp size of peaks to consider
#' @param ... pass additional arguments to findpeaks, or change the default arguments
#' we set. minimum_peak_signal above is passed to findpeaks as minpeakheight, and
#' peakpat has been set to '\[+\]\{6,\}\[0\]*\[-\]\{6,\}' so that peaks with flat tops are
#' still called, see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
#'
#'
#' @return a list of fragments_repeats objects.
#' @export
#'
#' @importFrom pracma findpeaks
#' @importFrom pracma savgol
#'
#' @details
#' 
#' [find_fragments()] takes in a list of fragments_trace objects and returns a list of new fragments_repeats objects.
#' 
#' This function is basically a wrapper around pracma::findpeaks. As mentioned above,
#' the default arguments arguments of pracma::findpeaks can be changed by passing them
#' to find_fragments with ... .
#'
#' If too many and inappropriate peaks are being called, this may also be solved with the different repeat calling algorithms in [call_repeats()].
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list)
#'
#' fragments_list <- find_fragments(fsa_list,
#'   min_bp_size = 300
#' )
#'
#'
#' # Manually inspect the ladders
#' plot_traces(fragments_list,
#'   show_peaks = TRUE, n_facet_col = 1,
#'   xlim = c(400, 550), ylim = c(0, 1200)
#' )
find_fragments <- function(
    fragments_trace_list,
    smoothing_window = 21,
    minimum_peak_signal = 20,
    min_bp_size = 100,
    max_bp_size = 1000,
    ...) {
  find_fragment_peaks <- function(trace_bp_df,
                                  ...) {
    smoothed_signal <- pracma::savgol(
      trace_bp_df$signal,
      smoothing_window
    )

    # deals with cases of user overriding values
    if ("peakpat" %in% ...names()) {
      peaks <- pracma::findpeaks(smoothed_signal,
        minpeakheight = minimum_peak_signal,
        ...
      )
    } else if ("minpeakheight" %in% ...names()) {
      # user minpeakheight instead of minimum_peak_signal
      peaks <- pracma::findpeaks(smoothed_signal,
        peakpat = "[+]{6,}[0]*[-]{6,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
        ...
      )
    } else {
      peaks <- pracma::findpeaks(smoothed_signal,
        peakpat = "[+]{6,}[0]*[-]{6,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
        minpeakheight = minimum_peak_signal,
        ...
      )
    }

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

    # remove shoulder peaks
    df2 <- deshoulder(df, shoulder_window = 1.5)

    return(df2)
  }

  fragments_list <- lapply(fragments_trace_list, function(x) {
    # find peak table
    df <- find_fragment_peaks(x$trace_bp_df, ...)
    df$unique_id <- rep(x$unique_id, nrow(df))
    df <- df[which(df$size > min_bp_size & df$size < max_bp_size), ]

    # generate new class
    new_fragments_repeats <- fragments_repeats$new(unique_id = x$unique_id)
    new_fragments_repeats$trace_bp_df <- x$trace_bp_df
    new_fragments_repeats$peak_table_df <- df
    new_fragments_repeats <- transfer_metadata_helper(x, new_fragments_repeats)
    new_fragments_repeats$.__enclos_env__$private$min_bp_size <- min_bp_size
    new_fragments_repeats$.__enclos_env__$private$max_bp_size <- max_bp_size

    return(new_fragments_repeats)
  })

  return(fragments_list)
}

