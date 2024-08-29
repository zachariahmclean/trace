# find fragments ------------------------------------------------------------




#' Find fragment peaks
#'
#' Find fragment peaks in continuous trace data and convert to fragments_repeats
#' class.
#'
#' @param fragments_trace_list A list of fragments_trace objects containing fragment data.
#' @param smoothing_window numeric: signal smoothing window size passed to pracma::savgol()
#' @param minimum_peak_signal numeric: minimum height of peak from smoothed trace
#' @param min_bp_size numeric: minimum bp size of peaks to consider
#' @param max_bp_size numeric: maximum bp size of peaks to consider
#' @param ... pass additional arguments to findpeaks, or change the default arguments
#' we set. minimum_peak_signal above is passed to findpeaks as minpeakheight, and
#' peakpat has been set to '\[+\]\{6,\}\[0\]*\[-\]\{6,\}' so that peaks with flat tops are
#' still called, see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
#'
#'
#' @return a list of fragments_repeats objects, equal length and names to the list input
#' @export
#'
#' @importFrom pracma findpeaks
#' @importFrom pracma savgol
#'
#' @details
#' This function is basically a wrapper around pracma::findpeaks. As mentioned above,
#' the default arguments arguments of pracma::findpeaks can be changed by passing them
#' to find_fragments with ... .
#'
#' If too many and inappropriate peaks are being called, this may also be solved with the different repeat calling algorithms in [call_repeats()].
#'
#' @examples
#' file_list <- trace::cell_line_fsa_list
#'
#' test_ladders <- find_ladders(file_list[1])
#'
#' fragments_list <- find_fragments(test_ladders,
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
    # it will also deal with cases where the scans have the same height (which.max will chose first)
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
    colnames(df) <- c("scan", "size", "height", "off_scale")

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
}




## set class from peak table ----------------------------------------------------------------

#' Convert Peak Table to Fragments_repeats class
#'
#' This function converts a peak table data frame into a list of fragments_repeats objects.
#'
#' @param df A data frame containing the peak data.
#' @param data_format The format that the data frame is in (for example, a genemapper peak table). Choose between: genemapper5, generic.
#' @param unique_id A character string specifying column name giving the unique sample id (often the file name).
#' @param peak_size_col A character string specifying column name giving the peak size.
#' @param peak_height_col A character string specifying column name giving the peak height.
#' @param min_size_bp Numeric value indicating the minimum size of the peak table to import.
#' @param max_size_bp Numeric value indicating the maximum size of the peak table to import.
#' @param dye_col Genemapper specific. A character string specifying column name indicating the dye channel.
#' @param dye_channel Genemapper specific. A character string indicating the channel to extract data from. For example, 6-FAM is often "B".
#' @param allele_col Genemapper specific. A character string specifying column name indicating the called alleles. This is often used when the peaks have been called in genemapper.
#'
#' @return A list of fragments_repeats. objects.
#'
#' @details This function takes a peak table data frame (eg. Genemapper output) and converts it into a list of fragment objects.
#' The function supports different data formats and allows specifying column names for various attributes.
#'
#' @seealso \code{\link{repeat_table_to_repeats}}
#'
#' @examples
#'
#' gm_raw <- trace::example_data
#'
#' test_fragments <- peak_table_to_fragments(
#'   gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B",
#'   min_size_bp = 400
#' )
#'
#' @export
peak_table_to_fragments <- function(df,
                                    data_format = NULL,
                                    peak_size_col = NULL,
                                    peak_height_col = NULL,
                                    unique_id = NULL,
                                    dye_col = NULL,
                                    dye_channel = NULL,
                                    allele_col = NULL,
                                    min_size_bp = 100,
                                    max_size_bp = 1000) {
  # check to make sure that if the user supplies a column name, that it's actually in the dataframe
  if (any(!is.null(peak_size_col), !is.null(peak_height_col), !is.null(unique_id))) {
    function_input_vector <- c(peak_size_col, peak_height_col, unique_id)
    function_input_name_vector <- c("peak_size_col", "peak_height_col", "unique_id")
    for (i in seq_along(function_input_vector)) {
      if (!any(names(df) == function_input_vector[[i]])) {
        stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ", function_input_name_vector[[i]], " input"),
          call. = FALSE
        )
      }
    }
  }

  # chose the tidying function
  # Use the supplied user column names if given
  if (data_format == "genemapper5") {
    df2 <- clean_genemapper5(df,
      peak_size_col = ifelse(length(peak_size_col) == 0, "Size", peak_size_col),
      peak_height_col = ifelse(length(peak_height_col) == 0, "Height", peak_height_col),
      unique_id = ifelse(length(unique_id) == 0, "Sample.File.Name", unique_id),
      dye_col = ifelse(length(dye_col) == 0, "Dye.Sample.Peak", dye_col),
      dye_channel = ifelse(length(dye_channel) == 0, "B", dye_channel),
      allele_col = ifelse(length(allele_col) == 0, "Allele", allele_col)
    )
  } else if (data_format == "generic") {
    df2 <- clean_generic(df,
      peak_size_col = peak_size_col,
      peak_height_col = peak_height_col,
      unique_id = unique_id
    )
  } else {
    stop("Data format not recognised. Choose between: genemapper5, generic",
      call. = FALSE
    )
  }

  # filter size and split up into a list of fragments
  fragments_list <-
    lapply(
      split(df2, df2$unique_id),
      function(x) {
        # filter size
        df <- x[x$size > min_size_bp & x$size < max_size_bp & !is.na(x$size), , drop = FALSE]
        # check to see if all rows removed and give warning
        if (nrow(df) == 0) {
          warning(paste0("Size filtering removed all rows for ", unique(x$unique_id)),
            call. = FALSE
          )
        }

        new_fragments_repeats <- fragments_repeats$new(unique_id = unique(x$unique_id))
        new_fragments_repeats$peak_table_df <- df

        return(new_fragments_repeats)
      }
    )

  return(fragments_list)
}



## set class from repeats table ----------------------------------------------------------------

#' Convert Repeat Table to Repeats Fragments
#'
#' This function converts a repeat table data frame into a list of fragments_repeats. class.
#'
#' @param df A data frame containing the repeat data.
#' @param unique_id A character string indicating the column name for unique identifiers.
#' @param repeat_col A character string indicating the column name for the repeats.
#' @param frequency_col A character string indicating the column name for the repeat frequencies.
#'
#' @return A list of fragments_repeats.
#'
#' @details This function takes a repeat table data frame and converts it into a list of repeats fragments.
#' The function allows specifying column names for repeats, frequencies, and unique identifiers.
#' @export
#'
#' @examples
#' repeat_table <- trace::example_data_repeat_table
#' test_fragments <- repeat_table_to_repeats(
#'   repeat_table,
#'   repeat_col = "repeats",
#'   frequency_col = "height",
#'   unique_id = "unique_id"
#' )
repeat_table_to_repeats <- function(df,
                                    unique_id,
                                    repeat_col,
                                    frequency_col) {
  # validate inputs to give good errors to user
  ## check to make sure that if the user supplies a column name, that it's actually in the dataframe
  function_input_vector <- c(repeat_col, frequency_col, unique_id)
  function_input_name_vector <- c("repeat_col", "frequency_col", "unique_id")
  for (i in seq_along(function_input_vector)) {
    if (!any(names(df) == function_input_vector[[i]])) {
      stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ", function_input_name_vector[[i]], " input"),
        call. = FALSE
      )
    }
  }

  names(df)[names(df) == repeat_col] <- "repeats"
  names(df)[names(df) == frequency_col] <- "height"
  names(df)[names(df) == unique_id] <- "unique_id"

  repeats_list <- lapply(
    split(df, df$unique_id),
    function(x) {
      new_fragments_repeats <- fragments_repeats$new(unique_id = unique(x$unique_id))
      new_fragments_repeats$repeat_table_df <- x
      return(new_fragments_repeats)
    }
  )


  return(repeats_list)
}

