# genemapper5_tidier ------------------------------------------------------

clean_genemapper5 <- function(
    df,
    peak_size_col,
    peak_height_col,
    unique_id,
    dye_col,
    dye_channel,
    allele_col) {
  df <- as.data.frame(df)
  # rename cols based on used supplied name
  names(df)[names(df) == peak_size_col] <- "size"
  names(df)[names(df) == peak_height_col] <- "height"
  names(df)[names(df) == unique_id] <- "unique_id"
  names(df)[names(df) == dye_col] <- "Dye.Sample.Peak"
  names(df)[names(df) == allele_col] <- "allele"

  # extract channel
  df$dye <- substr(df$Dye.Sample.Peak, 1, 1)

  # tidy dataframe names to snake case
  names(df) <- tolower(names(df))
  names(df) <- gsub("\\.", "_", names(df))

  # filter for only the dye of interest
  df2 <- df[which(df$dye == dye_channel), , drop = FALSE]

  return(df2)
}

# for generic peak table

clean_generic <- function(
    df,
    peak_size_col,
    peak_height_col,
    unique_id) {
  df <- as.data.frame(df)
  # rename cols based on used supplied name
  names(df)[names(df) == peak_size_col] <- "size"
  names(df)[names(df) == peak_height_col] <- "height"
  names(df)[names(df) == unique_id] <- "unique_id"

  # tidy dataframe names to snake case
  names(df) <- tolower(names(df))
  names(df) <- gsub("\\.", "_", names(df))

  return(df)
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
peak_table_to_fragments <- function(
    df,
    data_format = NULL,
    peak_size_col = NULL,
    peak_height_col = NULL,
    unique_id = NULL,
    dye_col = NULL,
    dye_channel = NULL,
    allele_col = NULL,
    min_size_bp = 200,
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
    stop("Data format not recognized. Choose between: genemapper5, generic",
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
repeat_table_to_repeats <- function(
    df,
    unique_id,
    repeat_col,
    frequency_col) {
  # need to make sure table is dataframe (an not a tibble)
  df <- as.data.frame(df)

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
