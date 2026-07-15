# read_fsa ----------------------------------------------------------------

#' Read fsa file
#'
#' Read fsa file into memory and create fragments object
#'
#' @param files a chr vector of fsa file names. For example, return all the fsa files in a directory with 'list.files("example_directory/", full.names = TRUE, pattern = ".fsa")'.
#' @param unique_id_from_fsa Logical. If `TRUE` (default), the `unique_id` of each sample is built from the fsa file itself as the sample name (the `SpNm` tag) joined to the run id (the `RunN` tag) with an underscore, e.g. "20230413_A07_Run_MRSPHILLIPS2_2023-04-14_11-30_0063". This guarantees ids are unique across fragment analysis runs, even when the same file name is reused in different batches. If `FALSE`, the file name is used (the previous behaviour). When a sample name or run id is missing from the fsa file the file name is used as a fallback, and a numeric suffix is added to any ids that are still not unique.
#'
#' @details
#' read_fsa is a wrapper around [seqinr::read.abif()] that reads the fsa file into memory and stores it inside a fragments object. That enables you to use the next function [find_ladders()]. It also parses useful header fields from the fsa file (see [extract_fsa_metadata()]) and, by default, uses them to build a `unique_id` that is unique across runs.
#'
#' @return A list of fragments objects
#' @seealso [find_ladders()], [plot_data_channels()], [extract_fsa_metadata()]
#' @export
#' @importFrom seqinr read.abif
#'
#' @examples
#'
#' fsa_file <- read_fsa(system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr"))
#' plot_data_channels(fsa_file)
#'
read_fsa <- function(
  files,
  unique_id_from_fsa = TRUE) {
  # make sure file extension is fsa
  unique_file_ext <- unique(tools::file_ext(files))
  if (length(unique_file_ext) > 1) {
    stop("Files must be only be .fsa")
  }

  # read each file once and parse its header metadata
  abif_list <- vector("list", length(files))
  metadata_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    abif_list[[i]] <- seqinr::read.abif(files[i])
    metadata_list[[i]] <- parse_fsa_metadata(abif_list[[i]])
  }

  # build unique ids from the fsa (sample name + run id), falling back to the
  # file name when those tags are missing
  if (unique_id_from_fsa) {
    unique_names <- vapply(seq_along(files), function(i) {
      sample_name <- metadata_list[[i]]$sample_name
      run_id <- metadata_list[[i]]$run_id
      if (!is.na(sample_name) && !is.na(run_id)) {
        paste(sample_name, run_id, sep = "_")
      } else if (!is.na(sample_name)) {
        sample_name
      } else {
        basename(files[i])
      }
    }, character(1))
  } else {
    unique_names <- basename(files)
  }
  # final safety net so ids are always unique (e.g. a sample name repeated within a run)
  unique_names <- make.unique(unique_names)

  fragments_list <- vector("list", length(files))
  names(fragments_list) <- unique_names

  for (i in seq_along(files)) {
    fragments_list[[i]] <- fragments$new(unique_names[i], "fsa", abif_list[[i]])
    fragments_list[[i]]$fsa_metadata <- metadata_list[[i]]
  }
  return(fragments_list)
}

 
 ## set class from size table ----------------------------------------------------------------

#' Convert Size Table to Fragments
#'
#' This function converts a size table data frame into a list of fragments. class.
#'
#' @param df A data frame containing the size data with the columns "unique_id", "size", "signal".
#' @param min_size_bp Numeric value indicating the minimum size of the peak table to import.
#' @param max_size_bp Numeric value indicating the maximum size of the peak table to import.
#'
#' @return A list of fragments objects.
#'
#' @details This function takes a size table data frame and converts it into a list of fragments objects. The column names must be "unique_id" (unique sample id), "size" (base pair size), and "signal" (the signal associated with fragment). The dataframe should be long (eg bind rows if the data are separate).
#' @export
#' @seealso [repeat_table_to_fragments()], [size_table_to_fragments()], [read_fsa()]
#'
#' @examples
#' size_table <- trace::example_data
#' colnames(size_table)[c(2, 5, 6)] <- c("unique_id", "size", "signal")
#' fragments_list <- size_table_to_fragments(size_table)
size_table_to_fragments <- function(
    df,
    min_size_bp = 200,
    max_size_bp = 1000) {
  # need to make sure table is dataframe (an not a tibble)
  if(class(df)[1] == 'tbl_df'){
    df <- as.data.frame(df)
    }
  
  # make sure df has the required columns
  if(!all(c("unique_id", "size", "signal") %in% names(df))){
    stop(call. = FALSE, "Dataframe must have the column names 'unique_id', 'size', and 'signal'")
    }

  # filter size and split up into a list of fragments
  fragments_list <-lapply(
    split(df, df$unique_id),
    function(x) {
      # filter size
      df <- x[x$size > min_size_bp & x$size < max_size_bp & !is.na(x$size), , drop = FALSE]
      # check to see if all rows removed and give warning
      if (nrow(df) == 0) {
        warning(paste0("Size filtering removed all rows for ", unique(x$unique_id)),
          call. = FALSE
        )
      }

      new_fragments <- fragments$new(unique(x$unique_id), "size", df)

      return(new_fragments)
    }
  )

  return(fragments_list)
}

## set class from genemapper peak table ----------------------------------------------------------------

#' Convert Genemapper Peak Table to fragments class
#'
#' This function converts a genemapper peak table data frame into a list of fragments objects.
#'
#' @param df A data frame from genemapper containing the peak data. 
#' @param min_size_bp Numeric value indicating the minimum size of the peak table to import.
#' @param max_size_bp Numeric value indicating the maximum size of the peak table to import.
#' @param dye_channel A character string indicating the Genemapper channel to extract data from. This is used to filter 'Dye.Sample.Peak' for the channel containing the data. For example, 6-FAM is often "B" while ladder is "O".
#'
#' @return A list of fragments objects.
#'
#' @details This function takes a peak table data frame (eg. Genemapper 5 output) and converts it into a list of fragment objects. It uses the "Sample.File.Name" as the unique id, so make sure that each file as a unique name. Column names should contain: "Dye.Sample.Peak", "Sample.File.Name", "Allele", "Size", "Height".
#'
#' @seealso [repeat_table_to_fragments()], [size_table_to_fragments()], [read_fsa()]
#' @export
#'
#' @examples
#'
#' gm_raw <- trace::example_data
#'
#' test_fragments <- genemapper_table_to_fragments(
#'   gm_raw,
#'   dye_channel = "B",
#'   min_size_bp = 400
#' )
#'
genemapper_table_to_fragments <- function(
  df,
  dye_channel,
  min_size_bp = 200,
  max_size_bp = 1000) {
  
  # make sure df has the required columns
  if(!all(c("Dye.Sample.Peak", "Sample.File.Name", "Allele", "Size", "Height") %in% names(df))){
    stop(call. = FALSE, "Dataframe must have the column names 'Dye.Sample.Peak', 'Sample.File.Name', 'Allele', 'Size', 'Height'")
    }
  
  # extract channel
  df$dye <- substr(df$Dye.Sample.Peak, 1, 1)

  # filter for only the dye of interest
  df <- df[which(df$dye == dye_channel), , drop = FALSE]
  names(df)[names(df) == "Size"] <- "size"
  names(df)[names(df) == "Height"] <- "signal"
  names(df)[names(df) == "Sample.File.Name"] <- "unique_id"

  fragments_list <- size_table_to_fragments(
    df,
    min_size_bp = min_size_bp,
    max_size_bp = max_size_bp
  )

  return(fragments_list)
}


## set class from repeats table ----------------------------------------------------------------

#' Convert Repeat Table to Repeats Fragments
#'
#' This function converts a repeat table data frame into a list of fragments. class.
#'
#' @param df A data frame containing the repeat data with the columns "unique_id", "repeats", "signal".
#' @param min_repeat minimum repeat size
#' @param max_repeat maximum repeat size
#'
#' @return A list of fragments objects.
#'
#' @details This function takes a repeat table data frame and converts it into a list of fragments objects. The column names must be "unique_id" (unique sample id), "repeats" (number of repeats), and "signal" (the signal associated with the number of repeats. for example a frequency count of reads.). The dataframe should be long (eg bind rows if the data are separate).
#' @export
#'
#' @examples
#' repeat_table <- trace::example_data_repeat_table
#' test_fragments <- repeat_table_to_fragments(repeat_table)
repeat_table_to_fragments <- function(df,
  min_repeat = 0,
  max_repeat = 1000) {
  # need to make sure table is dataframe (an not a tibble)
  if(class(df)[1] == 'tbl_df'){
    df <- as.data.frame(df)
    }
  
  # make sure df has the required columns
  if(!all(c("unique_id", "repeats",   "signal") %in% names(df))){
    stop(call. = FALSE, "Dataframe must have the column names 'unique_id', 'repeats', and 'signal'")
    }

  repeats_list <- lapply(
    split(df, df$unique_id),
    function(x) {
      # filter size
      x <- x[x$repeats > min_repeat & x$repeats < max_repeat & !is.na(x$repeats), , drop = FALSE]
      new_fragments <- fragments$new(unique(x$unique_id), "repeats", x)
      return(new_fragments)
    }
  )

  return(repeats_list)
}
