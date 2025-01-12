# read_fsa ----------------------------------------------------------------

#' Read fsa file
#'
#' Read fsa file into memory and create fragments_trace object
#'
#' @param files a chr vector of fsa file names. For example, return all the fsa files in a directory with 'list.files("example_directory/", full.names = TRUE, pattern = ".fsa")'.
#' 
#' @details
#' read_fsa is just a wrapper around [seqinr::read.abif()] that reads the fsa file into memory and stores it inside a fragments_trace object. That enables you to use the next function [find_ladders()].
#'
#' @return A list of fragments_trace objects
#' @seealso [find_ladders()], [plot_data_channels()]
#' @export
#' @importFrom seqinr read.abif
#'
#' @examples
#' 
#' fsa_file <- read_fsa(system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr"))
#' plot_data_channels(fsa_file)
#'
read_fsa <- function(
  files) {
  # make sure file extension is fsa
  unique_file_ext <- unique(tools::file_ext(files))
  if (length(unique_file_ext) > 1) {
    stop("Files must be only be .fsa")
  }
  if (unique_file_ext != "fsa") {
    stop("Files must be .fsa")
  }

  # read in samples
  unique_names <- make.unique(basename(files))
  fragments_list <- vector("list", length(files))
  names(fragments_list) <- unique_names

  for (i in seq_along(files)) {
    fragments_list[[i]] <- fragments_trace$new(
      unique_id = names(fragments_list[i]),
      fsa_file = seqinr::read.abif(files[i])
    )
  }
  return(fragments_list)
}
