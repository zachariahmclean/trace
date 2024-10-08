# read_fsa ----------------------------------------------------------------

#' Read fsa file
#'
#' This function is just a wrapper of seqinr::read.abif
#'
#' @param files a chr vector of fsa file names.
#'
#' @return A list of fragments_trace objects
#' @seealso [seqinr::read.abif()]
#' @export
#' @importFrom seqinr read.abif
#'
#' @examples
#' # files <- list.files("example_directory/", full.names = TRUE, pattern = ".fsa")
#' # read_fsa(files)
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
