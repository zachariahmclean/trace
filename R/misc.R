print_helper <- function(fragment,
                         exclude = NULL,
                         sample_attrs) {
  unique_id <- fragment$unique_id

  cat(paste0("\033[1;34m< ", class(fragment)[1], " object >\033[0m\n"))
  cat("\033[1;36m-----------------------------\033[0m\n")

  # Section: Sample Attributes
  for (attr in sample_attrs) {
    if (attr %in% names(fragment)) {
      value <- fragment[[attr]]
      cat(sprintf("\033[1;33m%-30s\033[0m", attr))

      if (is.null(value)) {
        cat("NULL\n")
      } else if (is.logical(value)) {
        cat(ifelse(value, "\033[1;32mTRUE\033[0m", "\033[1;31mFALSE\033[0m"), "\n")
      } else if (is.numeric(value)) {
        if (length(value) == 1) {
          cat(format(value), "\n")
        } else {
          cat(format(paste("numeric vector length", length(value))), "\n")
        }
      } else if (is.character(value)) {
        if (length(value) == 1) {
          cat(format(value), "\n")
        } else {
          cat(format(paste("character vector length", length(value))), "\n")
        }
      } else {
        cat(class(value), "\n")
      }
    }
  }

  cat("\033[1;36m-----------------------------\033[0m\n")

  # Check and print called alleles in private fields of fragments
  private_names <- ls(fragment$.__enclos_env__$private, all.names = TRUE)
  alleles_names <- private_names[which(grepl("allele_", private_names))]
  #should allele_2 be included? remove if is na since it's only an edge case
  if(is.na(fragment$.__enclos_env__$private$allele_2_signal)){
    alleles_names <- alleles_names[-grep("allele_2_", alleles_names)]
  }
  alleles_names <- c(alleles_names, "index_repeat")
  for (name in alleles_names) {
    value <- fragment$.__enclos_env__$private[[name]]
    class_value <- class(value)

    cat(sprintf("\033[1m%-30s\033[0m", name))

    if (is.null(value)) {
      cat("NULL\n")
    } else if (is.numeric(value)) {
      cat(format(value), "\n")
    }
  }
  

  slot_names <- ls(fragment, all.names = TRUE)
  all_exclusions <- c(
    exclude,
    sample_attrs,
    "input_method", "raw_data", "raw_ladder", "scan", "off_scale_scans",
    slot_names[which(sapply(slot_names, function(name) class(fragment[[name]]) == "function"))],
    ".__active__",
    ".__enclos_env__"
  )

  if(fragment$input_method != "fsa"){
    all_exclusions <- c(all_exclusions, "fsa", "ladder_df", "trace_bp_df")
  }

  slot_names <- setdiff(slot_names, all_exclusions)

  for (name in slot_names) {
    value <- fragment[[name]]
    class_value <- class(value)

    cat(sprintf("\033[1m%-30s\033[0m", name))

    if (is.null(value)) {
      cat("NULL\n")
    } else if (is.logical(value)) {
      cat(ifelse(value, "\033[1;32mTRUE\033[0m", "\033[1;31mFALSE\033[0m"), "\n")
    } else if (is.numeric(value)) {
      if (length(value) == 1) {
        cat(format(value), "\n")
      } else {
        cat(format(paste("numeric vector length", length(value))), "\n")
      }
    } else if (is.character(value)) {
      if (length(value) == 1) {
        cat(format(value), "\n")
      } else {
        cat(format(paste("character vector length", length(value))), "\n")
      }
    } else if (is.data.frame(value)) {
      cat(sprintf("data.frame: %d rows x %d cols\n", nrow(value), ncol(value)))
    } else {
      cat(class_value, "\n")
    }
  }

  invisible(fragment)
}


# remove fragments -------------------------------------------------------

#' Remove Samples from List
#'
#' A convenient function to remove specific samples from a list of fragments.
#'
#' @param fragments_list A list of fragments objects containing fragment data.
#' @param samples_to_remove A character vector containing the unique IDs of the samples to be removed.
#'
#' @return A modified list of fragments with the specified samples removed.
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#'
#' test_fragments <- genemapper_table_to_fragments(
#'   gm_raw,
#'   dye_channel = "B",
#'   min_size_bp = 300
#' )
#'
#' all_fragment_names <- names(test_fragments)
#'
#' # pull out unique ids of samples to remove
#' samples_to_remove <- all_fragment_names[c(1, 5, 10)]
#'
#' samples_removed <- remove_fragments(test_fragments, samples_to_remove)
#'
remove_fragments <- function(
    fragments_list,
    samples_to_remove) {
  unique_ids <- vector("numeric", length = length(fragments_list))
  for (i in seq_along(fragments_list)) {
    unique_ids[[i]] <- fragments_list[[i]]$unique_id
  }
  samples_removed <- fragments_list
  suppressWarnings(
    samples_removed[which(unique_ids %in% samples_to_remove)] <- NULL
  )
  return(samples_removed)
}