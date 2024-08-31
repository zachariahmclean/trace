
# metrics_override_helper ---------------------------------------------------------





#' Assign index peaks
#'
#' Assign index peaks in preparation for calculation of instability metrics
#'
#' @param fragments_list A list of "fragments_repeats" class objects representing
#' fragment data.
#' @param grouped Logical value indicating whether samples should be grouped to
#' share a common index peak. `FALSE` will assign the sample's own modal allele as the index peak. `TRUE` will use metadata to assign the index peak based on the modal peak of another sample. This is useful for cases like inferring repeat size of inherited alleles from mouse tail data. Requires metadata via \code{link{add_metadata()}}.
#' @param index_override_dataframe A data.frame to manually set index peaks.
#' Column 1: unique sample IDs, Column 2: desired index peaks (the order of the
#' columns is important since the information is pulled by column position rather
#' than column name). Closest peak in each sample is selected.
#'
#' @return A list of \code{"fragments_repeats"} objects with index_repeat and index_height added.
#' @details
#' A key part of several instability metrics is the index peak. This is the repeat
#' length used as the reference peak for relative instability metrics calculations, like expansion index or average repeat gain.
#' For example, this is the the inherited repeat length of a mouse, or the modal repeat length for the cell line at a starting time point.
#'
#'
#' If `grouped` is set to `TRUE`, this function groups the samples by their metrics_group_id and uses the samples set as metrics_baseline_control to set the index peak. Use \code{link{add_metadata()}} to set these variables. For mice, if just a few samples have the inherited repeat height shorter than the expanded population, you could not worry about this and instead use the `index_override_dataframe`. This can be used to manually override these assigned index repeat values (irrespective of whether `grouped` is TRUE or FALSE).
#'
#' As a final option, the index peak could be manually assigned directly to a \code{link{fragments_repeats}} using the internal setter function \code{link{fragments_repeats$set_index_peak()}}.
#'
#' @export
#'
#' @examples
#'
#'
#' file_list <- trace::cell_line_fsa_list
#'
#' ladder_list <- find_ladders(file_list)
#'
#' fragments_list <- find_fragments(ladder_list,
#'   min_bp_size = 300
#' )
#'
#' allele_list <- find_alleles(
#'   fragments_list = fragments_list
#' )
#' repeats_list <- call_repeats(
#'   fragments_list = allele_list
#' )
#'
#' metadata_added_list <- add_metadata(
#'   fragments_list = repeats_list,
#'   metadata_data.frame = trace::metadata
#' )
#'
#'index_assigned <- assign_index_peaks(metadata_added_list,
#'                                     grouped = TRUE)
#'
#' plot_traces(index_assigned[1], xlim = c(100,150))
#'
#'
#'
#'
#'
#'
assign_index_peaks <- function(
    fragments_list,
    grouped = FALSE,
    index_override_dataframe = NULL) {
  if (grouped == TRUE) {
    # what we're doing here is pulling out the key data for all the samples that are metrics controls
    # each sample will then have the data for their appropriate control inserted inside
    # that can then be used in the calculation of instability metrics

    # we need to insert the whole peak table of the control because the calculation of the weighted mean
    # looks at a specific subset of the table, which is not set until the calculate metrics function

    # make a list of dataframes and alleles for each of the controls for the groups
    metrics_group_ids <- sapply(fragments_list, function(x) x$metrics_group_id)
    unique_metrics_group_ids <- unique(metrics_group_ids)

    baseline_control_list <- vector("list", length(unique_metrics_group_ids))
    names(baseline_control_list) <- unique_metrics_group_ids

    for (i in seq_along(fragments_list)) {
      if (fragments_list[[i]]$metrics_baseline_control == TRUE) {
        # since there can be more than one control, make a list of them
        baseline_control_list[[fragments_list[[i]]$metrics_group_id]] <- c(
          baseline_control_list[[fragments_list[[i]]$metrics_group_id]],
          list(
            list(
              fragments_list[[i]]$get_alleles()$allele_1_repeat,
              fragments_list[[i]]$repeat_table_df
            )
          )
        )
      }
    }

    # do some quality control

    for (i in seq_along(baseline_control_list)) {
      controls_missing_allele <- all(sapply(baseline_control_list[[fragments_list[[i]]$metrics_group_id]], function(x) is.na(x[[1]])))

      if (length(baseline_control_list[[i]]) == 0) {
        stop(paste0("Group '", names(baseline_control_list)[[i]], "' has no 'metrics_baseline_control'. Go back to metadata to check that each group has a baseline control, or remove samples from the list for analysis with 'remove_fragments()' if it doesn't make sense to include them beyond this point (eg size standards or no template controls)"),
          call. = FALSE
        )
      } else if (controls_missing_allele == TRUE) {
        stop(paste0("Group '", names(baseline_control_list)[[i]], "' control has no allele called. Grouped analysis won't work for these samples."),
          call. = FALSE
        )
      } else if (length(baseline_control_list[[i]]) > 1) {
        message(paste0("Group '", names(baseline_control_list)[[i]], "' has more than one 'metrics_baseline_control'. The median repeat of the assigned samples will be used to assign the index peak"))
      }
    }

    # loop over each sample and put data inside
    for (i in seq_along(fragments_list)) {
      fragments_list[[i]]$.__enclos_env__$private$index_samples <- baseline_control_list[[fragments_list[[i]]$metrics_group_id]]

      control_index_median_repeat <- median(sapply(baseline_control_list[[fragments_list[[i]]$metrics_group_id]], function(x) x[[1]]))


      # since the repeat size may not be an integer, need to find what the closest peak is to the control sample
      # delta between repeat of index sample and all repeats of sample.

      # skip samples with no data

      if (nrow(fragments_list[[i]]$repeat_table_df) > 0) {
        index_delta <- fragments_list[[i]]$repeat_table_df$repeats - control_index_median_repeat
        closest_peak <- which(abs(index_delta) == min(abs(index_delta)))

        if (length(closest_peak) == 1) {
          fragments_list[[i]]$set_index_peak(fragments_list[[i]]$repeat_table_df$repeats[closest_peak])
        } else {
          tallest_candidate <- closest_peak[which(fragments_list[[i]]$repeat_table_df$height[closest_peak] == max(fragments_list[[i]]$repeat_table_df$height[closest_peak]))]
          fragments_list[[i]]$set_index_peak(fragments_list[[i]]$repeat_table_df$repeats[tallest_candidate])
        }
      } else {
        fragments_list[[i]]$set_index_peak(NA_real_)
      }
    }
  } else {
    # otherwise just use the modal peak as the index peak
    fragments_list <- lapply(fragments_list, function(x) {
      x$set_index_peak(x$get_alleles()$allele_1_repeat)
      return(x)
    })
  }

  # override index peak with manually supplied values
  if (!is.null(index_override_dataframe)) {
    index_override_dataframe <- as.data.frame(index_override_dataframe)

    if (!any(index_override_dataframe[, 1] %in% names(fragments_list))) {
      missing_unique_ids <- which(!index_override_dataframe[, 1] %in% names(fragments_list))

      warning(
        call. = FALSE,
        paste0(
          "The following unique ids from the index override data frame are not in the repeats list:",
          paste0(index_override_dataframe[, 1], collapse = ", ")
        )
      )
    }

    lapply(fragments_list, function(x) {
      # if there is nothing to override, then just return the existing index values
      if (any(index_override_dataframe[, 1] == x$unique_id)) {
        index_delta <- x$repeat_table_df$repeats - index_override_dataframe[which(index_override_dataframe[, 1] == x$unique_id), 2]

        closest_peak <- which(abs(index_delta) == min(abs(index_delta)))
        if (length(closest_peak) == 1) {
          x$set_index_peak(x$repeat_table_df$repeats[closest_peak])
        } else {
          tallest_candidate <- closest_peak[which(x$repeat_table_df$height[closest_peak] == max(x$repeat_table_df$height[closest_peak]))]
          x$set_index_peak(x$repeat_table_df$repeats[tallest_candidate])
        }
      }
      return(x)
    })
  }

  return(fragments_list)
}

