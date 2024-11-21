
# metrics_override_helper ---------------------------------------------------------





#' Assign index peaks
#'
#' Assign index peaks in preparation for calculation of instability metrics
#'
#' @param fragments_list A list of "fragments_repeats" class objects representing
#' fragment data.
#' @param grouped Logical value indicating whether samples should be grouped to
#' share a common index peak. `FALSE` will assign the sample's own modal allele as the index peak. `TRUE` will use metadata to assign the index peak based on the modal peak of another sample (see below for more details).
#' @param index_override_dataframe A data.frame to manually set index peaks.
#' Column 1: unique sample IDs, Column 2: desired index peaks (the order of the
#' columns is important since the information is pulled by column position rather
#' than column name). Closest peak in each sample is selected so the number needs to just be approximate.
#'
#' @return This function modifies list of fragments_repeats objects in place with index_repeat and index_signal added.
#' @details
#' A key part of instability metrics is the index peak. This is the repeat
#' length used as the reference peak for relative instability metrics calculations, like expansion index.
#' This is usually the the inherited repeat length of a mouse, or the modal repeat length for the cell line at a starting time point.
#'
#' If `grouped` is set to `TRUE`, this function groups the samples by their `metrics_group_id` and uses the samples set as `metrics_baseline_control` to set the index peak. Use [add_metadata()] to set these variables. This is useful for cases like inferring repeat size of inherited alleles from mouse tail data. If the samples that are going to be used to assign index peak are from different fragment analysis runs, use `correction = "batch"` in [call_repeats()] to make sure the systematic differences between runs are corrected and the correct index peak is assigned. If there are multiple samples used as baseline control, the median value will be used to assign index peak to corresponding samples.
#' 
#' For mice, if just a few samples have the inherited repeat signal shorter than the expanded population, you could not worry about this and instead use the `index_override_dataframe`. This can be used to manually override these assigned index repeat values (irrespective of whether `grouped` is TRUE or FALSE).
#'
#' As a final option, the index peak could be manually assigned directly to a [fragments_repeats] class using the internal setter function fragments_repeats$set_index_peak().
#'
#' @export
#'
#' @examples
#'
#'
#' fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' fragments_list <- find_fragments(fsa_list,
#'   min_bp_size = 300
#' )
#'
#' find_alleles(
#'   fragments_list
#' )
#' call_repeats(
#'   fragments_list
#' )
#'
#' add_metadata(
#'   fragments_list,
#'   metadata_data.frame = trace::metadata
#' )
#'
#'assign_index_peaks(
#'   fragments_list,
#'   grouped = TRUE
#' )
#'
#' plot_traces(fragments_list[1], xlim = c(100,150))
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
      if (!is.na(fragments_list[[i]]$metrics_group_id) && fragments_list[[i]]$metrics_baseline_control == TRUE) {
        # since there can be more than one control, make a list of them
        baseline_control_list[[fragments_list[[i]]$metrics_group_id]] <- c(
          baseline_control_list[[fragments_list[[i]]$metrics_group_id]],
          list(
            list(
              fragments_list[[i]]$get_allele_peak()$allele_repeat,
              fragments_list[[i]]$repeat_table_df,
              fragments_list[[i]]$batch_run_id
            )
          )
        )
      }
    }

    # do some quality control

    for (i in seq_along(baseline_control_list)) {
      controls_missing_allele <- all(sapply(baseline_control_list[[fragments_list[[i]]$metrics_group_id]], function(x) is.na(x[[1]])))

      if (length(baseline_control_list[[i]]) == 0) {
        warning(paste0("Group '", names(baseline_control_list)[[i]], "' has no 'metrics_baseline_control'. Instability metrics won't be calculated for this group in subsequent calculations."),
          call. = FALSE
        )
      }  else if (controls_missing_allele == TRUE) {
        warning(paste0("Group '", names(baseline_control_list)[[i]], "' control has no allele called. Instability metrics won't be calculated for this group in subsequent calculations."),
          call. = FALSE
        )
      } 
    }

    # loop over each sample and put data inside
    for (i in seq_along(fragments_list)) {
      # if the group has no metrics_baseline_control it will be NULL so length == 0
      if(length(baseline_control_list[[fragments_list[[i]]$metrics_group_id]]) > 0){
              control_index_median_repeat <- median(sapply(baseline_control_list[[fragments_list[[i]]$metrics_group_id]], function(x) x[[1]]), na.rm = TRUE)
      } else{
        control_index_median_repeat <- NA_real_
      }
      # set index peak
        # samples with no data are skipped inside set_index_peaks and if NA value is provided index peak will be set to NA
      fragments_list[[i]]$set_index_peak(control_index_median_repeat)
      fragments_list[[i]]$.__enclos_env__$private$index_samples <- baseline_control_list[[fragments_list[[i]]$metrics_group_id]]
      fragments_list[[i]]$.__enclos_env__$private$assigned_index_peak_grouped <- TRUE

      # check if the index samples are from a different batch and the samples were not batch corrected
      index_sample_batch_ids <- unique(sapply(baseline_control_list[[fragments_list[[i]]$metrics_group_id]], function(x) x[[3]]))
      if(length(index_sample_batch_ids) > 0 && !fragments_list[[i]]$batch_run_id %in% index_sample_batch_ids){
        # so we've established that the index samples are from different run batch. 
        # now check if they are they were batch corrected
        if(is.na(fragments_list[[i]]$.__enclos_env__$private$batch_correction_factor)){
          warning(
            call. = FALSE,
            paste0(fragments_list[[i]]$unique_id, " was grouped for index assignment, but its 'metrics_baseline_control' appears to be from a different 'batch_run_id'. ",
              "Please run use 'batch_correction' in 'call_repeats()' to correct systematic differences between runs that may impact correct index peak assignment.")
          )
        }        
      }
    }
  } else {
    # otherwise just use the modal peak as the index peak
    fragments_list <- lapply(fragments_list, function(x) {
      x$set_index_peak(x$get_allele_peak()$allele_repeat)
      x$.__enclos_env__$private$assigned_index_peak_grouped <- FALSE
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
        x$set_index_peak(
          as.numeric(index_override_dataframe[which(index_override_dataframe[, 1] == x$unique_id), 2])
        )
      }
      return(x)
    })
  }

  invisible()
}

