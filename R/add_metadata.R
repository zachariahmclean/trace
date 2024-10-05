transfer_metadata_helper <- function(old_fragment,
                                     new_fragment) {
  metadata_names <- c(
    "unique_id",
    "metrics_group_id",
    "metrics_baseline_control",
    "batch_run_id",
    "batch_sample_id",
    "batch_sample_modal_repeat"
  )


  for (name in metadata_names) {
    eval(parse(
      text = paste0(
        "new_fragment$",
        name,
        "<- old_fragment$",
        name
      )
    ))
  }
  return(new_fragment)
}

# add metadata ------------------------------------------------------------

#' Add Metadata to Fragments List
#'
#' This function adds metadata information to a list of fragments.
#'
#' @param fragments_list A list of fragment objects to which metadata will be added.
#' @param metadata_data.frame A data frame containing the metadata information.
#' @param unique_id (required) A character string indicating the column name for unique sample identifiers in the metadata.
#' @param metrics_group_id (optional) A character string indicating the column name for sample group identifiers in the metadata. This is for the logical grouping of samples for metrics calculations (see [assign_index_peaks()]). To skip, provide NA.
#' @param metrics_baseline_control (optional) A character string indicating the column name for baseline control indicators in the metadata. This is to identify the baseline control samples for grouping of samples for metrics calculations (see [assign_index_peaks()]). To skip, provide NA.
#' @param batch_run_id (optional) A character string indicating the column name for the batch run identifiers in the metadata. This is for either batch correction or repeat correction in [call_repeats()]. To skip, provide NA.
#' @param batch_sample_id (optional) A character string indicating the column name for an id of the size standard. For example, a sample code. This is for either batch correction or repeat correction in [call_repeats()]. To skip, provide NA.
#' @param batch_sample_modal_repeat (optional) A character string indicating column name for the validated modal repeat length of size standard sample. This is for either repeat correction in [call_repeats()]. To skip, provide NA.
#'
#' @return This function modifies list of fragments objects in place with metadata added.
#'
#' @details 
#' This function adds specified metadata attributes to each fragment in the list. It matches the unique sample identifiers from the fragments list with those in the metadata data frame. To skip any of the optional columns, make parameter NA.
#'
#' There are two key things metadata are required for. First is the grouping of samples (metrics_group_id & metrics_baseline_control) for the calculation of metrics and is used in [assign_index_peaks()]. For example, specifying a sample where the modal allele is the inherited repeat length (eg a mouse tail sample) or sample(s) at the start of a time-course experiment. This is indicated with a TRUE in the metrics_baseline_control column of the metadata. Samples are then grouped together with the metrics_group_id column of the metadata. Multiple samples can be metrics_baseline_control, which can be helpful for the average repeat gain metric to have a more accurate representation of the average repeat at the start of the experiment.
#' 
#' The second key thing metadata can be used for is corrections in [call_repeats()]. There are two main correction approaches in [call_repeats()] that are somewhat related: either 'batch' or 'repeat'. Batch correction is relatively none and just requires you to link samples across batches to correct batch-batch variation in repeat sizes. However, even though the repeat size that is return will be precise, it will not be accurate and underestimates the real repeat length. By contrast, repeat correction can be used to accurately call repeat lengths (which also corrects the batch effects). However, the repeat correction will only be as good as your sample used to call the repeat length so this is a challenging and advanced feature. You need to use a sample that reliably returns the same peak as the modal peak, or you need to be willing to understand the shape of the distribution and manually validate the repeat length of each batch_sample_id for each run. 
#' 
#' Batch correction uses common sample(s) across fragment analysis runs to correct systematic batch effects that occur with repeat-containing amplicons in capillary electrophoresis. There are slight fluctuations of size across runs for amplicons containing repeats that result in systematic differences around 1-3 base pairs. So, if samples are to be analyzed for different runs, the absolute bp size is not comparable unless this batch effect is corrected. This is only relevant when the absolute size of a amplicons are compared for grouping metrics as described above (otherwise instability metrics are all relative and it doesn’t matter that there’s systematic batch effects across runs) or when plotting traces from different runs. This correction can be achieved by running a couple of samples in every fragment analysis run, or having a single run that takes a couple of samples from every run together, thereby linking them. These samples are then indicated in the metadata with batch_run_id (to group samples by fragment analysis run) and batch_sample_id (to enable linking samples across batches).
#' 
#' Finally, samples with known and validated repeat size can be used to accurately call the repeat length (and therefore also correct batch effects) in [call_repeats()]. Similar to batch correction, batch_run_id (to group samples by fragment analysis run) and batch_sample_id (to enable linking samples across batches) are used, but importantly batch_sample_modal_repeat is also set. The batch_sample_modal_repeat is the validated repeat length of the modal repeat of the sample. This validated repeat length is then used to call the repeat length of the modal repeat for each sample (by each batch_run_id). Importantly, this correction requires you to know with confidence the repeat length of the modal peak of the sample. Therefore it's important that the sample used for repeat correction has a clear and prominent modal peak. If the repeat length is very long, it's common for the modal peak of a sample to change so if you use this feature you're going to have to understand the shape of the distribution of your sample and double check that the correct peak has been called as the modal peak after you have used [find_alleles()]. If a different peak is selected as the modal peak than usual, you need to go back to the metadata and adjust the repeat size of the size standard (For example, your size standard sample has been validated to have 120 repeats. You run [find_alleles()] and look at the distribution of peaks and notice that the peak one repeat unit higher is the modal peak this time. Therefore, you're going to need to set the batch_sample_modal_repeat as 121 in the metadata just for that batch_run_id. In the other runs you would keep the batch_sample_modal_repeat as 120.).
#' 
#' @export
#'
#' @examples
#'
#' gm_raw <- trace::example_data
#' metadata <- trace::metadata
#'
#' test_fragments <- peak_table_to_fragments(gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B",
#'   min_size_bp = 300
#' )
#'
#' add_metadata(
#'   fragments_list = test_fragments,
#'   metadata_data.frame = metadata,
#'   unique_id = "unique_id",
#'   metrics_group_id = "metrics_group_id",
#'   metrics_baseline_control = "metrics_baseline_control",
#'   batch_run_id = "batch_run_id",
#'   batch_sample_id = "batch_sample_id",
#'   batch_sample_modal_repeat = "batch_sample_modal_repeat"
#' )
#'
#' # skip unwanted metadata by using NA
#'
#' add_metadata(
#'   fragments_list = test_fragments,
#'   metadata_data.frame = metadata,
#'   unique_id = "unique_id",
#'   metrics_group_id = "metrics_group_id",
#'   metrics_baseline_control = "metrics_baseline_control",
#'   batch_run_id = NA,
#'   batch_sample_id = NA,
#'   batch_sample_modal_repeat = NA
#' )
#'
add_metadata <- function(
    fragments_list,
    metadata_data.frame,
    unique_id = "unique_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = "metrics_baseline_control",
    batch_run_id = "batch_run_id",
    batch_sample_id = "batch_sample_id",
    batch_sample_modal_repeat = "batch_sample_modal_repeat") {
  # validate inputs to give good errors to user
  ## check to make sure that if the user supplies a column name, that it's actually in the dataframe they supplied
  function_input_vector <- c(
    unique_id,  metrics_group_id, metrics_baseline_control, 
    batch_run_id, batch_sample_id, batch_sample_modal_repeat
    
  )
  function_input_name_vector <- c(
    "unique_id", "metrics_group_id", "metrics_baseline_control", 
    "batch_run_id", "batch_sample_id", "batch_sample_modal_repeat"
  )
  for (i in seq_along(function_input_vector)) {
    if (!any(names(metadata_data.frame) == function_input_vector[[i]]) & !is.na(function_input_vector[[i]])) {
      stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ", function_input_name_vector[[i]], " input. If you don't want to add this metadata category, set '", function_input_name_vector[[i]], " = NA'"),
        call. = FALSE
      )
    }
  }
  ## check if user has any duplicated unique ids
  supplied_ids <- metadata_data.frame[, unique_id, drop = TRUE]
  if (anyDuplicated(supplied_ids) != 0) {
    stop(paste0(unique_id, " does not contain unique sample ids. The metadata must have one row per unique sample id."),
      call. = FALSE
    )
  }
  ## Give warning if samples don't have metadata
  not_in_metadata <- which(!names(fragments_list) %in% supplied_ids)
  if (length(not_in_metadata) > 0) {
    warning(
      paste0(
        "The following samples do not have a corresponding unique id in the metadata: ",
        paste0(names(fragments_list)[not_in_metadata], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## Give warning if user tries to give a df metadata but it's not in sample list
  not_in_samples <- which(!supplied_ids %in% names(fragments_list))
  if (length(not_in_samples) > 0) {
    warning(
      paste0(
        "The following unique ids in the metadata file do not have a corresponding sample: ",
        paste0(supplied_ids[not_in_samples], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # make sure dataframe, not tibble
  metadata_data.frame <- as.data.frame(metadata_data.frame)
  metadata_data.frame[metadata_data.frame == ''] <- NA

  metadata_added <- lapply(
    fragments_list,
    function(fragments) {
    
      # filter for row of sample
      sample_metadata <- metadata_data.frame[which(metadata_data.frame[unique_id] == fragments$unique_id), , drop = FALSE]
    
      # add metadata to slots, checking if parameters are NA
      fragments$metrics_group_id <- if (!is.na(metrics_group_id)) as.character(sample_metadata[[metrics_group_id]]) else NA_character_
      fragments$metrics_baseline_control <- if (!is.na(metrics_baseline_control)) {
        ifelse(is.na(sample_metadata[[metrics_baseline_control]]) || !as.logical(sample_metadata[[metrics_baseline_control]]), FALSE, TRUE)
      } else {
        FALSE
      }
      fragments$batch_run_id <- if (!is.na(batch_run_id)) as.character(sample_metadata[[batch_run_id]]) else NA_character_
      fragments$batch_sample_id <- if (!is.na(batch_sample_id)) as.character(sample_metadata[[batch_sample_id]]) else NA_character_
      fragments$batch_sample_modal_repeat <- if(!is.na(batch_sample_modal_repeat)) as.numeric(sample_metadata[[batch_sample_modal_repeat]]) else NA_real_
    
      return(fragments)
    }
  )

  invisible()
}

