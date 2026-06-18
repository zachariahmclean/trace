# metadata template -------------------------------------------------------

#' Generate a metadata template from fsa files
#'
#' Create a metadata template with one row per sample and the `batch_run_id`
#' column already populated from each fsa file, ready to be filled in and
#' passed to [trace()] (via [add_metadata()]).
#'
#' @param input One of: a list of fragments objects (from [read_fsa()]), a
#'   character vector of `.fsa` file paths, or a single directory path
#'   containing `.fsa` files (in which case the files are read with
#'   [read_fsa()]).
#' @param output_csv Optional file path. If supplied, the template is written
#'   to this csv with `write.csv()` (blank cells for values to be filled in).
#'
#' @return A data frame with one row per sample. The six columns recognised by
#'   [add_metadata()] are included: `unique_id` and `batch_run_id` are
#'   pre-filled from the fsa file, while `metrics_group_id`,
#'   `metrics_baseline_control`, `batch_sample_id`, and
#'   `batch_sample_modal_repeat` are left blank (`NA`) for the user to
#'   complete. Additional read-only `fsa_*` provenance columns (run date, well,
#'   plate, instrument, capillary, sample name) are included to help fill in
#'   the rest; these are ignored by [add_metadata()].
#'
#' @details
#' This is the recommended starting point for building a metadata file.
#' Historically users typed the fragment analysis run for each sample by hand,
#' which is error prone. Here `batch_run_id` is read directly from the fsa file
#' (the `RunN` tag, or a date/time fallback) so samples are grouped by the run
#' they were actually processed in. Do not edit `batch_run_id` unless you know
#' a run was mislabelled. Rows are ordered by run then well so the template
#' reads like the original plate.
#'
#' @export
#' @seealso [read_fsa()], [extract_fsa_metadata()], [add_metadata()]
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#' # import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
#'
#' template <- generate_metadata_template(fsa_list)
#' head(template)
#'
generate_metadata_template <- function(input, output_csv = NULL) {
  # resolve input to a list of fragments objects
  if (is.list(input) && length(input) > 0 && inherits(input[[1]], "fragments")) {
    fragments_list <- input
  } else if (is.character(input) && length(input) > 0) {
    files <- input
    if (length(files) == 1 && dir.exists(files)) {
      files <- list.files(files, pattern = "\\.fsa$", full.names = TRUE, ignore.case = TRUE)
    }
    if (length(files) == 0) {
      stop(call. = FALSE, "No .fsa files found for 'input'")
    }
    fragments_list <- read_fsa(files)
  } else {
    stop(call. = FALSE, "'input' must be a list of fragments objects (from read_fsa()), a vector of .fsa file paths, or a directory path")
  }

  meta <- extract_fsa_metadata(fragments_list)

  template <- data.frame(
    unique_id = meta$unique_id,
    metrics_group_id = NA_character_,
    metrics_baseline_control = NA,
    batch_run_id = meta$run_id,
    batch_sample_id = NA_character_,
    batch_sample_modal_repeat = NA_real_,
    fsa_run_date = meta$run_datetime,
    fsa_well = meta$well,
    fsa_plate = meta$plate,
    fsa_instrument = meta$instrument,
    fsa_capillary = meta$capillary,
    fsa_sample_name = meta$sample_name,
    stringsAsFactors = FALSE
  )

  # order by run then well so the template reads like the plate
  template <- template[order(template$batch_run_id, template$fsa_well), , drop = FALSE]
  rownames(template) <- NULL

  if (!is.null(output_csv)) {
    utils::write.csv(template, output_csv, row.names = FALSE, na = "")
    message("Wrote metadata template to ", output_csv)
  }

  template
}
