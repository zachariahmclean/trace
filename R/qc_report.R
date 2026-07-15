# qc_report ---------------------------------------------------------------

#' Sample quality control report
#'
#' Summarise per-sample quality control metrics for a processed list of
#' fragments objects, flagging samples that may need attention before
#' calculating instability metrics.
#'
#' @param fragments_list A list of fragments objects that have been processed
#'   (e.g. with [trace()]).
#' @param config_file Optional file path to a YAML config (see [trace()]).
#' @param ... QC threshold parameters that overwrite the config, namely
#'   `qc_min_rsq`, `qc_min_peaks`, `qc_min_modal_signal`,
#'   `qc_saturation_ceiling`, `qc_window`, and `qc_prominence_min`.
#'
#' @return A data frame with one row per sample containing provenance
#'   (`unique_id`, `batch_run_id`, `run_date`, `instrument`, `capillary`,
#'   `well`, `plate`), quality metrics (`ladder_avg_rsq`, `ladder_min_rsq`,
#'   `n_peaks`, `modal_size`, `modal_signal`, `modal_prominence`,
#'   `modal_saturated`, `saturation_in_window`), and a `qc_flags` column
#'   (a `;`-separated list of tripped checks) plus a logical `qc_pass`.
#'
#' @details
#' No single metric identifies every problematic sample, so this is a panel of
#' checks and a sample fails if it trips any of them:
#'   \itemize{
#'     \item `low_ladder_rsq`: worst ladder segment R-squared (`ladder_min_rsq`)
#'       below `qc_min_rsq`. The worst segment is used rather than the average
#'       because a single broken region can be hidden by a good average.
#'     \item `few_peaks`: fewer than `qc_min_peaks` peaks detected (a failed or
#'       empty sample).
#'     \item `no_modal_peak`: no modal/allele peak could be identified
#'       (`modal_size` or `modal_signal` is `NA`), regardless of `n_peaks`.
#'     \item `low_signal`: modal peak signal below `qc_min_modal_signal`.
#'     \item `saturated_modal`: the modal peak is off-scale or at/above
#'       `qc_saturation_ceiling`, so peak heights are unreliable.
#'     \item `saturation_in_window`: one or more off-scale peaks within
#'       `qc_window` bp of the modal peak.
#'     \item `low_prominence`: modal peak prominence (modal signal divided by
#'       the mean of its immediate neighbours) below `qc_prominence_min`.
#'   }
#'
#' Note that modal peak prominence decreases as repeat length increases (long
#' repeats have broad distributions with no single dominant peak), so a
#' low-prominence flag is expected for very long repeats and should be
#' interpreted relative to the repeat size rather than as an absolute failure.
#'
#' @export
#' @seealso [trace()], [extract_ladder_summary()], [extract_fsa_metadata()]
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#' # import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
#'
#' processed <- trace(fsa_list)
#' qc_report(processed)
#'
qc_report <- function(fragments_list, config_file = NULL, ...) {
  config <- load_config()
  if (!is.null(config_file)) {
    config <- update_config(config, load_config(config_file))
  }
  config <- update_config(config, list(...))

  rows <- lapply(fragments_list, function(x) {
    # ladder quality (only when a ladder was fit)
    avg_rsq <- NA_real_
    min_rsq <- NA_real_
    if (!is.null(x$ladder_df)) {
      rsq <- tryCatch(
        sapply(ladder_fit_cor(x), function(z) z$rsq),
        error = function(e) NA_real_
      )
      if (length(rsq) > 0 && !all(is.na(rsq))) {
        avg_rsq <- mean(rsq)
        min_rsq <- min(rsq)
      }
    }

    ap <- tryCatch(x$get_allele_peak(), error = function(e) NULL)
    modal_size <- if (!is.null(ap)) ap$allele_size else NA_real_
    modal_signal <- if (!is.null(ap)) ap$allele_signal else NA_real_

    pt <- x$peak_table_df
    n_peaks <- if (is.null(pt)) 0L else nrow(pt)

    modal_prominence <- NA_real_
    saturation_in_window <- NA_integer_
    modal_off_scale <- NA

    if (!is.null(pt) && nrow(pt) > 0 && !is.na(modal_size)) {
      pt <- pt[order(pt$size), ]

      below <- pt$signal[pt$size < modal_size & pt$size > modal_size - 4]
      above <- pt$signal[pt$size > modal_size & pt$size < modal_size + 4]
      neighbours <- c(
        if (length(below) > 0) below[length(below)] else NULL,
        if (length(above) > 0) above[1] else NULL
      )
      if (length(neighbours) > 0 && !is.na(modal_signal)) {
        modal_prominence <- modal_signal / mean(neighbours, na.rm = TRUE)
      }

      if ("off_scale" %in% names(pt)) {
        near <- pt[abs(pt$size - modal_size) < config$qc_window, , drop = FALSE]
        saturation_in_window <- sum(near$off_scale, na.rm = TRUE)
        modal_row <- pt[which.min(abs(pt$size - modal_size)), , drop = FALSE]
        modal_off_scale <- isTRUE(as.logical(modal_row$off_scale))
      }
    }

    modal_saturated <- isTRUE(modal_off_scale) ||
      (!is.na(modal_signal) && modal_signal >= config$qc_saturation_ceiling)

    # build flags
    flags <- character()
    if (!is.na(min_rsq) && min_rsq < config$qc_min_rsq) flags <- c(flags, "low_ladder_rsq")
    if (n_peaks < config$qc_min_peaks) flags <- c(flags, "few_peaks")
    if (is.na(modal_size) || is.na(modal_signal)) flags <- c(flags, "no_modal_peak")
    if (!is.na(modal_signal) && modal_signal < config$qc_min_modal_signal) flags <- c(flags, "low_signal")
    if (isTRUE(modal_saturated)) flags <- c(flags, "saturated_modal")
    if (!is.na(saturation_in_window) && saturation_in_window > 0) flags <- c(flags, "saturation_in_window")
    if (!is.na(modal_prominence) && modal_prominence < config$qc_prominence_min) flags <- c(flags, "low_prominence")

    m <- x$fsa_metadata
    if (is.null(m)) m <- parse_fsa_metadata(x$fsa)

    data.frame(
      unique_id = x$unique_id,
      batch_run_id = x$batch_run_id,
      run_date = m$run_datetime,
      instrument = m$instrument,
      capillary = m$capillary,
      well = m$well,
      plate = m$plate,
      ladder_avg_rsq = avg_rsq,
      ladder_min_rsq = min_rsq,
      n_peaks = n_peaks,
      modal_size = modal_size,
      modal_signal = modal_signal,
      modal_prominence = modal_prominence,
      modal_saturated = modal_saturated,
      saturation_in_window = saturation_in_window,
      qc_flags = paste(flags, collapse = ";"),
      qc_pass = length(flags) == 0,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
