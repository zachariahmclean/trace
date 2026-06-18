# fsa metadata parsing -----------------------------------------------------

#' Parse header fields from an ABIF (fsa) object
#'
#' Internal helper that extracts provenance and QC-relevant fields from the
#' ABIF directory returned by [seqinr::read.abif()]. Missing tags yield `NA`
#' (or length-zero) rather than an error, so files from different instruments
#' are tolerated.
#'
#' @param abif An ABIF object (the `fsa` field of a fragments object).
#'
#' @return A named list with elements: `run_id`, `run_name`, `run_datetime`,
#'   `instrument`, `instrument_serial`, `plate`, `well`, `capillary`,
#'   `sample_name`, `dye_names`, `n_offscale`, `n_saturated`.
#'
#' @details
#' `run_id` is taken from the `RunN` tag (a ready-made unique run name, e.g.
#' "Run_MRSPHILLIPS2_2023-04-14_11-30_0063"). If `RunN` is absent it falls back
#' to a zero-padded "YYYYMMDD_HHMMSS" built from the `RUND`/`RUNT` tags.
#'
#' @keywords internal
parse_fsa_metadata <- function(abif) {
  if (is.null(abif) || is.null(abif$Data)) {
    return(.empty_fsa_metadata())
  }
  d <- abif$Data
  list(
    run_id            = .fsa_run_id(d),
    run_name          = .fsa_chr(d$RunN.1),
    run_datetime      = .fsa_datetime(d),
    instrument        = .fsa_chr(d$MCHN.1),
    instrument_serial = .fsa_serial(d$HCFG.4),
    plate             = .fsa_chr(d$CTNM.1),
    well              = .fsa_chr(d$TUBE.1),
    capillary         = .fsa_num(d$LANE.1),
    sample_name       = .fsa_chr(d$SpNm.1),
    dye_names         = .fsa_dye_names(d),
    n_offscale        = .fsa_len(d$OfSc.1),
    n_saturated       = .fsa_len(d$Satd.1)
  )
}

.empty_fsa_metadata <- function() {
  list(
    run_id = NA_character_, run_name = NA_character_, run_datetime = NA_character_,
    instrument = NA_character_, instrument_serial = NA_character_, plate = NA_character_,
    well = NA_character_, capillary = NA_real_, sample_name = NA_character_,
    dye_names = character(0), n_offscale = 0L, n_saturated = 0L
  )
}

# coerce first element to a clean string, treating empty/blank as NA
.fsa_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  v <- trimws(as.character(x[[1]]))
  if (length(v) == 0 || is.na(v) || !nzchar(v)) NA_character_ else v
}

.fsa_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x[[1]]))
}

.fsa_len <- function(x) if (is.null(x)) 0L else length(x)

# pull SerialNumber=... out of the HCFG hardware-config string
.fsa_serial <- function(x) {
  s <- .fsa_chr(x)
  if (is.na(s) || !grepl("SerialNumber=", s)) return(NA_character_)
  trimws(sub(".*SerialNumber=([^;]+).*", "\\1", s))
}

# RUND is a list (year, month, day); RUNT is a list (hour, min, sec, ...)
.fsa_date_parts <- function(d) suppressWarnings(as.integer(unlist(d$RUND.1)))
.fsa_time_parts <- function(d) suppressWarnings(as.integer(unlist(d$RUNT.1)))

.fsa_datetime <- function(d) {
  dt <- .fsa_date_parts(d)
  tm <- .fsa_time_parts(d)
  if (length(dt) < 3 || length(tm) < 3 || anyNA(dt[1:3]) || anyNA(tm[1:3])) return(NA_character_)
  sprintf("%04d-%02d-%02d %02d:%02d:%02d", dt[1], dt[2], dt[3], tm[1], tm[2], tm[3])
}

.fsa_run_id <- function(d) {
  runn <- .fsa_chr(d$RunN.1)
  if (!is.na(runn)) return(runn)
  dt <- .fsa_date_parts(d)
  tm <- .fsa_time_parts(d)
  if (length(dt) < 3 || length(tm) < 3 || anyNA(dt[1:3]) || anyNA(tm[1:3])) return(NA_character_)
  sprintf("%04d%02d%02d_%02d%02d%02d", dt[1], dt[2], dt[3], tm[1], tm[2], tm[3])
}

# named vector of dye names, ordered by channel (e.g. dye_1 = "6-FAM", dye_5 = "LIZ")
.fsa_dye_names <- function(d) {
  nm <- grep("^DyeN\\.[0-9]+$", names(d), value = TRUE)
  if (!length(nm)) return(character(0))
  nm <- nm[order(as.integer(sub("^DyeN\\.", "", nm)))]
  v <- vapply(nm, function(n) .fsa_chr(d[[n]]), character(1))
  names(v) <- sub("^DyeN\\.", "dye_", nm)
  v
}

# pick the run id for a parsed fsa_metadata list according to the configured tag
.run_id_for_tag <- function(m, tag) {
  if (identical(tag, "RUND_RUNT")) {
    if (is.na(m$run_datetime)) return(NA_character_)
    # "2023-04-14 11:30:50" -> "20230414_113050"
    parts <- strsplit(m$run_datetime, " ", fixed = TRUE)[[1]]
    paste0(gsub("-", "", parts[1]), "_", gsub(":", "", parts[2]))
  } else {
    m$run_id
  }
}

#' Set batch_run_id from fsa metadata
#'
#' Internal step that derives `batch_run_id` for fsa-imported samples directly
#' from the fsa file, so the fragment analysis run does not have to be entered
#' by hand. The fsa file is treated as authoritative: if a sample already has a
#' user-supplied `batch_run_id` that disagrees with the fsa-derived value, the
#' fsa value is used and a warning lists the mismatches so run mix-ups are
#' caught.
#'
#' @param fragments_list A list of fragments objects.
#' @param config A trace_config object (uses `batch_run_id_tag`).
#'
#' @return A trace_output status object (warnings for any mismatches). Modifies
#'   the fragments objects in place.
#'
#' @keywords internal
set_batch_run_id_from_fsa <- function(fragments_list, config) {
  output <- trace_output$new("set_batch_run_id_from_fsa")
  tag <- if (!is.null(config$batch_run_id_tag)) config$batch_run_id_tag else "RunN"

  mismatches <- character()
  for (x in fragments_list) {
    if (!identical(x$input_method, "fsa")) next

    m <- x$fsa_metadata
    if (is.null(m)) m <- parse_fsa_metadata(x$fsa)

    fsa_run <- .run_id_for_tag(m, tag)
    if (is.na(fsa_run)) next

    if (!is.na(x$batch_run_id) && x$batch_run_id != fsa_run) {
      mismatches <- c(
        mismatches,
        sprintf("%s (metadata='%s', fsa='%s')", x$unique_id, x$batch_run_id, fsa_run)
      )
    }
    x$batch_run_id <- fsa_run
  }

  if (length(mismatches) > 0) {
    output$set_status(
      "warning",
      paste0(
        "User-supplied batch_run_id disagreed with the fsa-derived run id (the fsa value was used). ",
        "Check for run mix-ups:\n  ",
        paste(mismatches, collapse = "\n  ")
      )
    )
  }

  output
}
