# helper functions for dashboard ---------------------------------------------

.dashboard_trace_plot <- function(fragment, x_axis, show_peaks, xlim = NULL, ylim = NULL) {
  if (is.null(fragment$trace_bp_df)) {
    return(plotly::plotly_empty(type = "scatter", mode = "lines"))
  }

  data <- fragment$trace_bp_df
  has_repeats <- !is.null(fragment$repeat_table_df)

  if (is.null(x_axis) && !has_repeats) {
    data$x <- data$size
    x_axis_label <- "Size"
  } else if (is.null(x_axis) && has_repeats) {
    data$x <- data$calculated_repeats
    x_axis_label <- "Repeats"
  } else if (x_axis == "size") {
    data$x <- data$size
    x_axis_label <- "Size"
  } else {
    data$x <- data$calculated_repeats
    x_axis_label <- "Repeats"
  }

  p <- plotly::plot_ly(
    data,
    x = ~x, y = ~signal,
    type = "scatter", mode = "lines",
    line = list(color = "black", width = 1),
    hoverinfo = "x+y",
    height = 260
  )

  shapes <- list()

  if (any(data$off_scale)) {
    shapes <- c(shapes, lapply(data$x[data$off_scale], function(v) {
      list(
        type = "line", x0 = v, x1 = v, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "rgba(217,83,79,0.35)", width = 1)
      )
    }))
  }

  # index peak, drawn as a dotted line behind the trace when it's been assigned
  # and this plot's x-axis is in repeat units (index_repeat is only meaningful there)
  index_axis_repeats <- has_repeats && (is.null(x_axis) || x_axis != "size")
  if (index_axis_repeats) {
    index_repeat <- tryCatch(fragment$get_index_peak()$index_repeat, error = function(e) NA_real_)
    if (!is.na(index_repeat)) {
      shapes <- c(shapes, list(list(
        type = "line", x0 = index_repeat, x1 = index_repeat, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "black", width = 1.5, dash = "dot"),
        layer = "below"
      )))
    }
  }

  if (length(shapes) > 0) {
    p <- plotly::layout(p, shapes = shapes)
  }

  if (show_peaks) {
    peak_table <- if (has_repeats) fragment$repeat_table_df else fragment$peak_table_df
    if (!is.null(peak_table) && nrow(peak_table) > 0) {
      peak_table$x <- if (has_repeats) peak_table$repeats else peak_table$size
      p <- plotly::add_markers(p,
        data = peak_table, x = ~x, y = ~signal,
        marker = list(color = "blue", size = 5), inherit = FALSE
      )

      allele_peak <- tryCatch(fragment$get_allele_peak(), error = function(e) NULL)
      if (!is.null(allele_peak) && !is.na(allele_peak$allele_signal)) {
        allele_x <- if (has_repeats) allele_peak$allele_repeat else allele_peak$allele_size
        p <- plotly::add_markers(p,
          x = allele_x, y = allele_peak$allele_signal,
          marker = list(color = "green", size = 8), inherit = FALSE
        )
      }
    }
  }

  p <- plotly::layout(p,
    showlegend = FALSE,
    xaxis = list(title = x_axis_label, range = xlim),
    yaxis = list(title = "Signal", range = ylim),
    margin = list(t = 10, b = 30, l = 45, r = 10)
  )
  plotly::config(p, displaylogo = FALSE)
}


.dashboard_format_value <- function(x) {
  if (is.numeric(x) && !is.na(x)) {
    return(format(round(x, 3), big.mark = ""))
  }
  if (is.logical(x)) {
    return(as.character(x))
  }
  as.character(x)
}


.dashboard_qc_table <- function(qc) {
  header <- htmltools::tags$thead(
    htmltools::tags$tr(lapply(names(qc), htmltools::tags$th))
  )

  body_rows <- lapply(seq_len(nrow(qc)), function(i) {
    row_class <- if (isTRUE(qc$qc_pass[i])) NULL else "qc-fail-row"
    cells <- lapply(qc[i, ], .dashboard_format_value)
    htmltools::tags$tr(class = row_class, lapply(cells, htmltools::tags$td))
  })

  htmltools::tags$table(class = "qc-table", header, htmltools::tags$tbody(body_rows))
}


.dashboard_legend <- function() {
  legend_item <- function(swatch, label) {
    htmltools::tags$span(class = "legend-item", swatch, label)
  }

  htmltools::tags$div(
    class = "dashboard-legend",
    legend_item(htmltools::tags$span(class = "legend-dot", style = "background: blue;"), "Called peak"),
    legend_item(htmltools::tags$span(class = "legend-dot", style = "background: green;"), "Modal / allele peak"),
    legend_item(htmltools::tags$span(class = "legend-line"), "Index peak")
  )
}


.dashboard_tile <- function(fragment, qc_row, x_axis, show_peaks, xlim = NULL, ylim = NULL) {
  tile_class <- if (isTRUE(qc_row$qc_pass)) "dashboard-tile" else "dashboard-tile qc-fail"

  htmltools::tags$div(
    class = tile_class,
    htmltools::tags$h4(fragment$unique_id),
    if (!isTRUE(qc_row$qc_pass)) htmltools::tags$p(class = "qc-flags", qc_row$qc_flags),
    .dashboard_trace_plot(fragment, x_axis = x_axis, show_peaks = show_peaks, xlim = xlim, ylim = ylim)
  )
}


# Main dashboard function -----------------------------------------------------

#' Generate a static HTML trace dashboard
#'
#' Build a single self-contained-ish HTML page showing every sample's trace as
#' a small, zoomable/hoverable plot in a grid, alongside a [qc_report()]
#' summary table, so you can quickly scan a whole batch for anything that
#' looks visually wrong.
#'
#' @param fragments_list A list of fragments objects that have been processed
#'   (e.g. with [trace()]), and must have trace data (i.e. come from fsa files).
#' @param output_file File path for the output HTML file. Defaults to
#'   `"trace_dashboard.html"` in the current working directory.
#' @param config_file Optional file path to a YAML config, passed to
#'   [qc_report()] for QC thresholds (see [trace()]).
#' @param sample_subset A character vector of unique ids for a subset of
#'   samples to plot.
#' @param x_axis A character indicating what should be plotted on the x-axis,
#'   chose between `size` or `repeats`. If neither is selected, an assumption
#'   is made per-sample based on whether repeats have been called.
#' @param xlim the x limits applied to every trace plot. A numeric vector of
#'   length two.
#' @param ylim the y limits applied to every trace plot. A numeric vector of
#'   length two.
#' @param show_peaks If peak data are available, TRUE will plot the peaks on
#'   top of the trace (blue dots, with the modal/allele peak in green). A
#'   dotted vertical line is drawn behind the trace at the index peak
#'   whenever one has been assigned (see [assign_index_peaks()]), regardless
#'   of `show_peaks`. A legend for these is shown above the trace grid.
#' @param n_col A numeric value indicating the number of columns in the trace
#'   grid.
#' @param open_browser If TRUE, open the generated HTML file in a browser
#'   after it's written. Defaults to `interactive()`.
#' @param ... QC threshold parameters passed on to [qc_report()], namely
#'   `qc_min_rsq`, `qc_min_peaks`, `qc_min_modal_signal`,
#'   `qc_saturation_ceiling`, `qc_window`, and `qc_prominence_min`.
#'
#' @return Invisibly, the `output_file` path.
#'
#' @details
#' Samples that fail [qc_report()] are sorted to the front of the grid and
#' outlined in red with their tripped QC flags shown, so problems are easy to
#' spot without having to scan every tile equally.
#'
#' The traces are interactive (plotly) so you can zoom/pan/hover for exact
#' values, but the page itself is static (no R process needs to stay running
#' to view it) and can be opened directly in a browser or shared as a file.
#' Because Pandoc is not used, the traces are saved with their supporting
#' JavaScript/CSS libraries in a folder named `<output_file>_files` next to
#' the HTML file rather than embedded as one fully self-contained file --
#' keep the two together when copying or sharing.
#'
#' This is intended for visually scanning a batch, not for publication
#' figures; for that see [plot_traces()].
#'
#' @export
#' @seealso [qc_report()], [plot_traces()], [trace()]
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
#' # import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
#' processed <- trace(fsa_list)
#'
#' dashboard_file <- tempfile(fileext = ".html")
#' generate_dashboard(processed, output_file = dashboard_file, open_browser = FALSE)
#'
generate_dashboard <- function(
    fragments_list,
    output_file = "trace_dashboard.html",
    config_file = NULL,
    sample_subset = NULL,
    x_axis = NULL,
    xlim = NULL,
    ylim = NULL,
    show_peaks = TRUE,
    n_col = 3,
    open_browser = interactive(),
    ...) {
  if (!is.null(sample_subset)) {
    fragments_list <- fragments_list[which(names(fragments_list) %in% sample_subset)]
  }

  if (length(fragments_list) == 0) {
    stop(call. = FALSE, "No samples to plot after applying 'sample_subset'.")
  }

  qc <- qc_report(fragments_list, config_file = config_file, ...)

  # sort failed-QC samples first so problems are easy to spot
  ord <- order(qc$qc_pass, qc$unique_id)
  fragments_list <- fragments_list[ord]
  qc <- qc[ord, ]

  tiles <- lapply(seq_along(fragments_list), function(i) {
    tryCatch(
      .dashboard_tile(fragments_list[[i]], qc[i, ], x_axis = x_axis, show_peaks = show_peaks, xlim = xlim, ylim = ylim),
      error = function(e) {
        htmltools::tags$div(
          class = "dashboard-tile qc-fail",
          htmltools::tags$h4(fragments_list[[i]]$unique_id),
          htmltools::tags$p(class = "qc-flags", paste("Error plotting trace:", e$message))
        )
      }
    )
  })

  css <- sprintf("
    body { font-family: sans-serif; margin: 20px; }
    .qc-table { border-collapse: collapse; font-size: 0.85em; margin-bottom: 20px; }
    .qc-table th, .qc-table td { border: 1px solid #ddd; padding: 4px 8px; text-align: right; }
    .qc-table thead th { position: sticky; top: 0; background: #f7f7f7; z-index: 2; box-shadow: 0 1px 0 #ccc; }
    .qc-fail-row { background: #fff5f5; }
    .dashboard-grid { display: grid; grid-template-columns: repeat(%d, 1fr); gap: 12px; }
    .dashboard-tile { border: 1px solid #ccc; border-radius: 4px; padding: 6px; }
    .dashboard-tile h4 { margin: 2px 0; font-size: 0.9em; }
    .dashboard-tile.qc-fail { border: 2px solid #d9534f; background: #fffafa; }
    .qc-flags { color: #d9534f; font-size: 0.8em; margin: 2px 0; }
    .dashboard-legend { display: flex; gap: 20px; align-items: center; margin: 4px 0 16px; font-size: 0.85em; }
    .legend-item { display: inline-flex; align-items: center; gap: 6px; }
    .legend-dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%%; }
    .legend-line { display: inline-block; width: 24px; height: 0; border-top: 2px dotted black; }
  ", n_col)

  page <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$title("Trace dashboard"),
      htmltools::tags$style(htmltools::HTML(css))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("Trace dashboard"),
      htmltools::tags$p(sprintf(
        "%d samples, %d failed QC.",
        nrow(qc), sum(!qc$qc_pass)
      )),
      .dashboard_qc_table(qc),
      htmltools::tags$h2("Traces"),
      .dashboard_legend(),
      htmltools::tags$div(class = "dashboard-grid", tiles)
    )
  )

  htmltools::save_html(page, file = output_file, background = "white")

  if (isTRUE(open_browser)) {
    utils::browseURL(output_file)
  }

  invisible(output_file)
}
