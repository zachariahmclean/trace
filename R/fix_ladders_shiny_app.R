# modules -----------------------------------------------------------------
## select sample module
sample_selection_ui <- function(id) {
  shiny::tagList(
    shiny::selectInput(shiny::NS(id, "unique_id_selection"), "Sample selection", NULL),
    shiny::checkboxInput(shiny::NS(id, "warning_checkbox"), label = "Select only samples with ladder warnings", value = FALSE)
  )
}



sample_selection_module <- function(id, fragment_trace_list) {
  shiny::moduleServer(id, function(input, output, session) {
    ladder_warning_samples <- shiny::reactive({
      sapply(
        shiny::reactiveValuesToList(fragment_trace_list),
        function(x) {
          if (is.null(tryCatch(ladder_rsq_warning_helper(x, 0.998),
            warning = function(w) w
          ))) {
            FALSE
          } else {
            TRUE
          }
        }
      )
    })

    choices <- shiny::reactive({
      if (input$warning_checkbox) {
        names(fragment_trace_list)[which(ladder_warning_samples())]
      } else {
        names(fragment_trace_list)
      }
    })


    shiny::observe({
      shiny::updateSelectInput(session, "unique_id_selection",
        choices = choices()
      )
    })


    selected_fragments_trace <- shiny::reactive({
      # if(is.null(input$unique_id_selection)){
      if (input$unique_id_selection == "") {
        fragment_trace_list[[choices()[1]]]
      } else {
        fragment_trace_list[[input$unique_id_selection]]
      }
    })

    return(list(
      sample = selected_fragments_trace,
      input_unique_id_selection = shiny::reactive(input$unique_id_selection)
    ))
  })
}


## plot module
plot_module_ui <- function(id) {
  shiny::tagList(
    plotly::plotlyOutput(shiny::NS(id, "plot"))
  )
}

plot_module_server <- function(id, fragment_ladder, input_unique_id_selection,
                               find_scan_max) {
  shiny::moduleServer(id, function(input, output, session) {
    # Initialize ladders as NULL
    ladders <- shiny::reactiveValues(scan = NULL, size = NULL)
    relayout_data <- shiny::reactiveVal(NULL) # Initialize relayout_data

    # Reset ladders and relayout_data when unique_id_selection changes
    shiny::observeEvent(input_unique_id_selection(), {
      ladders$scan <- NULL
      ladders$size <- NULL
      relayout_data(NULL)
    })

    shiny::observe({
      ladders$scan <- fragment_ladder()$ladder_df$scan
      ladders$size <- fragment_ladder()$ladder_df$size
    })

    output$plot <- plotly::renderPlotly({
      if (is.null(ladders$scan) || is.null(ladders$size)) {
        # Return a blank plot if ladders are not initialized
        return(plotly::plot_ly())
      }

      shapes_with_labels <- list()
      text_annotations <- list()
      for (i in 1:length(ladders$scan)) {
        shapes_with_labels[[i]] <- list(
          type = "line",
          x0 = ladders$scan[i], # Adjust as needed for the positions of your shapes
          x1 = ladders$scan[i], # Adjust as needed for the positions of your shapes
          y0 = 0.05,
          y1 = 0.45,
          yref = "paper",
          fillcolor = "rgba(0,0,0,0)", # Transparent fill
          line = list(
            color = "black",
            width = 1
          ),
          editable = TRUE # Allow shape editing
        )

        # Add text annotation
        text_annotations[[i]] <- list(
          x = ladders$scan[i], # X-position of the text
          y = max(fragment_ladder()$trace_bp_df$ladder_signal) / 2, # Adjust Y-position as needed
          text = ladders$size[i],
          showarrow = FALSE, # Remove arrow if not desired
          textanchor = "end", # Horizontal text alignment
          yanchor = "middle", # Vertical text alignment
          font = list(
            color = "black",
            size = 10
          ),
          textangle = 270
        )
      }

      p <- plotly::plot_ly(fragment_ladder()$trace_bp_df, x = ~scan, y = ~ladder_signal, type = "scatter", mode = "lines")
      p <- plotly::layout(p, shapes = shapes_with_labels, annotations = text_annotations, title = fragment_ladder()$unique_id)
      # allow to edit plot by dragging lines
      plotly::config(p, edits = list(shapePosition = TRUE))
    })

    # Reset relayout_data when plot is clicked or dragged
    shiny::observeEvent(plotly::event_data("plotly_relayout"), {
      relayout_data(plotly::event_data("plotly_relayout"))
    })

    # Capture relayout_data
    shiny::observe({
      if (!is.null(relayout_data())) {
        ed <- relayout_data()
        scan_positions <- ed[grepl("^shapes.*x.*", names(ed))]
        if (length(scan_positions) != 2) {
          return()
        }
        row_index <- unique(as.numeric(sub(".*\\[(.*?)\\].*", "\\1", names(scan_positions)[1])) + 1)

        # find maximal signal in the user defined region
        selected_scan <- round(as.numeric(scan_positions))[1]
        window_df <- fragment_ladder()$trace_bp_df[which(fragment_ladder()$trace_bp_df$scan > selected_scan - find_scan_max() & fragment_ladder()$trace_bp_df$scan < selected_scan + find_scan_max()), ]
        new_scan <- window_df[which(window_df$ladder_signal == max(window_df$ladder_signal)), "scan"]

        # assign scan
        ladders$scan[row_index] <- new_scan[1]
      }
    })

    return(list(ladders = shiny::reactive(ladders)))
  })
}

## export ladder fixes

ladder_export_ui <- function(id) {
  shiny::tagList(
    shiny::downloadButton(shiny::NS(id, "download"), "Download ladder corrections")
  )
}

ladder_export_server <- function(id, manual_ladder_list) {
  shiny::moduleServer(id, function(input, output, session) {
    output$download <- shiny::downloadHandler(
      filename = function() {
        paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "_ladder_df_list", ".rds")
      },
      content = function(file) {
        saveRDS(shiny::reactiveValuesToList(manual_ladder_list), file)
      }
    )
  })
}


## r squared table


rsq_table_ui <- function(id) {
  shiny::tagList(
    shiny::tableOutput(shiny::NS(id, "rsq_table"))
  )
}

rsq_table_server <- function(id, fragment_ladder, fragment_ladder_trigger) {
  shiny::moduleServer(id, function(input, output, session) {
    rsq_table <- shiny::reactive({
      fragment_ladder_trigger()  # Trigger reactivity with fragment_ladder_trigger

      rsq <- sapply(fragment_ladder()$local_southern_mod, function(y) suppressWarnings(summary(y$mod)$r.squared))
      size_ranges <- sapply(fragment_ladder()$local_southern_mod, function(y) y$mod$model$yi)
      size_ranges_vector <- vector("numeric", ncol(size_ranges))
      for (j in seq_along(size_ranges_vector)) {
        size_ranges_vector[j] <- paste0(size_ranges[1, j], ", ", size_ranges[2, j], ", ", size_ranges[3, j])
      }

      data.frame(
        sizes = size_ranges_vector,
        r_squared = as.character(round(rsq, digits = 4))
      )
    })


    output$rsq_table <- shiny::renderTable({
      rsq_table()
    })
  })
}

# Shiny App ---------------------------------------------------------------

ui <- shiny::fluidPage(
  shiny::titlePanel("Interactive ladder fixing"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      sample_selection_ui("sample_selection"),
      shiny::sliderInput("find_scan_max", "Snap to tallest scan window",
        min = 0, max = 50, value = 10
      ),
      ladder_export_ui("ladder_df_list_download")
    ),
    shiny::mainPanel(
      plot_module_ui("plot_module"),
      rsq_table_ui("rsq_table")
    )
  )
)




###
server_function <- function(input, output, session, fragment_trace_list) {
  fragment_trace_list_reactive <- shiny::reactiveValues()
  for (sample_name in names(fragment_trace_list)) {
    fragment_trace_list_reactive[[sample_name]] <- fragment_trace_list[[sample_name]]
  }
  manual_ladder_list <- shiny::reactiveValues()

  selected_fragments_trace <- sample_selection_module("sample_selection", fragment_trace_list_reactive)


  plot_output <- plot_module_server(
    "plot_module",
    selected_fragments_trace$sample,
    selected_fragments_trace$input_unique_id_selection,
    shiny::reactive(input$find_scan_max)
  )

  # provide a reactive trigger that class has been updated
  # this is for the rsq_table_server that needs to see that selected_fragments_trace$sample has been updated
  fragment_ladder_trigger <- shiny::reactiveVal(0)

  rsq_table_server("rsq_table", selected_fragments_trace$sample, fragment_ladder_trigger)


  # have a reactive list that gets updated when you change the stuff
  shiny::observe({
    sample_unique_id <- selected_fragments_trace$sample()$unique_id

    selected_ladder_df <- selected_fragments_trace$sample()$ladder_df
    selected_sample_scans <- selected_ladder_df[which(!is.na(selected_ladder_df$size)), "scan"]

    plot_ladder_df <- as.data.frame(shiny::reactiveValuesToList(plot_output$ladders()))
    plot_scans <- plot_ladder_df[which(!is.na(plot_ladder_df$size)), "scan"]

    # skip if ladder info hasn't been updated
    if (identical(selected_sample_scans, plot_scans)) {
      return()
    } else if (nrow(as.data.frame(shiny::reactiveValuesToList(plot_output$ladders()))) == 0) {
      return()
    }

    manual_ladder_list[[sample_unique_id]] <- as.data.frame(shiny::reactiveValuesToList(plot_output$ladders()))
    fix_ladders_manual(
      shiny::reactiveValuesToList(fragment_trace_list_reactive)[sample_unique_id],
      shiny::reactiveValuesToList(manual_ladder_list)[sample_unique_id]
    )

    fragment_ladder_trigger(fragment_ladder_trigger() + 1)

  })

  # export data
  ladder_export_server("ladder_df_list_download", manual_ladder_list)
}


#' Fix ladders interactively
#'
#' An app for fixing ladders
#'
#' @param fragment_trace_list A list of fragments_trace objects containing fragment data
#'
#' @return interactive shiny app for fixing ladders
#' @export
#' @details
#' This function helps you fix ladders that are incorrectly assigned. Run `fix_ladders_interactive()`
#' and provide output from `find_ladders`. In the app, for each sample, click on
#' line for the incorrect ladder size and drag it to the correct peak.
#'
#' Once you are satisfied with the ladders for all the broken samples, click the download
#' button to generate a file that has the ladder correction data. Read this file
#' back into R using readRDS, then use [fix_ladders_manual()] and supply the ladder
#' correction data as `ladder_df_list`. This allows the manually corrected data to
#' be saved and used within a script so that the correct does not need to be done
#' every time.
#'
#' @seealso [fix_ladders_manual()], [find_ladders()]
#'
#'
#' @examples
#' fsa_list <- lapply(cell_line_fsa_list["20230413_B03.fsa"], function(x) x$clone())
#' 
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' # to create an example, lets brake one of the ladders
#' brake_ladder_list <- list(
#'   "20230413_B03.fsa" = data.frame(
#'     size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#'     scan = c(1555, 1633, 1783, 1827, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792,
#'              4085, 4322, 4370)
#'   )
#' )
#'
#' fix_ladders_manual(
#'   fsa_list,
#'   brake_ladder_list
#' )
#'
#' plot_ladders(fsa_list)
#'
#'
#' if (interactive()) {
#'   fix_ladders_interactive(test_ladders_broken)
#' }
#'
#' # once you have corrected your ladders in the app,
#' # export the data we need to incorporate that into the script:
#' # ladder_df_list <- readRDS('path/to/exported/data.rds')
#' # test_ladders_fixed <- fix_ladders_manual(test_ladders_broken, ladder_df_list)
#'
#' # plot_ladders(test_ladders_fixed["20230413_B03.fsa"],
#' #           n_facet_col = 1)
#'
fix_ladders_interactive <- function(fragment_trace_list) {
  message("To incorporate the manual corrections into your script you need to do the following:")
  message("1: read in the corrected ladder data using 'ladder_df_list <- readRDS('path/to/exported/data.rds')'")
  message("2: Run 'fix_ladders_manual(fragments_trace_list, ladder_df_list)'")

  # Launch the Shiny app with fragment_trace_list passed as a parameter
  shiny::shinyApp(
    ui = ui,
    server = function(input, output, session) {
      server_function(input, output, session, fragment_trace_list)
    }
  )
}




