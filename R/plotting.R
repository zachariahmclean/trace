# helper functions for class ------------------------------------------------------

plot_ladder_helper <- function(fragments_trace,
                               xlim, ylim,
                               plot_title) {
  plot(fragments_trace$trace_bp_df$scan, fragments_trace$trace_bp_df$ladder_signal,
    xlab = "Scan", ylab = "Ladder Signal",
    main = ifelse(is.null(plot_title), fragments_trace$unique_id, plot_title),
    type = "l",
    xlim = xlim,
    ylim = ylim
  )

  # Adding text
  text(fragments_trace$ladder_df$scan, rep(max(fragments_trace$trace_bp_df$ladder_signal) / 3, nrow(fragments_trace$ladder_df)),
    labels = fragments_trace$ladder_df$size,
    adj = 0.5, cex = 0.7, srt = 90
  )

  # Adding vertical lines with transparency
  for (i in 1:nrow(fragments_trace$ladder_df)) {
    abline(
      v = fragments_trace$ladder_df$scan[i],
      lty = 3,
      col = rgb(1, 0, 0, alpha = 0.3)
    )
  }
}


plot_fragments_helper <- function(fragment_repeats,
                                  ylim,
                                  xlim,
                                  plot_title) {
  if (is.null(fragment_repeats$repeat_table_df)) {
    data <- fragment_repeats$peak_table_df
    data$x <- data$size
  } else {
    data <- fragment_repeats$repeat_table_df
    data$x <- data$repeats
  }

  if (nrow(data) == 0) {
    plot.new()
    title(main = fragment_repeats$unique_id)
    return()
  }

  if (!is.null(xlim)) {
    if (length(xlim == 2) & is(xlim, "numeric")) {
      data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
    } else {
      stop(
        call. = FALSE,
        "xlim must be a numeric vector with length of 2"
      )
    }
  }


  allele_mode <- ifelse(is.null(fragment_repeats$repeat_table_df), round(fragment_repeats$get_allele_peak()$allele_size), round(fragment_repeats$get_allele_peak()$allele_repeat))

  # Fill missing y values with zeros
  rounded_x <- round(data$x)
  all_x_values <- seq(min(rounded_x), max(rounded_x))
  y_values <- rep(0, length(all_x_values))
  for (i in seq_along(rounded_x)) {
    y_values[which(all_x_values == rounded_x[i])] <- data[which(data$x == data$x[i]), "height"]
  }



  barplot(
    names.arg = all_x_values,
    height = y_values,
    main = ifelse(is.null(plot_title), fragment_repeats$unique_id, plot_title),
    xlab = ifelse(is.null(fragment_repeats$repeat_table_df), "Size", "Repeat"),
    ylab = "Signal",
    ylim = ylim,
    beside = TRUE,
    col = sapply(all_x_values, function(x) if (!is.na(allele_mode) && x == allele_mode) "red" else "gray")
  )
}


plot_trace_helper <- function(fragments,
                              show_peaks,
                              x_axis,
                              ylim,
                              xlim,
                              height_color_threshold,
                              plot_title) {
  if (is.null(fragments$trace_bp_df)) {
    stop(
      call. = FALSE,
      paste(fragments$unique_id, "This sample does not have trace data. Use fsa files as inputs to pipeline to plot trace.")
    )
  }

  # there must be a simpler way of the following if else below
  if (is.null(x_axis) && is.null(fragments$repeat_table_df)) {
    data <- fragments$trace_bp_df
    data$x <- data$size
    x_axis_label <- "Size"
  } else if (is.null(x_axis) && !is.null(fragments$repeat_table_df)) {
    data <- fragments$trace_bp_df
    data$x <- data$calculated_repeats
    x_axis_label <- "Repeats"
  } else if (x_axis == "size") {
    data <- fragments$trace_bp_df
    data$x <- data$size
    x_axis_label <- "Size"
  } else {
    data <- fragments$trace_bp_df
    data$x <- data$calculated_repeats
    x_axis_label <- "Repeats"
  }

  if (!is.null(xlim)) {
    data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
  }

  plot(data$x,
    data$signal,
    main = ifelse(is.null(plot_title), fragments$unique_id, plot_title),
    type = "l",
    xlab = x_axis_label,
    ylab = "Signal",
    ylim = ylim
  )


  if (any(data$off_scale)) {
    abline(v = data[which(data$off_scale), "x"], col = adjustcolor("red", alpha.f = 0.3), lwd = 2.5)
  }

  # add points onto plot showing peaks
  if (!is.null(fragments$peak_table_df) && show_peaks) {
    if (is.null(x_axis) && is.null(fragments$repeat_table_df)) {
      peak_table <- fragments$peak_table_df
      peak_table$x <- peak_table$size
    } else if (is.null(x_axis) && !is.null(fragments$repeat_table_df)) {
      peak_table <- fragments$repeat_table_df
      peak_table$x <- peak_table$repeats
    } else if (x_axis == "size") {
      peak_table <- fragments$peak_table_df
      peak_table$x <- peak_table$size
    } else {
      peak_table <- fragments$repeat_table_df
      peak_table$x <- peak_table$repeats
    }

    # exit early if the peak table is empty
    if (nrow(peak_table) == 0) {
      return()
    }

    if (!is.null(xlim)) {
      peak_table <- peak_table[which(peak_table$x < xlim[2] & peak_table$x > xlim[1]), ]
    }

    tallest_peak_height <- peak_table[which(peak_table$height == max(peak_table$height)), "height"]
    tallest_peak_x <- peak_table[which(peak_table$height == tallest_peak_height), "x"]
    if (!is.null(fragments$get_allele_peak()$allele_height) && !is.na(fragments$get_allele_peak()$allele_height)) {
      tallest_peak_height <- fragments$get_allele_peak()$allele_height
      # find the tallest peak x axis position
      if (is.null(x_axis) && is.na(fragments$get_allele_peak()$allele_repeat)) {
        tallest_peak_x <- fragments$get_allele_peak()$allele_size
      } else if (is.null(x_axis) && !is.na(fragments$get_allele_peak()$allele_repeat)) {
        tallest_peak_x <- fragments$get_allele_peak()$allele_repeat
      } else if (x_axis == "size") {
        tallest_peak_x <- fragments$get_allele_peak()$allele_size
      } else {
        tallest_peak_x <- fragments$get_allele_peak()$allele_repeat
      }
    }

    peaks_above <- peak_table[which(peak_table$height > tallest_peak_height * height_color_threshold), ]
    peaks_below <- peak_table[which(peak_table$height < tallest_peak_height * height_color_threshold), ]

    # Adding peaks
    points(peaks_above$x,
      peaks_above$height,
      col = "blue"
    )
    points(peaks_below$x,
      peaks_below$height,
      col = "purple"
    )
    points(tallest_peak_x,
      tallest_peak_height,
      col = "green"
    )

    # Draw horizontal dotted lines to connect repeats to their actual place on the plot
    if (!is.null(peak_table$repeats) && !is.null(peak_table$calculated_repeats)) {
      for (i in 1:nrow(peak_table)) {
        segments(
          x0 = peak_table$repeats[i],
          y0 = peak_table$height[i],
          x1 = peak_table$calculated_repeats[i],
          y1 = peak_table$height[i],
          lty = 2
        )
      }
    }
  }


  if (!is.null(fragments$get_index_peak()$index_repeat) && !is.na(fragments$get_index_peak()$index_repeat)) {
    abline(v = fragments$get_index_peak()$index_repeat, col = "black", lwd = 2, lty = 3)
  }
}


plot_data_channels_helper <- function(fragment){

  # Extract the names of the data channels that match "DATA."
  data_channels <- names(fragment$fsa$Data)[grep("DATA.", names(fragment$fsa$Data))]
  raw_data_list <- fragment$fsa$Data[data_channels]
    
  colors <- grDevices::hcl.colors(length(data_channels), palette = "Viridis")

  # Initialize the plot with the first data channel
  plot(raw_data_list[[1]], 
      type = "l", col = colors[1], lwd = 2.5, 
      ylim = range(raw_data_list),
      ylab = "Signal", xlab = "Scan", 
      main = fragment$unique_id)

  # Overlay the rest of the data channels
  for (i in 2:length(data_channels)) {
    graphics::lines(raw_data_list[[i]], col = colors[i], lwd = 2)
  }

  # Add a legend to the plot
  graphics::legend("topright", 
        legend = data_channels, 
        col = colors, 
        lty = 1, lwd = 2.5,
        cex = 0.8,        # reduce the size of the text
        ncol = 2,         # Display the legend in 2 columns
        inset = c(0.1, 0),  # Adjust position to be inside the plot, but in the corner
        x.intersp = 0.25,   # Adjust horizontal spacing between symbols and text
        y.intersp = 0.6,   # Adjust vertical spacing to align text under the lines
        text.width = 100,  #make columns closer for some reason
        bty = "n")        # Remove the box around the legend

}



# Main plotting functions -------------------------------------------------------------

#' Plot ladder
#'
#' Plot the ladder signal
#'
#' @param fragments_trace_list A list of fragments_trace objects containing fragment data.
#' @param n_facet_col A numeric value indicating the number of columns for faceting in the plot.
#' @param sample_subset A character vector of unique ids for a subset of samples to plot
#' @param xlim the x limits of the plot. A numeric vector of length two.
#' @param ylim the y limits of the plot. A numeric vector of length two.
#'
#' @return a plot of ladders
#' @export
#'
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' # Manually inspect the ladders
#' plot_ladders(fsa_list[1])
#'
plot_ladders <- function(
    fragments_trace_list,
    n_facet_col = 1,
    sample_subset = NULL,
    xlim = NULL,
    ylim = NULL) {
  if (!is.null(sample_subset)) {
    fragments_trace_list <- fragments_trace_list[which(names(fragments_trace_list) %in% sample_subset)]
  }

  graphics::par(mfrow = c(ceiling(length(fragments_trace_list) / n_facet_col), n_facet_col)) # Adjust layout as needed
  for (i in seq_along(fragments_trace_list)) {
    fragments_trace_list[[i]]$plot_ladder(
      xlim = xlim,
      ylim = ylim
    )
  }
  graphics::par(mfrow = c(1, 1)) # Reset the layout
}

# plot traces -------------------------------------------------------------

#' Plot sample traces
#'
#' Plot the raw trace data
#'
#' @param fragments_list A list of fragments_repeats or fragments_trace objects containing fragment data.
#' @param show_peaks If peak data are available, TRUE will plot the peaks on top of the trace as dots.
#' @param n_facet_col A numeric value indicating the number of columns for faceting in the plot.
#' @param sample_subset A character vector of unique ids for a subset of samples to plot
#' @param xlim the x limits of the plot. A numeric vector of length two.
#' @param ylim the y limits of the plot. A numeric vector of length two.
#' @param x_axis A character indicating what should be plotted on the x-axis, chose between `size` or `repeats`. If neither is selected, an assumption is made based on if repeats have been called.
#' @param height_color_threshold Threshold relative to tallest peak to color the dots (blue above, purple below).
#'
#' @return plot traces from fragments object
#' @export
#'
#' @importFrom methods is
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices rgb
#' @importFrom graphics plot.new
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom graphics barplot
#'
#'
#' @details
#' A plot of the raw signal by bp size. Red vertical line indicates the scan was
#' flagged as off-scale. This is in any channel, so use your best judgment to determine
#' if it's from the sample or ladder channel.
#'
#' If peaks are called, green is the tallest peak, blue is peaks above the height threshold (default 5%), purple is below the height threshold. If `force_whole_repeat_units` is used within [call_repeats()], the called repeat will be connected to the peak in the trace with a horizontal dashed line.
#'
#' The index peak will be plotted as a vertical dashed line when it has been set using `assign_index_peaks()`.
#'
#'
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
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
#'
#' # Simple conversion from bp size to repeat size
#' call_repeats(
#'   fragments_list
#' )
#'
#' plot_traces(fragments_list, xlim = c(105, 150))
#'
plot_traces <- function(
    fragments_list,
    show_peaks = TRUE,
    n_facet_col = 1,
    sample_subset = NULL,
    xlim = NULL,
    ylim = NULL,
    x_axis = NULL,
    height_color_threshold = 0.05) {
  if (!is.null(sample_subset)) {
    fragments_list <- fragments_list[which(names(fragments_list) %in% sample_subset)]
  }

  graphics::par(mfrow = c(ceiling(length(fragments_list) / n_facet_col), n_facet_col)) # Adjust layout as needed
  for (i in seq_along(fragments_list)) {
    fragments_list[[i]]$plot_trace(
      show_peaks = show_peaks,
      xlim = xlim,
      ylim = ylim,
      x_axis = x_axis,
      height_color_threshold = height_color_threshold
    )
  }
  graphics::par(mfrow = c(1, 1)) # Reset the layout
}

# plot fragment data -------------------------------------------------------

#' Plot Peak Data
#'
#' Plots peak data from a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#' @param n_facet_col A numeric value indicating the number of columns for faceting in the plot.
#' @param sample_subset A character vector of unique ids for a subset of samples to plot
#' @param xlim the x limits of the plot. A numeric vector of length two.
#' @param ylim the y limits of the plot. A numeric vector of length two.
#'
#' @return A plot object displaying the peak data.
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#'
#' fragments_list <- peak_table_to_fragments(gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B",
#'   min_size_bp = 300
#' )
#'
#' find_alleles(
#'   fragments_list
#' )
#'
#' plot_fragments(fragments_list[1])
plot_fragments <- function(
    fragments_list,
    n_facet_col = 1,
    sample_subset = NULL,
    xlim = NULL,
    ylim = NULL) {
  if (!is.null(sample_subset)) {
    fragments_list <- fragments_list[which(names(fragments_list) %in% sample_subset)]
  }

  graphics::par(mfrow = c(ceiling(length(fragments_list) / n_facet_col), n_facet_col)) # Adjust layout as needed
  for (i in seq_along(fragments_list)) {
    fragments_list[[i]]$plot_fragments(
      xlim = xlim,
      ylim = ylim
    )
  }
  graphics::par(mfrow = c(1, 1)) # Reset the layout
}


# plot batch correction sample groups ----------------------------------------

#' Plot correction samples
#'
#' Plot the overlapping traces of the batch control samples
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data. must have trace information.
#' @param selected_sample A character vector of batch_sample_id for a subset of samples to plot. Or alternatively supply a number to select batch sample by position in alphabetical order.
#' @param xlim the x limits of the plot. A numeric vector of length two.
#'
#' @return plot of batch corrected samples
#' @export
#' @importFrom grDevices recordPlot replayPlot
#' 
#' @details
#' A plot of the raw signal by bp size or repeats for the batch correction samples. 
#'
#' When plotting the traces before repeat correction, we do not expect the samples to be closely overlapping due to run-to-run variation. After repeat correction, the traces should be basically overlapping. 
#' 
#' These plots are made using base R plotting. Sometimes these fail to render in the viewing panes of IDEs (eg you get the error 'Error in `plot.new()`: figure margins too large)'. If this happens, try saving the plot as a pdf using traditional approaches (see grDevices::pdf). 
#'
#' @seealso [call_repeats()] for more info on batch correction.
#' @examples
#'
#' fsa_list <- lapply(cell_line_fsa_list[91:94], function(x) x$clone())
#'
#' find_ladders(fsa_list, show_progress_bar = FALSE)
#'
#' fragments_list <- find_fragments(fsa_list, min_bp_size = 300)
#'
#' test_alleles <- find_alleles(
#'   fragments_list 
#' )
#' 
#' add_metadata(
#'   fragments_list,
#'   metadata
#' )
#'
#'
#' call_repeats(
#'   fragments_list = fragments_list,
#'   repeat_calling_algorithm = "none",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3,
#'   batch_correction = TRUE
#' )
#'
#' # traces of bp size shows traces at different sizes
#' plot_batch_correction_samples(
#'   fragments_list,
#'   selected_sample = "S-21-212", xlim = c(100, 120)
#' )
#'
#'
plot_batch_correction_samples <- function(
  fragments_list,
  selected_sample,
  xlim = NULL) {
# first check to see if repeats have been called or batches corrected and give some warnings
  
if(all(sapply(fragments_list, function(x) is.null(x$repeat_table_df)))){
  warning(call. = FALSE,
    "Repeats not detected. Only samples before correction values will be plotted. Use call_repeats(correction = 'batch') or call_repeats(correction = 'repeat') to see the effect of batch correction on selected sample."
  )
} else if(!any(sapply(fragments_list, function(x) !is.na(x$.__enclos_env__$private$batch_correction_factor))) & !any(sapply(fragments_list, function(x) !is.na(x$.__enclos_env__$private$repeat_correction_factor)))){
  warning(call. = FALSE,
    "Batch or repeat correction not detected. Only samples before correction values will be plotted. Use call_repeats(correction = 'batch') or call_repeats(correction = 'repeat') to see the effect of batch correction on selected sample."
  )
}
  
size_standard_fragments <- sapply(fragments_list, function(x) x$batch_sample_id)
controls_fragments_list <- fragments_list[which(!is.na(size_standard_fragments))]

if (length(unique(na.omit(size_standard_fragments))) == 0) {
  stop(
    call. = FALSE,
    "There are no samples with batch_sample_id assigned. Check that the batch_sample_id has been added to the samples via add_metadata()."
  )
}

if (!is.null(selected_sample) & is.character(selected_sample)) {
  selected_sample_list <- sapply(controls_fragments_list, function(x) x$batch_sample_id %in% selected_sample)
  controls_fragments_list <- controls_fragments_list[which(selected_sample_list)]

  if (length(controls_fragments_list) == 0) {
    stop(
      call. = FALSE,
      "After subsetting the samples with the provided id, no samples were left. Check that you provided the correct id or that the batch_sample_id has been added to the samples."
    )
  }
} else if(is.numeric(selected_sample)){
  batch_sample_ids <- sapply(controls_fragments_list, function(x) x$batch_sample_id)
  unique_batch_sample_ids <- unique(na.omit(batch_sample_ids))
  if(selected_sample > length(unique_batch_sample_ids)){
    stop(
      call. = FALSE,
      paste0("Invalid selected_sample number, there are only ", length(unique_batch_sample_ids), "unique batch_sample_ids")
    )
  }
  controls_fragments_list <- controls_fragments_list[which(batch_sample_ids == unique_batch_sample_ids[selected_sample])]
} 

overlapping_plot <- function(sample_fragments, before_or_after_correction) {


  sample_traces <- vector("list", length(sample_fragments))
  for (i in seq_along(sample_fragments)) {
    sample_traces[[i]] <- sample_fragments[[i]]$trace_bp_df

    # check to see if repeats have been called to see what should be plotted on the x_axis
    if(all(sapply(sample_fragments, function(x) is.null(x$repeat_table_df)))){
      x_axis <- "size"
      sample_traces[[i]]$x <- sample_traces[[i]]$size
    } else{
      x_axis <- "repeats"
      if(before_or_after_correction == "before"){
        # re-calculate repeats using uncorrected size
        sample_traces[[i]]$x <- (sample_traces[[i]]$size - sample_fragments[[i]]$.__enclos_env__$private$assay_size_without_repeat) / sample_fragments[[i]]$.__enclos_env__$private$repeat_size
      } else{
        sample_traces[[i]]$x <- sample_traces[[i]]$calculated_repeats
      }
    }
  }
  # Generate colors dynamically
  n_dfs <- length(sample_traces)
  colors <- rainbow(n_dfs, alpha = 0.5) # Generates n colors with alpha for transparency


  if (!is.null(xlim)) {
    sample_traces <- lapply(sample_traces, function(df) {
      df <- df[which(df$x < xlim[2] & df$x > xlim[1]), ]
      return(df)
    })
  }

  # normalize signal to samples have the same maximum
  sample_traces <- lapply(sample_traces, function(x){
    x$signal <- x$signal - min(x$signal)
    x$rel_signal <- x$signal / max(x$signal)
    return(x)
  })

  plot_title <- ifelse(before_or_after_correction == "before", "Before correction", "After correction")
  plot_title <- paste0(plot_title, paste0("\n", unique(sapply(sample_fragments, function(x) x$batch_sample_id))))

  # Plot the first dataframe
  plot(sample_traces[[1]]$x, sample_traces[[1]]$rel_signal,
    type = "l",
    col = colors[1],
    xlab = ifelse(x_axis == "size", "Size", "Repeats"),
    ylab = "Signal",
    ylim = range(sapply(sample_traces, function(df) range(df$rel_signal))),
    main = plot_title
  )

  # Add lines for remaining dataframes
  for (i in 2:n_dfs) {
    graphics::lines(sample_traces[[i]]$x, sample_traces[[i]]$rel_signal, col = colors[i])

  }

  graphics::legend("topright", 
  legend = sapply(sample_fragments, function(x) x$unique_id), 
  col = colors,
  lty = 1, lwd = 1,
  cex = 0.8, # reduce the size of the text
  inset = c(-0.5, 0),  # Adjust position to be inside the plot, but in the corner
  x.intersp = 0.25,   # Adjust horizontal spacing between symbols and text
  y.intersp = 0.6,   # Adjust vertical spacing to align text under the lines
  bty = "n"  # Remove the box around the legend
)       

}
graphics::par(mfrow = c(1, 2)) 
# for some reason we need to use the recordPlot() strategy below.
# just looping over the plots only rendered one for some reason
recorded_plots <- vector("list", 2)

#before correction plot
overlapping_plot(controls_fragments_list, before_or_after_correction = "before")
recorded_plots[[1]] <- grDevices::recordPlot()
  
#after correction plot
if(!all(sapply(controls_fragments_list, function(x) is.null(x$repeat_table_df))) && any(sapply(controls_fragments_list, function(x) !is.na(x$.__enclos_env__$private$batch_correction_factor))) | any(sapply(controls_fragments_list, function(x) !is.na(x$.__enclos_env__$private$repeat_correction_factor)))){
  overlapping_plot(controls_fragments_list, before_or_after_correction = "after")
  recorded_plots[[2]] <- grDevices::recordPlot()
} else{
  recorded_plots[[2]] <- plot.new()
}

for (i in seq_along(recorded_plots)) {
  grDevices::replayPlot(recorded_plots[[i]])
}

graphics::par(mfrow = c(1, 1)) # Reset the layout
  
}



#' plot_data_channels
#'
#' Plot the raw data from the fsa file
#'
#' @param fragments_list A list of fragments_trace objects.
#' @param n_facet_col A numeric value indicating the number of columns for faceting in the plot.
#' @param sample_subset A character vector of unique ids for a subset of samples to plot
#'
#' @return a plot of the raw data channels
#' @export
#' @details
#' A plot of the raw data channels in the fsa file.
#' 
#' These plots are made using base R plotting. Sometimes these fail to render in the viewing panes of IDEs (eg you get the error 'Error in `plot.new()`: figure margins too large)'. If this happens, try saving the plot as a pdf using traditional approaches (see grDevices::pdf). To get it to render in the IDE pane, trying matching `n_facet_col` to the number of samples you're attempting to plot, or using `sample_subset` to limit it to a single sample.
#'
#' @examples
#' 
#' plot_data_channels(cell_line_fsa_list[1])
#'
plot_data_channels <- function(
    fragments_list,
    sample_subset = NULL,
    n_facet_col = 1) {
  if (!is.null(sample_subset)) {
    fragments_list <- fragments_list[which(names(fragments_list) %in% sample_subset)]
  }
  graphics::par(mfrow = c(ceiling(length(fragments_list) / n_facet_col), n_facet_col)) # Adjust layout as needed
  # for some reason we need to use the recordPlot() strategy below.
  # just looping over the plots only rendered one for some reason
  recorded_plots <- vector("list", length(fragments_list))
  for (i in seq_along(fragments_list)) {
    fragments_list[[i]]$plot_data_channels()
    # Record the plot
    recorded_plots[[i]] <- grDevices::recordPlot()
  }

  for (i in seq_along(recorded_plots)) {
    grDevices::replayPlot(recorded_plots[[i]])
  }

  graphics::par(mfrow = c(1, 1)) # Reset the layout
}



#' Plot Repeat Correction Model
#'
#' Plots the results of the repeat correction model for a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats class objects obtained from the 'call_repeats' function when the 'repeat_length_correction' was either 'from_metadata' or 'from_genemapper'.
#'
#' @return A base R graphic object displaying the repeat correction model results.
#' @export
#'
#' @examples
#'
#'
#' gm_raw <- instability::example_data
#' metadata <- instability::metadata
#'
#' test_fragments <- peak_table_to_fragments(
#'   gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B"
#' )
#'
#' test_alleles <- find_alleles(
#'   fragments_list = test_fragments,
#'   number_of_peaks_to_return = 2,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1
#' )
#'
#' test_metadata <- add_metadata(
#'   fragments_list = test_alleles,
#'   metadata_data.frame = metadata
#' )
#'
#' test_repeats_corrected <- call_repeats(
#'   fragments_list = test_metadata,
#'   repeat_calling_algorithm = "none",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3,
#'   repeat_length_correction = "from_metadata"
#' )
#'
#' plot_repeat_correction_model(test_repeats_corrected)
#'
plot_repeat_correction_model <- function(fragments_list, batch_run_id_subset = NULL) {
  # Check if all models in the list are the same
  first_model_df <- fragments_list[[1]]$.__enclos_env__$private$repeat_correction_mod
  identical_model_test <- logical(length(fragments_list))
  for (i in seq_along(fragments_list)) {
    identical_model_test[i] <- identical(first_model_df, fragments_list[[i]]$.__enclos_env__$private$repeat_correction_mod)
  }

  if (!all(identical_model_test)) {
    stop("The supplied fragments list must come from the same 'call_repeats' function output", call. = FALSE)
  }

  controls_repeats_df <- fragments_list[[1]]$.__enclos_env__$private$repeat_correction_mod$model
  controls_repeats_df$unique_id <- sub("\\.[0-9]+$", "", row.names(controls_repeats_df))
  # add back in batch_run_id if it's not there (because a different lm is made when just one run)
  # assume that all the samples are the same batch since they have identical model
  if(!"batch_run_id" %in% names(controls_repeats_df)){
    controls_repeats_df$batch_run_id <- rep(fragments_list[[1]]$batch_run_id, nrow(controls_repeats_df))
  }

  # Plotting
  unique_batch_run_ids <- unique(controls_repeats_df$batch_run_id)

  if(!is.null(batch_run_id_subset) && is.numeric(batch_run_id_subset)){
    if(batch_run_id_subset > length(unique_batch_run_ids)){
      stop(call. = FALSE, paste0("The 'batch_run_id_subset' number was too large. There are only ",length(unique_batch_run_ids), " 'batch_run_id'."))
    }
    unique_batch_run_ids <- unique_batch_run_ids[batch_run_id_subset]
  } else if(is.character(batch_run_id_subset)){
    unique_batch_run_ids <- unique_batch_run_ids[which(unique_batch_run_ids %in% batch_run_id_subset)]
  }
# TOODO
  # add n_facet_col

  graphics::par(mfrow = c(1, length(unique_batch_run_ids)))
  recorded_plots <- vector("list", length(unique_batch_run_ids))
  for (i in 1:length(unique_batch_run_ids)) {
    plate_data <- controls_repeats_df[which(controls_repeats_df$batch_run_id == unique_batch_run_ids[i]),]

    # Generating unique colors for each unique value in unique_id
    unique_ids <- unique(plate_data$unique_id)
    colors <- grDevices::rainbow(length(unique_ids))
    id_color_map <- setNames(colors, unique_ids)

    plot(plate_data$size, plate_data$validated_repeats,
      pch = 21, col = id_color_map[plate_data$unique_id],
      cex = 2, main = paste("Plate ID:", unique_batch_run_ids[i]), xlab = "Size", ylab = "User supplied repeat length"
    )

    lm_model <- lm(validated_repeats ~ size, data = plate_data)
    graphics::abline(lm_model, col = "blue")

    # Record the plot
    recorded_plots[[i]] <- grDevices::recordPlot()
  }

  for (i in 1:length(unique_batch_run_ids)) {
    grDevices::replayPlot(recorded_plots[[i]])
  }

}

