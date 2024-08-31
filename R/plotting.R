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


  allele_1_mode <- ifelse(is.null(fragment_repeats$repeat_table_df), round(fragment_repeats$get_alleles()$allele_1_size), round(fragment_repeats$get_alleles()$allele_1_repeat))
  allele_2_mode <- ifelse(is.null(fragment_repeats$repeat_table_df), round(fragment_repeats$get_alleles()$allele_2_size), round(fragment_repeats$get_alleles()$allele_2_repeat))

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
    col = sapply(all_x_values, function(x) if (!is.na(allele_1_mode) && x == allele_1_mode) "red" else if (!is.na(allele_2_mode) && x == allele_2_mode) "blue" else "gray")
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
    if (!is.null(fragments$get_alleles()$allele_1_height) && !is.na(fragments$get_alleles()$allele_1_height)) {
      tallest_peak_height <- fragments$get_alleles()$allele_1_height
      # find the tallest peak x axis position
      if (is.null(x_axis) && is.na(fragments$get_alleles()$allele_1_repeat)) {
        tallest_peak_x <- fragments$get_alleles()$allele_1_size
      } else if (is.null(x_axis) && !is.na(fragments$get_alleles()$allele_1_repeat)) {
        tallest_peak_x <- fragments$get_alleles()$allele_1_repeat
      } else if (x_axis == "size") {
        tallest_peak_x <- fragments$get_alleles()$allele_1_size
      } else {
        tallest_peak_x <- fragments$get_alleles()$allele_1_repeat
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
#' file_list <- trace::cell_line_fsa_list
#'
#' test_ladders <- find_ladders(file_list)
#'
#' # Manually inspect the ladders
#' plot_ladders(test_ladders[1], n_facet_col = 1)
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
#' file_list <- trace::cell_line_fsa_list[1]
#'
#' test_ladders <- find_ladders(file_list)
#'
#' fragments_list <- find_fragments(test_ladders,
#'   min_bp_size = 300
#' )
#'
#' test_alleles <- find_alleles(
#'   fragments_list = fragments_list
#' )
#'
#' # Simple conversion from bp size to repeat size
#' test_repeats <- call_repeats(
#'   fragments_list = test_alleles
#' )
#'
#' plot_traces(test_repeats[1], xlim = c(105, 150))
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
#' @return A base R plot object displaying the peak data.
#' @export
#'
#' @examples
#' gm_raw <- trace::example_data
#'
#' test_fragments <- peak_table_to_fragments(gm_raw,
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
#' plot_fragments(test_alleles[1:2])
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


# plot size standard sample groups ----------------------------------------

#' plot size standard samples
#'
#' Plot the overlapping traces of the size standards by their size standard ids
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#' @param sample_subset A character vector of batch_sample_id for a subset of samples to plot. Or alternativly supply a numeric vector.
#' @param xlim the x limits of the plot. A numeric vector of length two.
#' @param x_axis A character indicating what should be plotted on the x-axis, chose between `size` or `repeats`. Only use repeats if plotting after the repeat correction.
#' @param n_facet_col A numeric value indicating the number of columns for faceting in the plot.
#'
#' @return plot traces from fragments object
#' @export
#' @importFrom grDevices recordPlot replayPlot
#' 
#' @details
#' A plot of the raw signal by bp size or repeats for the size standard samples. The cicle at the top of the plot is for the called allele for that sample.
#'
#' When plotting the traces before repeat correction, we do not expect the samples to be closely overallping due to run-to-run variation. After repeat correction and plotting "repeats" on the x-axis, the traces should be bascially overlapping. It can be difficult from the "repeats" x-axis to figure out which sample is wrong because if one is wrong it will mess up the repeat size for all other samples in that same batch_run_id. Use the "size" x-axis to make sure all of the traces have the same distribution and modal peak."
#' 
#' These plots are made using base R plotting. Sometimes these fail to render in the viewing panes of IDEs (eg you get the error 'Error in `plot.new()`: figure margins too large)'. If this happens, try saving the plot as a pdf using traditional approaches (see grDevices::pdf). To get it to render in the IDE pane, trying matching `n_facet_col` to the number of samples you're attmpting to plot, or using `sample_subset` to limit it to a single sample.
#'
#' @examples
#'
#'
#' test_ladders <- find_ladders(cell_line_fsa_list, show_progress_bar = FALSE)
#'
#'
#'
#' # duplicate data to generate an example
#'
#' test_metadata <- add_metadata(
#'   fragments_list = test_ladders,
#'   metadata_data.frame = metadata
#' )
#'
#' test_metadata2 <- lapply(test_metadata, function(sample) {
#'   sample2 <- sample$clone()
#'   sample2$trace_bp_df$size <- sample2$trace_bp_df$size + 2
#'   sample2$unique_id <- paste0(sample$unique_id, "_2")
#'   sample2$batch_run_id <- paste0(sample2$batch_run_id, "_2")
#'   return(sample2)
#' })
#'
#' names(test_metadata2) <- sapply(test_metadata2, function(x) x$unique_id)
#'
#' metadata_added_combined <- c(test_metadata, test_metadata2)
#'
#' test_fragments <- find_fragments(metadata_added_combined, min_bp_size = 300)
#'
#' test_alleles <- find_alleles(
#'   fragments_list = test_fragments,
#'   number_of_peaks_to_return = 1,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1
#' )
#'
#'
#' test_repeats_corrected <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3,
#'   batch_correction = TRUE
#' )
#'
#' # traces of bp size shows traces at different sizes
#' plot_batch_correction_samples(test_repeats_corrected,
#'   x_axis = "size",
#'   sample_subset = "S-21-212", xlim = c(400, 450)
#' )
#'
#' # overlapping traces when looking at the corrected repeat length
#' plot_batch_correction_samples(test_repeats_corrected,
#'   x_axis = "repeats",
#'   sample_subset = "S-21-212", xlim = c(100, 130)
#' )
#'
plot_batch_correction_samples <- function(
  fragments_list,
  sample_subset = NULL,
  x_axis = "size",
  n_facet_col  = 1, 
  xlim = NULL) {
size_standard_fragments <- sapply(fragments_list, function(x) x$batch_sample_id)
controls_fragments_list <- fragments_list[which(!is.na(size_standard_fragments))]

if (length(unique(na.omit(size_standard_fragments))) == 0) {
  stop(
    call. = FALSE,
    "There are no samples with batch_sample_id assigned. Check that the batch_sample_id has been added to the samples via add_metadata()."
  )
}

if (!is.null(sample_subset)) {
  sample_subset <- sapply(controls_fragments_list, function(x) x$batch_sample_id %in% sample_subset)
  controls_fragments_list <- controls_fragments_list[which(sample_subset)]

  if (length(controls_fragments_list) == 0) {
    stop(
      call. = FALSE,
      "After subsetting the samples with the provided id, no samples were left. Check that you provided the correct id or that the batch_sample_id has been added to the samples."
    )
  }
}

size_standard_fragments_sample_groups <- sapply(controls_fragments_list, function(x) x$batch_sample_id)

split_by_sample <- split(controls_fragments_list, size_standard_fragments_sample_groups)

overlapping_plot <- function(sample_fragments) {
  sample_traces <- lapply(sample_fragments, function(y) {
    df <- y$trace_bp_df
    if (x_axis == "size") {
      df$x <- df$size
    } else if (x_axis == "repeats") {
      df$x <- df$calculated_repeats
    } else {
      stop(
        call. = FALSE,
        "Please provide valid x-axis, either: 'size' or 'repeats'"
      )
    }
    return(df)
  })

  # Generate colors dynamically
  n_dfs <- length(sample_traces)
  colors <- rainbow(n_dfs, alpha = 0.5) # Generates n colors with alpha for transparency


  if (!is.null(xlim)) {
    sample_traces <- lapply(sample_traces, function(df) {
      df <- df[which(df$x < xlim[2] & df$x > xlim[1]), ]
      return(df)
    })
  }

  # normalize signal to samples have the same maxium
  sample_traces <- lapply(sample_traces, function(x){
    x$signal <- x$signal - min(x$signal)
    x$rel_signal <- x$signal / max(x$signal)
    return(x)
  })


  # Plot the first dataframe
  plot(sample_traces[[1]]$x, sample_traces[[1]]$rel_signal,
    type = "l",
    col = colors[1],
    xlab = ifelse(x_axis == "size", "Size", "Repeats"),
    ylab = "Signal",
    ylim = range(sapply(sample_traces, function(df) range(df$rel_signal))),
    main = sample_fragments[[1]]$batch_sample_id
  )
  # also add point for tallest peak. sample_traces and sample_fragments are in the same order
  points(
    ifelse(x_axis == "size", sample_fragments[[1]]$get_alleles()$allele_1_size, sample_fragments[[1]]$get_alleles()$allele_1_repeat),
    1,
    col = colors[1]
  )

  # Add lines for remaining dataframes
  for (i in 2:n_dfs) {
    graphics::lines(sample_traces[[i]]$x, sample_traces[[i]]$rel_signal, col = colors[i])
    points(
      ifelse(x_axis == "size", sample_fragments[[i]]$get_alleles()$allele_1_size, sample_fragments[[i]]$get_alleles()$allele_1_repeat),
      1,
      col = colors[i]
    )
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
  
graphics::par(mfrow = c(ceiling(length(split_by_sample) / n_facet_col), n_facet_col)) # Adjust layout as needed
# for some reason we need to use the recordPlot() strategy below.
# just looping over the plots only rendered one for some reason
recorded_plots <- vector("list", length(split_by_sample))
for (i in seq_along(split_by_sample)) {
  overlapping_plot(split_by_sample[[i]])
  # Record the plot
  recorded_plots[[i]] <- grDevices::recordPlot()
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
#' @return a plot of ladders
#' @export
#' @details
#' A plot of the raw data channels in the fsa file.
#' 
#' These plots are made using base R plotting. Sometimes these fail to render in the viewing panes of IDEs (eg you get the error 'Error in `plot.new()`: figure margins too large)'. If this happens, try saving the plot as a pdf using traditional approaches (see grDevices::pdf). To get it to render in the IDE pane, trying matching `n_facet_col` to the number of samples you're attmpting to plot, or using `sample_subset` to limit it to a single sample.
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


