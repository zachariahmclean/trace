



# plot fragments ----------------------------------------------------------


testthat::test_that("full pipeline", {

  config <- load_config()

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  suppressWarnings(
    find_ladders(fsa_list,
      config,
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 5,
      show_progress_bar = FALSE
    )
  )

  # plot_ladders(test_ladders[1:9], n_facet_col = 3,
  #              xlim = c(1000, 4800),
  #              ylim = c(0, 15000))


  # # Start a PDF device
  # pdf(file = "C:/Users/zlm2/Downloads/ladder.pdf", width = 12, height = 6) # Set width and height as desired
  #
  # # Loop through the list of plots
  # for (i in seq_along(test_ladders)) {
  #   test_ladders[[i]]$plot_ladder(xlim = c(1400, 4500))
  # }
  #
  # # Close the PDF device
  # dev.off()


  find_fragments(fsa_list,
    config,
    minimum_peak_signal = 20,
    min_bp_size = 300
  )

add_metadata(
  fsa_list,
    metadata_data.frame = metadata
  )

find_alleles(
  fsa_list,
  config
  )

  suppressMessages(
    suppressWarnings(
      call_repeats(
        fsa_list,
        config
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(400, 550),
  #             ylim = c(0,2000))


  # plot_fragments(test_repeats[1:4])


  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        fsa_list,
        config,
        grouped = TRUE
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fsa_list,
        peak_threshold = 0.05,
        window_around_index_peak = c(-40, 40)
      )
    )
  )


  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_signal > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$metrics_group_id)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_change / median(plot_data[[i]]$average_repeat_change[which(plot_data[[i]]$treatment == 0)])
  }

  plot_data <- do.call(rbind, plot_data)

  # Revise genotype levels

  plot_data$genotype <- factor(plot_data$genotype, levels = c("non-edited", "edited"))



  # ggplot2::ggplot(plot_data,
  #                 ggplot2::aes(as.factor(treatment), rel_gain,
  #            colour = as.factor(treatment))) +
  #   ggplot2::geom_boxplot(outlier.shape = NA) +
  #   ggplot2::geom_jitter() +
  #   ggplot2::facet_wrap(ggplot2::vars(genotype)) +
  #   ggplot2::labs(y = "Average repeat gain\n(relative to DMSO)",
  #        x = "Branaplam (nM)") +
  #   ggplot2::theme(legend.position = "none")


  medians <- aggregate(rel_gain ~ treatment + genotype, plot_data, median, na.rm = TRUE)

  expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.86154, 0.73268, 0.55720)))
})

