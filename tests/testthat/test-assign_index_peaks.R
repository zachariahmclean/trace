test_that("index assignment", {

  suppressWarnings(
    test_ladders <- find_ladders(cell_line_fsa_list,
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


  peak_list <- find_fragments(test_ladders,
                              minimum_peak_signal = 20,
                              min_bp_size = 300
  )

  fragment_metadata <- add_metadata(
    fragments_list = peak_list,
    metadata_data.frame = metadata
  )

  fragment_alleles <- find_alleles(
    fragments_list = fragment_metadata,
    number_of_peaks_to_return = 1
  )

  suppressMessages(
    suppressWarnings(
      test_repeats <- call_repeats(
        fragments_list = fragment_alleles
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(100, 150),
  #             ylim = c(0,2000))




  suppressMessages(
    suppressWarnings(
      test_assignment <- assign_index_peaks(
        test_repeats,
        grouped = TRUE
      )
    )
  )
  # plot_traces(test_assignment[1:9], n_facet_col = 3,
  #             xlim = c(100, 150),
  #             ylim = c(0,2000))



  # plot_fragments(test_repeats[1:4])
  # plot_size_standard_model(test_repeats)

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_assignment,
        peak_threshold = 0.05,
        window_around_index_peak = c(-40, 40)
      )
    )
  )

testthat::expect_true(all(sapply(test_assignment, function(x) x$.__enclos_env__$private$assigned_index_peak_used)))




})
