testthat::test_that("find_fragments", {
  config <- load_config()

  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())

  suppressWarnings(
    find_ladders(fsa_list,
      config,
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )

  suppressWarnings(
    find_fragments(fsa_list,
      config,
      minimum_peak_signal = 20,
      min_bp_size = 100
    )
  )



  # plot_traces(fsa_list[which(names(fsa_list) =="20230413_B03.fsa")],
  #             xlim = c(450,650),
  #             ylim = c(0, 1000))


  extracted_fragments <- extract_fragments(fsa_list)


  tall_peaks <- extracted_fragments[which(extracted_fragments$size > 500 & extracted_fragments$signal > 700), ]


  testthat::expect_true(all(round(tall_peaks$size, 4) == c(500.7141, 503.5075, 506.1264)))
})


