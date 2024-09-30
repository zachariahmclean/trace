testthat::test_that("find_fragments", {
  # none option:

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  fsa_list <- fsa_list["20230413_B03.fsa"]
  suppressWarnings(
    find_ladders(fsa_list,
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )

  suppressWarnings(
    peak_list <- find_fragments(fsa_list,
      minimum_peak_signal = 20,
      min_bp_size = 100
    )
  )



  # plot_traces(peak_list[which(names(peak_list) =="20230413_B03.fsa")],
  #             xlim = c(450,650),
  #             ylim = c(0, 1000))


  extracted_fragments <- extract_fragments(peak_list[which(names(peak_list) == "20230413_B03.fsa")])


  tall_peaks <- extracted_fragments[which(extracted_fragments$size > 500 & extracted_fragments$height > 700), ]


  testthat::expect_true(all(round(tall_peaks$size, 4) == c(547.5879, 550.0142, 552.4404, 554.8666, 557.2929)))
})


