testthat::test_that("find_fragments", {
  # simple option:

  file_list <- trace::cell_line_fsa_list
  suppressWarnings(
    test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )

  suppressWarnings(
    peak_list <- find_fragments(test_ladders,
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



# metadata transfer



testthat::test_that("metadata transfer", {
  file_list <- trace::cell_line_fsa_list


  suppressWarnings(
    test_ladders <- find_ladders(file_list[1],
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )
  suppressWarnings(
    metadata_added <- add_metadata(test_ladders, metadata
    )
  )

  peak_list <- find_fragments(metadata_added,
    minimum_peak_signal = 20,
    min_bp_size = 300
  )

  testthat::expect_true(metadata_added[[1]]$unique_id == peak_list[[1]]$unique_id)
  testthat::expect_true(metadata_added[[1]]$batch_run_id == peak_list[[1]]$batch_run_id)
  testthat::expect_true(metadata_added[[1]]$metrics_group_id == peak_list[[1]]$metrics_group_id)
  testthat::expect_true(metadata_added[[1]]$metrics_baseline_control == peak_list[[1]]$metrics_baseline_control)
})
