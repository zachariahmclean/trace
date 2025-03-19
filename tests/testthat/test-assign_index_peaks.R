test_that("index assignment", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  config <- load_config()
  # Save raw data as a fragment class

  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(gm_raw,
      dye_channel = "B",
      min_size_bp = 400
    )
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments,
    config
  )

  suppressMessages(
    suppressWarnings(
      call_repeats(
        fragments_list = test_fragments,
        config
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(100, 150),
  #             ylim = c(0,2000))




  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        test_fragments,
        config,
        grouped = TRUE
      )
    )
  )
  # plot_traces(test_assignment[1:9], n_facet_col = 3,
  #             xlim = c(100, 150),
  #             ylim = c(0,2000))



  # plot_fragments(test_repeats[1:4])


  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_fragments,
        peak_threshold = 0.05,
        window_around_index_peak = c(-40, 40)
      )
    )
  )

  # come up with new test
# testthat::expect_true(all(sapply(test_fragments, function(x) x$.__enclos_env__$private$assigned_index_peak_used)))




})



# test warning that index was assigned without batch correction



testthat::test_that("calculate metrics", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  config <- load_config()
  # Save raw data as a fragment class

  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(gm_raw,
      dye_channel = "B",
      min_size_bp = 400
    )
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments,
    config
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      config,
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )


  #purposely messs up one batch run id to generate warning

  test_fragments[[1]]$batch_run_id <- "wrong_run"

suppressMessages(
    assign_index_peaks_output <- assign_index_peaks(
      test_fragments,
      config,
      grouped = TRUE
    )
)

testthat::expect_true(assign_index_peaks_output$status == "warning")
testthat::expect_true(grepl("20230413_A07.fsa", assign_index_peaks_output$warning_message))
testthat::expect_true(grepl("batch_run_id", assign_index_peaks_output$warning_message))

})



testthat::test_that("test situation where some samples have NA in grouped", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  config <- load_config()
  # Save raw data as a fragment class

  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(gm_raw,
      dye_channel = "B",
      min_size_bp = 400
    )
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments,
    config
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      config,
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )


  #purposely messs up one batch run id to generate warning

  test_fragments[[1]]$metrics_group_id <- NA_character_

suppressMessages(
    assign_index_peaks_output <- assign_index_peaks(
      test_fragments,
      config,
      grouped = TRUE
    )
)

testthat::expect_true(assign_index_peaks_output$status == "warning")
testthat::expect_true(grepl("The following 'metrics_group_id' have no corresponding 'metrics_baseline_control': NA", assign_index_peaks_output$warning_message[2]))

})




