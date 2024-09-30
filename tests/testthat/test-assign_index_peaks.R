test_that("index assignment", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  suppressWarnings(
    test_fragments <- peak_table_to_fragments(gm_raw,
      data_format = "genemapper5",
      dye_channel = "B",
      min_size_bp = 400
    )
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments
  )

  suppressMessages(
    suppressWarnings(
      call_repeats(
        fragments_list = test_fragments
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

testthat::expect_true(all(sapply(test_fragments, function(x) x$.__enclos_env__$private$assigned_index_peak_used)))




})



# test warning that index was assigned without batch correction



testthat::test_that("calculate metrics", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  suppressWarnings(
    test_fragments <- peak_table_to_fragments(gm_raw,
      data_format = "genemapper5",
      dye_channel = "B",
      min_size_bp = 400
    )
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      repeat_calling_algorithm = "none",
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )


  #purposely messs up one batch run id to generate warning

  test_fragments[[1]]$batch_run_id <- "wrong_run"

suppressMessages(
  tryCatch({
    assign_index_peaks(
      test_fragments,
      grouped = TRUE
    )
  },
    warning = function(w){
      assignment_warning <<- w
    }
  )
)

testthat::expect_true(class(assignment_warning)[1] == "simpleWarning")
testthat::expect_true(grepl("20230413_A01.fsa", assignment_warning))
testthat::expect_true(grepl("batch_run_id", assignment_warning))

})



testthat::test_that("test situation where some samples have NA in grouped", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  suppressWarnings(
    test_fragments <- peak_table_to_fragments(gm_raw,
      data_format = "genemapper5",
      dye_channel = "B",
      min_size_bp = 400
    )
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      repeat_calling_algorithm = "none",
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )


  #purposely messs up one batch run id to generate warning

  test_fragments[[1]]$metrics_group_id <- NA_character_

suppressMessages(
  tryCatch({
    assign_index_peaks(
      test_fragments,
      grouped = TRUE
    )
  },
    warning = function(w){
      assignment_warning <<- w
    }
  )
)

testthat::expect_true(class(assignment_warning)[1] == "simpleWarning")
testthat::expect_true(grepl("Group 'NA' has no 'metrics_baseline_control'", assignment_warning))

})




