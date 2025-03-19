# percentiles

testthat::test_that("percentiles", {
  gm_raw <- trace::example_data
  test_sample <- unique(gm_raw$Sample.File.Name)[1]
  config <- load_config()

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- genemapper_table_to_fragments(test_df,
    dye_channel = "B",
    min_size_bp = 350
  )

  find_alleles(test_fragments[1], config)

  call_repeats(
    test_fragments,
    config,
    assay_size_without_repeat = 87,
    repeat_size = 3
  )


  test_distribution_df <- test_fragments[[1]]$repeat_table_df
  test_distribution_df <- test_distribution_df[which(test_distribution_df$repeats > 113), ]


  percentiles <- find_percentiles(
    repeats = test_distribution_df$repeats,
    signals = test_distribution_df$signal,
    index_peak_repeat = test_fragments[[1]]$get_allele_peak()$allele_repeat,
    type = "percentile", # "percentile" or "repeat"
    range = seq(0.1, 0.99, .10),
    col_prefix = "percentile"
  )

  repeat_test <- find_percentiles(
    repeats = test_distribution_df$repeats,
    signals = test_distribution_df$signal,
    index_peak_repeat = test_fragments[[1]]$get_allele_peak()$allele_repeat,
    type = "repeat", # "percentile" or "repeat"
    range = percentiles$percentile_0.2,
    col_prefix = "repeat"
  )

  # the values go full cirle
  testthat::expect_true(repeat_test[[1]] == "0.2")
})



# metrics ---------------------------------------

testthat::test_that("calculate metrics", {
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(example_data,
      dye_channel = "B",
      min_size_bp = 400
    )
  )
  config <- load_config()



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments, config
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      config,
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )



  # ungrouped

  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        test_fragments,
        config,
        grouped = FALSE
      )
    )
  )

  suppressMessages(
    test_metrics_ungrouped <- calculate_instability_metrics(
      fragments_list = test_fragments,
      peak_threshold = 0.05,
      # note the lower lim should be a negative value
      window_around_index_peak = c(-40, 40),
      percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
      repeat_range = c(1, 2, 3, 4, seq(6, 20, 2))
    )
  )

  testthat::expect_true(round(mean(test_metrics_ungrouped$expansion_index, na.rm = TRUE), 3) == 4.635)
  testthat::expect_true(all(is.na(test_metrics_ungrouped$average_repeat_change)))
  testthat::expect_true(round(mean(test_metrics_ungrouped$skewness, na.rm = TRUE), 5) == 0.10263)
  testthat::expect_true(test_fragments[[1]]$get_allele_peak()$allele_repeat == test_fragments[[1]]$get_index_peak()$index_repeat)
  
  
  #test override


  mock_override_df <- data.frame(
    unique_id = names(test_fragments)[1],
    override = c(120)

  )

  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        test_fragments,
        config,
        grouped = FALSE,
        index_override_dataframe = mock_override_df
      )
    )
  )


  suppressMessages(
    test_metrics_ungrouped_override <- calculate_instability_metrics(
      fragments_list = test_fragments,
      peak_threshold = 0.05,
      # note the lower lim should be a negative value
      window_around_index_peak = c(-40, 40),
      percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
      repeat_range = c(1, 2, 3, 4, seq(6, 20, 2))
    )
  )



  # grouped

  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        test_fragments,
        config,
        grouped = TRUE
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_fragments,
        peak_threshold = 0.05,
        # note the lower lim should be a negative value
        window_around_index_peak = c(-40, 40),
        percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
        repeat_range = c(1, 2, 3, 4, seq(6, 20, 2))
      )
    )
  )


  testthat::expect_true(round(mean(test_metrics_grouped$expansion_index, na.rm = TRUE), 3) == 5.729)
  testthat::expect_true(round(mean(test_metrics_grouped$average_repeat_change, na.rm = TRUE), 3) == 3.348)
  testthat::expect_true(round(mean(test_metrics_grouped$skewness, na.rm = TRUE), 5) == 0.10263)
  testthat::expect_true(test_fragments[[1]]$get_allele_peak()$allele_repeat != test_fragments[[1]]$get_index_peak()$index_repeat)







})

testthat::test_that("test index signal filters", {
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(example_data,
      dye_channel = "B",
      min_size_bp = 400
    )
  )
  config <- load_config()



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments, config
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      config,
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )

  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        test_fragments,
        config,
        grouped = TRUE
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_fragments,
        peak_threshold = 0.05,
        # note the lower lim should be a negative value
        window_around_index_peak = c(-40, 40),
        percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
        repeat_range = c(1, 2, 3, 4, seq(6, 20, 2)),
        index_modal_signal_threshold = 3000
      )
    )
  )
  testthat::expect_true(is.na(test_metrics_grouped[18,2]))
  testthat::expect_true(grepl("filter removed all index samples", test_metrics_grouped[19,2]))
  
})


testthat::test_that("test index sum filters", {
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(example_data,
      dye_channel = "B",
      min_size_bp = 400
    )
  )
  config <- load_config()



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  find_alleles(
    fragments_list = test_fragments, config
  )


  suppressWarnings(
    call_repeats(
      fragments_list = test_fragments,
      config,
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )

  suppressMessages(
    suppressWarnings(
      assign_index_peaks(
        test_fragments,
        config,
        grouped = TRUE
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_fragments,
        peak_threshold = 0.05,
        # note the lower lim should be a negative value
        window_around_index_peak = c(-40, 40),
        percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
        repeat_range = c(1, 2, 3, 4, seq(6, 20, 2)),
        index_signal_sum_threshold = 10000
      )
    )
  )
  testthat::expect_true(is.na(test_metrics_grouped[18,2]))
  testthat::expect_true(grepl("filter removed all index samples", test_metrics_grouped[19,2]))
  
})
