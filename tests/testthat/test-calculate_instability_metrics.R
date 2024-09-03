# percentiles

testthat::test_that("percentiles", {
  gm_raw <- trace::example_data
  test_sample <- "20230413_A01.fsa"

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
    data_format = "genemapper5",
    # peak_size_col = "size",
    # peak_height_col = "signal",
    dye_channel = "B",
    min_size_bp = 350
  )

  test_main_peaks <- find_alleles(test_fragments[1],
                                  number_of_peaks_to_return = 1,
                                  peak_region_size_gap_threshold = 6,
                                  peak_region_height_threshold_multiplier = 1)

  test_fragments_repeats_simple <- call_repeats(
    test_main_peaks,
    repeat_calling_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3
  )


  test_distribution_df <- test_fragments_repeats_simple[[1]]$repeat_table_df
  test_distribution_df <- test_distribution_df[which(test_distribution_df$repeats > 113), ]


  percentiles <- find_percentiles(
    repeats = test_distribution_df$repeats,
    heights = test_distribution_df$height,
    index_peak_repeat = test_fragments_repeats_simple[[1]]$get_alleles()$allele_1_repeat,
    type = "percentile", # "percentile" or "repeat"
    range = seq(0.1, 0.99, .10),
    col_prefix = "percentile"
  )

  repeat_test <- find_percentiles(
    repeats = test_distribution_df$repeats,
    heights = test_distribution_df$height,
    index_peak_repeat = test_fragments_repeats_simple[[1]]$get_alleles()$allele_1_repeat,
    type = "repeat", # "percentile" or "repeat"
    range = percentiles$percentile_0.2,
    col_prefix = "repeat"
  )

  # the values go full cirle
  testthat::expect_true(repeat_test[[1]] == 0.2)
})



# metrics ---------------------------------------

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



  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )

  test_alleles <- find_alleles(
    fragments_list = test_metadata,
    number_of_peaks_to_return = 1,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )


  suppressWarnings(
    test_repeats <- call_repeats(
      fragments_list = test_alleles,
      repeat_calling_algorithm = "simple",
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )



  # ungrouped

  suppressMessages(
    suppressWarnings(
      test_assignment_ungrouped <- assign_index_peaks(
        test_repeats,
        grouped = FALSE
      )
    )
  )

  suppressMessages(
    test_metrics_ungrouped <- calculate_instability_metrics(
      fragments_list = test_assignment_ungrouped,
      peak_threshold = 0.05,
      # note the lower lim should be a negative value
      window_around_index_peak = c(-40, 40),
      percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
      repeat_range = c(1, 2, 3, 4, seq(6, 20, 2))
    )
  )

  testthat::expect_true(round(mean(test_metrics_ungrouped$expansion_index, na.rm = TRUE), 3) == 4.956)
  testthat::expect_true(all(is.na(test_metrics_ungrouped$average_repeat_gain)))
  testthat::expect_true(round(mean(test_metrics_ungrouped$skewness, na.rm = TRUE), 5) == -0.00905)
  testthat::expect_true(test_assignment_ungrouped[[1]]$get_alleles()$allele_1_repeat == test_assignment_ungrouped[[1]]$get_index_peak()$index_repeat)
  testthat::expect_true(all(sapply(test_assignment_ungrouped, function(x) x$.__enclos_env__$private$assigned_index_peak_used)))
  #test override


  mock_override_df <- data.frame(
    unique_id = names(test_repeats)[1],
    override = c(120)

  )

  suppressMessages(
    suppressWarnings(
      test_assignment_override <- assign_index_peaks(
        test_repeats,
        grouped = FALSE,
        index_override_dataframe = mock_override_df
      )
    )
  )


  suppressMessages(
    test_metrics_ungrouped_override <- calculate_instability_metrics(
      fragments_list = test_assignment_override,
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
      test_assignment_grouped <- assign_index_peaks(
        test_repeats,
        grouped = TRUE
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_assignment_grouped,
        peak_threshold = 0.05,
        # note the lower lim should be a negative value
        window_around_index_peak = c(-40, 40),
        percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
        repeat_range = c(1, 2, 3, 4, seq(6, 20, 2))
      )
    )
  )


  testthat::expect_true(round(mean(test_metrics_grouped$expansion_index, na.rm = TRUE), 3) == 6.727)
  testthat::expect_true(round(mean(test_metrics_grouped$average_repeat_gain, na.rm = TRUE), 3) == 4.14)
  testthat::expect_true(round(mean(test_metrics_grouped$skewness, na.rm = TRUE), 5) == -0.00905)
  testthat::expect_true(test_assignment_grouped[[1]]$get_alleles()$allele_1_repeat != test_assignment_grouped[[1]]$get_index_peak()$index_repeat)
  testthat::expect_true(all(sapply(test_assignment_grouped, function(x) x$.__enclos_env__$private$assigned_index_peak_used)))




})
