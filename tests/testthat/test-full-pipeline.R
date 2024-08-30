# peak_table_to_fragments -------------------------------------------------

testthat::test_that("peak_table_to_fragments", {
  gm_raw <- trace::example_data

  test_fragments <- peak_table_to_fragments(
    gm_raw,
    data_format = "genemapper5",
    peak_size_col = "Size",
    peak_height_col = "Height",
    unique_id = "Sample.File.Name",
    dye_col = "Dye.Sample.Peak",
    dye_channel = "B",
    allele_col = "Allele",
    min_size_bp = 100,
    max_size_bp = 1000
  )

  test_sample <- test_fragments[[1]]

  test_fragments_classes <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_classes[i] <- class(test_fragments[[i]])[1]
  }

  test_fragments_ids <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_ids[i] <- test_fragments[[i]]$unique_id
  }

  # it's a list
  testthat::expect_true(class(test_fragments) == "list")
  # we only expect one class
  testthat::expect_true(length(unique(test_fragments_classes)) == 1)
  # everything in the list is the class we expect
  testthat::expect_true(unique(test_fragments_classes) == "fragments_repeats")
  # no missing unique ids
  testthat::expect_false(any(is.na(test_fragments_ids)))
})

# repeat_table_to_fragments -------------------------------------------------

testthat::test_that("repeat_table_to_fragments", {
  repeat_table <- trace::example_data_repeat_table

  test_fragments <- repeat_table_to_repeats(
    repeat_table,
    repeat_col = "repeats",
    frequency_col = "height",
    unique_id = "unique_id"
  )

  test_sample <- test_fragments[[1]]

  test_fragments_classes <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_classes[i] <- class(test_fragments[[i]])[1]
  }

  test_fragments_ids <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_ids[i] <- test_fragments[[i]]$unique_id
  }

  # it's a list
  testthat::expect_true(class(test_fragments) == "list")
  # we only expect one class
  testthat::expect_true(length(unique(test_fragments_classes)) == 1)
  # everything in the list is the class we expect
  testthat::expect_true(unique(test_fragments_classes) == "fragments_repeats")
  # no missing unique ids
  testthat::expect_false(any(is.na(test_fragments_ids)))
})


# add metadata -------------------------------------------------

testthat::test_that("add_metadata", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
    data_format = "genemapper5",
    # peak_size_col = "size",
    # peak_height_col = "signal",
    dye_channel = "B"
  )

  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )


  # batch_run_id assigned
  test_fragments_batch_run_id <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_batch_run_id[i] <- test_metadata[[i]]$batch_run_id
  }

  testthat::expect_true(all(test_fragments_batch_run_id == 20230414))

  # metrics_group_id assigned
  test_fragments_metrics_group_id <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_metrics_group_id[i] <- test_metadata[[i]]$metrics_group_id
  }

  testthat::expect_true(all(!is.na(test_fragments_metrics_group_id)))

  # metrics_baseline_control assigned
  test_fragments_index <- vector("logical", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_index[i] <- test_metadata[[i]]$metrics_baseline_control
  }

  index_samples <- which(metadata$metrics_baseline_control == TRUE)

  testthat::expect_true(all(as.logical(test_fragments_index[index_samples])))
  testthat::expect_true(unique(test_fragments_index[which(!seq_along(test_fragments_index) %in% index_samples)]) == FALSE)

  # skip columns
  test_metadata_skip <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = NA,
    metrics_group_id = NA,
    metrics_baseline_control = NA
  )



  })



testthat::test_that("add_metadata missing", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
    data_format = "genemapper5",
    # peak_size_col = "size",
    # peak_height_col = "signal",
    dye_channel = "B"
  )

  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = NA
  )

  testthat::expect_true(is.na(test_metadata[[1]]$batch_run_id))

  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = NA
  )

  testthat::expect_true(is.na(test_metadata[[1]]$metrics_group_id))

  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = NA
  )

  testthat::expect_false(test_metadata[[1]]$metrics_baseline_control)


  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = "metrics_baseline_control"
  )



  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = "metrics_baseline_control"
  )

})


# find alleles ---------------------------------

testthat::test_that("find_alleles", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
    data_format = "genemapper5",
    dye_channel = "B"
  )

  test_alleles <- find_alleles(
    fragments_list = test_fragments,
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  allele_1_size <- vector("numeric", length(test_alleles))
  for (i in seq_along(test_alleles)) {
    allele_1_size[i] <- test_alleles[[i]]$get_alleles()$allele_1_size
  }

  testthat::expect_true(all(!is.na(allele_1_size)))
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


# remove fragments --------------------------------------------------------

testthat::test_that("remove fragments", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(
    gm_raw,
    data_format = "genemapper5",
    dye_channel = "B"
  )

  all_fragment_names <- names(test_fragments)
  samples_to_remove <- all_fragment_names[c(1, length(all_fragment_names))]

  samples_removed <- remove_fragments(test_fragments, samples_to_remove)

  testthat::expect_true(all(!names(samples_removed) %in% samples_to_remove))
})


# plot fragments ----------------------------------------------------------


testthat::test_that("full pipline", {
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
        fragments_list = fragment_alleles,
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(400, 550),
  #             ylim = c(0,2000))


  # plot_fragments(test_repeats[1:4])
  # plot_size_standard_model(test_repeats)

  suppressMessages(
    suppressWarnings(
      test_assignment <- assign_index_peaks(
        test_repeats,
        grouped = TRUE
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_assignment,
        peak_threshold = 0.05,
        window_around_index_peak = c(-40, 40)
      )
    )
  )


  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_height > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$metrics_group_id)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_gain / median(plot_data[[i]]$average_repeat_gain[which(plot_data[[i]]$treatment == 0)])
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

  expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.85831, 0.70219, 0.57056, 1.00000, 1.17666, 1.10977, 1.00459)))
})

