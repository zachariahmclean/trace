

testthat::test_that("call_repeats", {
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

  suppressMessages(
    test_repeats <- call_repeats(
      fragments_list = test_alleles,
      repeat_calling_algorithm = "simple",
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )


  test_repeats_class <- vector("numeric", length(test_repeats))
  for (i in seq_along(test_repeats)) {
    test_repeats_class[i] <- class(test_repeats[[i]])[1]
  }


  testthat::expect_true(all(unique(test_repeats_class) == "fragments_repeats"))


  # nearest peak algo

  suppressMessages(
    test_repeats_np <- call_repeats(
      fragments_list = test_alleles,
      force_whole_repeat_units = TRUE,
      assay_size_without_repeat = 87,
      repeat_size = 3
    )
  )


  test_repeats_np_dif <- vector("list", length(test_repeats_np))
  for (i in seq_along(test_repeats_np)) {
    repeat_sizes <- test_repeats_np[[i]]$repeat_table_df$repeats

    lag <- vector("numeric", length(repeat_sizes))
    for (j in 2:length(repeat_sizes)) {
      lag[j] <- repeat_sizes[j] - repeat_sizes[j - 1]
    }

    test_repeats_np_dif[[i]] <- lag
  }

  all_integers <- sapply(test_repeats_np_dif, function(x) all(round(x, 10) %in% 0:200))

  testthat::expect_true(all(all_integers))

  
})


testthat::test_that("fft", {

  suppressWarnings(
    test_ladders <- find_ladders(cell_line_fsa_list[1],
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


  fragment_alleles <- find_alleles(
    fragments_list = peak_list,
    number_of_peaks_to_return = 1
  )


  suppressMessages(
    suppressWarnings(
      test_repeats_size_fft <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "fft"
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_repeats_size_simple <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "simple"
      )
    )
  )


  testthat::expect_true(
    identical( test_repeats_size_fft[[1]]$repeat_table_df[which(test_repeats_size_fft[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_fft[[1]]$repeat_table_df$repeats <140),"repeats"],
               test_repeats_size_simple[[1]]$repeat_table_df[which(test_repeats_size_simple[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_simple[[1]]$repeat_table_df$repeats <140),"repeats"]
    )
  )

})



testthat::test_that("repeat period", {

  suppressWarnings(
    test_ladders <- find_ladders(cell_line_fsa_list[1],
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


  fragment_alleles <- find_alleles(
    fragments_list = peak_list,
    number_of_peaks_to_return = 1
  )

  suppressMessages(
    suppressWarnings(
      test_repeats_size_period <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "size_period"
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_repeats_size_simple <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "simple"
      )
    )
  )


  testthat::expect_true(
    identical( test_repeats_size_period[[1]]$repeat_table_df[which(test_repeats_size_period[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_period[[1]]$repeat_table_df$repeats <140),"repeats"],
               test_repeats_size_simple[[1]]$repeat_table_df[which(test_repeats_size_simple[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_simple[[1]]$repeat_table_df$repeats <140),"repeats"]
               )
  )


})



testthat::test_that("full pipline repeat size algo", {
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
        repeat_calling_algorithm = "size_period"
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(100, 200),
  #             ylim = c(0,2000))


  # plot_fragments(test_repeats[1:4])
  # plot_size_standard_model(test_repeats)

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_repeats,
        grouped = TRUE,
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

  testthat::expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.85696, 0.70208, 0.57056, 1.00000, 1.17918, 1.10739, 1.00046)))
})



testthat::test_that("size standards with ids", {


  H7_metadata <- lapply(1:10, function(x){
    df <- metadata[which(metadata$unique_id %in% c("20230413_H07.fsa", "20230413_H08.fsa")), ]
    df$unique_id <- paste0(df$unique_id, "_", x)
    df$batch_run_id <- x 

     return(df)
  })

  H7_metadata <- do.call(rbind, H7_metadata)

  H7_fsa_list <- lapply(1:10, function(x) {
    y <- cell_line_fsa_list[["20230413_H07.fsa"]]$clone()
    y$unique_id <-  paste0("20230413_H07.fsa", "_", x)
    return(y)
  })
  names(H7_fsa_list) <- paste0("20230413_H07.fsa", "_", 1:10)

  

  H8_fsa_list <- lapply(1:10, function(x) {
    y <- cell_line_fsa_list[["20230413_H08.fsa"]]$clone()
    y$unique_id <-  paste0("20230413_H08.fsa", "_", x)
    return(y)
  })
  names(H8_fsa_list) <- paste0("20230413_H08.fsa", "_", 1:10)

  ladder_list <- find_ladders(c(H7_fsa_list, H8_fsa_list),
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(ladder_list, min_bp_size = 300)

  metadata_list <- add_metadata(fragments_list,
    H7_metadata)

# add a random systematic error to each sample
  # will use the number appended on the end as that for convience
  metadata_list <- lapply(metadata_list, function(x){
    x$peak_table_df$size <- x$peak_table_df$size + as.numeric(substr(x$unique_id, 18,19))
    x$trace_bp_df$size <- x$trace_bp_df$size + as.numeric(substr(x$unique_id, 18,19))
    #chose one of them to have the peaks swapped
    if(as.numeric(substr(x$unique_id, 18,19)) == 5){
      x$peak_table_df$height[7] <- x$peak_table_df$height[7] + 390
      x$trace_bp_df$signal[which(x$trace_bp_df$size == x$peak_table_df$size[7])] <- x$trace_bp_df$signal[which(x$trace_bp_df$size == x$peak_table_df$size[7])] + 390
    }
    else{
      #judt below other peak
      x$peak_table_df$height[7] <- x$peak_table_df$height[7] + 380
      x$trace_bp_df$signal[which(x$trace_bp_df$size == x$peak_table_df$size[7])] <- x$trace_bp_df$signal[which(x$trace_bp_df$size == x$peak_table_df$size[7])] + 380
    }

    return(x)
  })
  

  allele_list <- find_alleles(metadata_list)

  repeats_list <- call_repeats(allele_list,
    batch_correction = TRUE)

  testthat::expect_true(all(sapply(allele_list, function(x) x$.__enclos_env__$private$batch_correction_factor) == c(0:9, 0:9)))
  
  # plot_size_standard_samples(repeats_list, x_axis = "size", xlim = c(400, 470), n_facet_col = 2)
  # plot_size_standard_samples(repeats_list, x_axis = "repeats", xlim = c(100, 130), n_facet_col = 2)
  # plot_size_standard_samples(repeats_list, x_axis = "repeats", xlim = c(100, 130), n_facet_col = 1, sample_subset = "S-21-211")


})

