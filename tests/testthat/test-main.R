test_that("main trace", {

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

  frag_list <- suppressMessages(trace_main(fsa_list, min_bp_size = 300, grouped = TRUE, metadata_data.frame = metadata, show_progress_bar = FALSE))


  test_metrics_grouped <- calculate_instability_metrics(
    fragments_list = frag_list,
    peak_threshold = 0.1,
    window_around_index_peak = c(-25, 25)
  )

  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_signal > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$metrics_group_id)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_change / median(plot_data[[i]]$average_repeat_change[which(plot_data[[i]]$treatment == 0)])
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

  expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.87985, 0.69696, 0.48567)))


})

test_that("ladder fixing", {

  # ladder fixing
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

  example_list <- list(
    "20230413_A07.fsa" = data.frame(
      size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)
    )
   )

   frag_list <- suppressMessages(trace_main(fsa_list, ladder_df_list = example_list, min_bp_size = 300, grouped = TRUE, metadata_data.frame = metadata, show_progress_bar = FALSE))

})


test_that("main fragments", {

  fragment_list <- genemapper_table_to_fragments(example_data, dye_channel = "B",  min_size_bp = 300)

  frag_list <- suppressMessages(trace_main(fragment_list, grouped = TRUE, metadata_data.frame = metadata))

  
  test_metrics_grouped <- calculate_instability_metrics(
    fragments_list = frag_list,
    peak_threshold = 0.05,
    window_around_index_peak = c(-40, 40)
  )

  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_signal > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$metrics_group_id)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_change / median(plot_data[[i]]$average_repeat_change[which(plot_data[[i]]$treatment == 0)])
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

  expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.85621 , 0.69507  , 0.50166)))



})


test_that("main repeats", {

  fragment_list <- repeat_table_to_fragments(example_data_repeat_table, min_repeat = 71)

  frag_list <- suppressMessages(trace_main(fragment_list, grouped = TRUE, metadata_data.frame = metadata[1:17,], peak_region_size_gap_threshold =2))

  
  test_metrics_grouped <- calculate_instability_metrics(
    fragments_list = frag_list,
    peak_threshold = 0.05,
    window_around_index_peak = c(-40, 40)
  )

  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_signal > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$metrics_group_id)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_change / median(plot_data[[i]]$average_repeat_change[which(plot_data[[i]]$treatment == 0)])
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

  expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.89499 , 0.67094, 0.51213)))



})

test_that("main error",{
  wrong_metadata <- metadata
  colnames(wrong_metadata)[1] <- "splunique_id"

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  
  test_error <- tryCatch(
    trace_main(fsa_list, metadata_data.frame = wrong_metadata),
  error = function(e) e
)

expect_true("error" %in% class(test_error))

  
#warning
  
# fsa_list <- lapply(cell_line_fsa_list[1:10], function(x) x$clone())

# frag_list <- trace_main(fsa_list, metadata_data.frame = metadata)
  
})
