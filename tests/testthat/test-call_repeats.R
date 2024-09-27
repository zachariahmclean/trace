

testthat::test_that("call_repeats", {

suppressWarnings(
  test_fragments <- peak_table_to_fragments(example_data,
    data_format = "genemapper5",
    dye_channel = "B",
    min_size_bp = 300
)
)

  

  test_alleles <- find_alleles(
    fragments_list = test_fragments
  )

  suppressWarnings(
    suppressMessages(
      test_repeats <- call_repeats(
        fragments_list = test_alleles,
        repeat_calling_algorithm = "simple",
        assay_size_without_repeat = 87,
        repeat_size = 3
      )
    )
  )

  


  test_repeats_class <- vector("numeric", length(test_repeats))
  for (i in seq_along(test_repeats)) {
    test_repeats_class[i] <- class(test_repeats[[i]])[1]
  }


  testthat::expect_true(all(unique(test_repeats_class) == "fragments_repeats"))


  # force_whole_repeat_units
  suppressWarnings(
    suppressMessages(
      test_repeats_np <- call_repeats(
        fragments_list = test_alleles,
        force_whole_repeat_units = TRUE,
        assay_size_without_repeat = 87,
        repeat_size = 3
      )
    )
  )


  #TOODO come up with test for whole repeat units here
  
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
    fragments_list = peak_list
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

  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())

  suppressWarnings(
    find_ladders(
      fsa_list,
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


  peak_list <- find_fragments(fsa_list,
                              minimum_peak_signal = 20,
                              min_bp_size = 300
  )


  fragment_alleles <- find_alleles(
    fragments_list = peak_list
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



testthat::test_that("full pipeline repeat size algo", {

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

  suppressWarnings(
    find_ladders(fsa_list,
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


  peak_list <- find_fragments(fsa_list,
                              minimum_peak_signal = 20,
                              min_bp_size = 300
  )

  add_metadata(
    fragments_list = peak_list,
    metadata_data.frame = metadata
  )

  find_alleles(peak_list)

  suppressMessages(
    suppressWarnings(
      call_repeats(
        fragments_list = peak_list,
        repeat_calling_algorithm = "size_period"
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(100, 200),
  #             ylim = c(0,2000))


  # plot_fragments(test_repeats[1:4])
suppressMessages(
  suppressWarnings(
    assign_index_peaks(
      peak_list,
      grouped = TRUE
    )
  )
)


  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = peak_list,
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

  fsa_list <- c(H7_fsa_list, H8_fsa_list)

  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    H7_metadata)

# add a random systematic error to each sample
  # will use the number appended on the end as that for convience
  fragments_list <- lapply(fragments_list, function(x){
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
  

  find_alleles(fragments_list)
  suppressWarnings(
    call_repeats(fragments_list,
      correction = "batch")
    )

  testthat::expect_true(all.equal(c(seq(-4.5, 4.5, 1), seq(-4.5, 4.5, 1)), as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor))))
  
  # plot_batch_correction_samples(repeats_list, x_axis = "size", xlim = c(400, 470), n_facet_col = 2)
  # plot_batch_correction_samples(repeats_list, x_axis = "repeats", xlim = c(100, 130), n_facet_col = 2)
  # plot_batch_correction_samples(repeats_list, x_axis = "repeats", xlim = c(100, 130), n_facet_col = 1, sample_subset = "S-21-211")


})


testthat::test_that("batch correction", {

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata)
  
  find_alleles(fragments_list)
  suppressWarnings(
   call_repeats(fragments_list,
      correction = "batch")
    )

  testthat::expect_true(all.equal(c(rep(0.78526, 92), rep(-0.78526, 2)), round(as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor)), 5)))
  
  # plot_batch_correction_samples(fragments_list, selected_sample = 1, xlim = c(100, 115))


})



testthat::test_that("batch correction with no data in one batch", {

  different_batch <- vector("list", 1)
  names(different_batch) <- "test1"
  different_batch[[1]] <- cell_line_fsa_list[[1]]$clone()
  different_batch[[1]]$unique_id <- "test1"
  different_batch_metadata <- metadata[which(metadata$unique_id == names(cell_line_fsa_list[1])), ]
  different_batch_metadata$unique_id <- "test1"
  different_batch_metadata$batch_run_id <- NA_character_
  different_batch_metadata$batch_sample_id <- NA_character_


  metadata_modification_df <- rbind(metadata,  different_batch_metadata)

  fsa_list <- c(lapply(cell_line_fsa_list, function(x) x$clone()), different_batch)
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata_modification_df)
  

  find_alleles(fragments_list)

  suppressMessages(
    tryCatch({
      call_repeats(fragments_list,
        correction = "batch")
    },
      warning = function(w){
        assignment_warning <<- w
      }
    )
  )


  
  testthat::expect_true(class(assignment_warning)[1] == "simpleWarning")
  testthat::expect_true(grepl("'batch_run_id' 'NA' had no batch correction carried out", assignment_warning))


})




testthat::test_that("batch correction with a single sample id", {



  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata)
  
    fragments_list <- lapply(fragments_list, function(x){
    if(x$batch_sample_id %in% "S-21-211"){
      x$batch_sample_id <- NA_character_ 
    } 
      
      return(x)
     
      })
  

  find_alleles(fragments_list)
  suppressWarnings(call_repeats(fragments_list, correction = "batch"))
  
  
  # plot_batch_correction_samples(fragments_list, selected_sample = 1, xlim = c(100, 120))

  testthat::expect_true(all.equal(c(rep(0.86519, 92), rep(-0.86519, 2)), round(as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor)), 5)))



})


testthat::test_that("repeat correction", {

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata)
  
  for (i in seq_along(fragments_list)) {
    if(!is.na(fragments_list[[i]]$batch_sample_id) && fragments_list[[i]]$batch_sample_id == "S-21-211" ){
      fragments_list[[i]]$batch_sample_modal_repeat_length <- 120
    } else if(!is.na(fragments_list[[i]]$batch_sample_id) && fragments_list[[i]]$batch_sample_id == "S-21-212" ){
      fragments_list[[i]]$batch_sample_modal_repeat_length <- 122
    }
  }
  
  find_alleles(fragments_list)
  call_repeats(fragments_list,
    correction = "repeat")
  
    sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor)

  testthat::expect_true(all.equal(c(rep(0.78526, 92), rep(-0.78526, 2)), round(as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor)), 5)))
  
  # plot_batch_correction_samples(fragments_list, selected_sample = 1, xlim = c(100, 130))


})