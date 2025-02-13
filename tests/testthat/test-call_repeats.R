

testthat::test_that("call_repeats", {

suppressWarnings(
  test_fragments <- genemapper_table_to_fragments(example_data,
    dye_channel = "B",
    min_size_bp = 300
)
)

  

  find_alleles(
    fragments_list = test_fragments
  )

  suppressWarnings(
    suppressMessages(
      call_repeats(
        fragments_list = test_fragments,
        assay_size_without_repeat = 87,
        repeat_size = 3
      )
    )
  )

  


  test_repeats_class <- vector("numeric", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_repeats_class[i] <- class(test_fragments[[i]])[1]
  }


  testthat::expect_true(all(unique(test_repeats_class) == "fragments"))


  # force_whole_repeat_units
  suppressWarnings(
    suppressMessages(
      call_repeats(
        fragments_list = test_fragments,
        force_whole_repeat_units = TRUE,
        assay_size_without_repeat = 87,
        repeat_size = 3
      )
    )
  )


  #TOODO come up with test for whole repeat units here
  
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


find_alleles(
    fragments_list = peak_list
  )

  test_repeats_size_period <- lapply(peak_list, function(x) x$clone())


  suppressMessages(
    suppressWarnings(
       call_repeats(
        fragments_list = test_repeats_size_period,
        force_repeat_pattern = TRUE,
        force_repeat_pattern_size_window = 0.5
      )
    )
  )

  test_repeats_size_none <- lapply(peak_list, function(x) x$clone())

  suppressMessages(
    suppressWarnings(
     call_repeats(
        fragments_list = test_repeats_size_none,
        force_repeat_pattern = FALSE )
    )
  )


  testthat::expect_true(
    identical( test_repeats_size_period[[1]]$repeat_table_df[which(test_repeats_size_period[[1]]$repeat_table_df$repeats > 118 & test_repeats_size_period[[1]]$repeat_table_df$repeats <140),"repeats"],
               test_repeats_size_none[[1]]$repeat_table_df[which(test_repeats_size_none[[1]]$repeat_table_df$repeats > 118 & test_repeats_size_none[[1]]$repeat_table_df$repeats <140),"repeats"]
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
        force_repeat_pattern = TRUE
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

  testthat::expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.86154, 0.73268, 0.55720)))
})



testthat::test_that("batch correction", {

  fsa_list <- lapply(cell_line_fsa_list[16:19], function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)


  add_metadata(fragments_list,
    metadata[16:19, ])

  
  
  find_alleles(fragments_list)
  suppressMessages(
  suppressWarnings(
   call_repeats(fragments_list,
      correction = "batch")
    )
  )
  testthat::expect_true(all.equal(
    c(rep(0.73299, 2), rep(-0.73299, 2)), 
    round(as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor)), 5)
  ))
  
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


  metadata_modification_df <- rbind(metadata[16:19, ],  different_batch_metadata)

  fsa_list <- c(lapply(cell_line_fsa_list[16:19], function(x) x$clone()), different_batch)
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata_modification_df)
  

  find_alleles(fragments_list)

  suppressMessages(
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
)

  
  testthat::expect_true(class(assignment_warning)[1] == "simpleWarning")
  testthat::expect_true(grepl("'batch_run_id' 'NA' had no batch correction carried out", assignment_warning))


})




testthat::test_that("batch correction with a single sample id", {



  fsa_list <- lapply(cell_line_fsa_list[16:19], function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata[16:19, ])
  
    fragments_list <- lapply(fragments_list, function(x){
    if(x$batch_sample_id %in% "S-21-211"){
      x$batch_sample_id <- NA_character_ 
    } 
      
      return(x)
     
      })
  

  find_alleles(fragments_list)
  suppressMessages(
  suppressWarnings(call_repeats(fragments_list, correction = "batch"))
  )
  
  # plot_batch_correction_samples(fragments_list, selected_sample = 1, xlim = c(100, 120))

  testthat::expect_true(all.equal(
    c(rep(0.8113, 2), rep(-0.8113, 2)), 
    round(as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$batch_correction_factor)), 5)
  ))



})


testthat::test_that("repeat correction", {

  fsa_list <- lapply(cell_line_fsa_list[16:19], function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)

  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata[16:19, ])
    
  find_alleles(fragments_list)

  suppressMessages(
    suppressWarnings(
      call_repeats(
        fragments_list = fragments_list,
        correction  = "repeat"
      )
    )
  )

  # plot_batch_correction_samples(fragments_list, 1, c(100,130))

  # plot_repeat_correction_model(fragments_list)


  testthat::expect_true(all.equal(
    c(10.51614, 10.58441, 11.09593, 11.17947), 
    round(as.numeric(sapply(fragments_list, function(x) x$.__enclos_env__$private$repeat_correction_factor)), 5)
  ))


  # extract model summary
  correction_summary <- extract_repeat_correction_summary(fragments_list)

  testthat::expect_true(is.data.frame(correction_summary))
  testthat::expect_true(all.equal(round(correction_summary$abs_avg_residual, 5), c(0.01585, 0.01585, 0.01988, 0.01704)))

})


testthat::test_that("repeat correction with one sample off warning", {

  fsa_list <- lapply(cell_line_fsa_list[16:19], function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)
  
  # make peak before just bigger 
  fsa_list[[1]]$trace_bp_df[which(round(fsa_list[[1]]$trace_bp_df$size, 2) == 418.25), "signal"] <- 4000
  
  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  add_metadata(fragments_list,
    metadata[16:19, ])
    
  find_alleles(fragments_list)


  suppressMessages(
    tryCatch({
      call_repeats(
        fragments_list = fragments_list,
        correction  = "repeat"
      )
    },
      warning = function(w){
        assignment_warning <<- w
      }
    )
  )

  # plot_batch_correction_samples(fragments_list, 1, c(100,115))

  # plot_repeat_correction_model(fragments_list)


  testthat::expect_true(class(assignment_warning)[1] == "simpleWarning")
  testthat::expect_true(grepl("had no batch correction carried out", assignment_warning))



})


testthat::test_that("repeat correction one run missing", {

  fsa_list <- lapply(cell_line_fsa_list[16:19], function(x) x$clone())
  find_ladders(fsa_list,
          show_progress_bar = FALSE)
  
  
  fragments_list <- find_fragments(fsa_list, min_bp_size = 300)

  metadata_2 <- metadata[16:19, ]

  metadata_2$batch_sample_modal_repeat <- ifelse(metadata_2$batch_run_id == "20230414", NA_real_,  metadata_2$batch_sample_modal_repeat)


  add_metadata(fragments_list,
    metadata_2)
    
  find_alleles(fragments_list)


  suppressMessages(
    tryCatch({
      call_repeats(
        fragments_list = fragments_list,
        correction  = "repeat"
      )
    },
      error = function(e){
        assignment_error <<- e
      }
    )
  )

  testthat::expect_true(class(assignment_error)[1] == "simpleError")
  testthat::expect_true(grepl("no samples with 'batch_sample_modal_repeat'", assignment_error))

})



