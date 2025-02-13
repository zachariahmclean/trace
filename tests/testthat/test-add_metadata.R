# add metadata -------------------------------------------------

testthat::test_that("add_metadata", {
  test_fragments <- genemapper_table_to_fragments(example_data,
    dye_channel = "B"
  )

  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata
  )


  # batch_run_id assigned
  test_fragments_batch_run_id <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_batch_run_id[i] <- test_fragments[[i]]$batch_run_id
  }

  testthat::expect_true(all(test_fragments_batch_run_id %in% c(20230414, 20220630)))

  # metrics_group_id assigned
  test_fragments_metrics_group_id <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_metrics_group_id[i] <- test_fragments[[i]]$metrics_group_id
  }

  testthat::expect_true(all(!is.na(test_fragments_metrics_group_id)))

  # metrics_baseline_control assigned
  test_fragments_index <- vector("logical", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_index[i] <- test_fragments[[i]]$metrics_baseline_control
  }

  index_samples <- which(metadata$metrics_baseline_control == TRUE)

  testthat::expect_true(all(as.logical(test_fragments_index[index_samples])))
  testthat::expect_true(unique(test_fragments_index[which(!seq_along(test_fragments_index) %in% index_samples)]) == FALSE)

  # skip columns
  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = NA,
    metrics_group_id = NA,
    metrics_baseline_control = NA
  )



  })



testthat::test_that("add_metadata missing", {
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(example_data,
      dye_channel = "B",
      min_size_bp = 300
    )
  )
  

  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = NA
  )

  testthat::expect_true(is.na(test_fragments[[1]]$batch_run_id))

  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = NA
  )

  testthat::expect_true(is.na(test_fragments[[1]]$metrics_group_id))

  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = NA
  )

  testthat::expect_false(test_fragments[[1]]$metrics_baseline_control)


  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = "metrics_baseline_control"
  )



  add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    batch_run_id = "batch_run_id",
    metrics_group_id = "metrics_group_id",
    metrics_baseline_control = "metrics_baseline_control"
  )

})



# metadata transfer



testthat::test_that("metadata transfer", {
  
  
  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())


  suppressWarnings(
    find_ladders(fsa_list,
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )
  suppressWarnings(
    add_metadata(fsa_list, metadata
    )
  )

  peak_list <- find_fragments(fsa_list,
    minimum_peak_signal = 20,
    min_bp_size = 300
  )

  testthat::expect_true(fsa_list[[1]]$unique_id == peak_list[[1]]$unique_id)
  testthat::expect_true(fsa_list[[1]]$batch_run_id == peak_list[[1]]$batch_run_id)
  testthat::expect_true(fsa_list[[1]]$metrics_group_id == peak_list[[1]]$metrics_group_id)
  testthat::expect_true(fsa_list[[1]]$metrics_baseline_control == peak_list[[1]]$metrics_baseline_control)
})

