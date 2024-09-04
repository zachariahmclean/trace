# add metadata -------------------------------------------------

testthat::test_that("add_metadata", {
  test_fragments <- peak_table_to_fragments(example_data,
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

  testthat::expect_true(all(test_fragments_batch_run_id %in% c(20230414, 20220630)))

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
  test_fragments <- peak_table_to_fragments(example_data,
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
