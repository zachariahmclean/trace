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

  })

