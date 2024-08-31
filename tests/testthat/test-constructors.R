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

