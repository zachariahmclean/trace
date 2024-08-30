

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
