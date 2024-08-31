

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
