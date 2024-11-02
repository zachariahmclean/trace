

# find alleles ---------------------------------

testthat::test_that("find_alleles", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class
  suppressWarnings(
    test_fragments <- peak_table_to_fragments(gm_raw,
      data_format = "genemapper5",
      dye_channel = "B",
      min_size_bp = 300
    )
  )
  

  find_alleles(
    fragments_list = test_fragments
  )

  testthat::expect_true(test_fragments[[1]]$get_allele_peak()$allele_size == 480.54)
  testthat::expect_true(test_fragments[[10]]$get_allele_peak()$allele_size == 476.53)

  allele_size <- vector("numeric", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    allele_size[i] <- test_fragments[[i]]$get_allele_peak()$allele_size
  }

  testthat::expect_true(all(!is.na(allele_size)))
})
