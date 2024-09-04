

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
  

  test_alleles <- find_alleles(
    fragments_list = test_fragments
  )

  testthat::expect_true(test_alleles[[1]]$get_allele_peak()$allele_size == 465.82)
  testthat::expect_true(test_alleles[[10]]$get_allele_peak()$allele_size == 505.21)

  allele_size <- vector("numeric", length(test_alleles))
  for (i in seq_along(test_alleles)) {
    allele_size[i] <- test_alleles[[i]]$get_allele_peak()$allele_size
  }

  testthat::expect_true(all(!is.na(allele_size[-c(4,47)])))
})
