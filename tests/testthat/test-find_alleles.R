

# find alleles ---------------------------------

testthat::test_that("find_alleles", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(gm_raw,
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

 
testthat::test_that("find_alleles two alleles", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(gm_raw,
      dye_channel = "B",
      min_size_bp = 100
    )
  )

  find_alleles(
    fragments_list = test_fragments,
    number_of_alleles = 2
  )

  # alleles <- extract_alleles(test_fragments)

  testthat::expect_true(test_fragments[[1]]$get_allele_peak()$allele_size == 480.54)
  testthat::expect_true(test_fragments[[1]]$get_allele_peak()$allele_2_size == 131.92)
  testthat::expect_true(test_fragments[[10]]$get_allele_peak()$allele_size == 476.53)
  testthat::expect_true(test_fragments[[10]]$get_allele_peak()$allele_2_size == 132.02)

  allele_size <- vector("numeric", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    allele_size[i] <- test_fragments[[i]]$get_allele_peak()$allele_size
  }

  testthat::expect_true(all(!is.na(allele_size)))
})

 
testthat::test_that("find_alleles warning", {
  gm_raw <- trace::example_data
  metadata <- trace::metadata
  # Save raw data as a fragment class
  suppressWarnings(
    test_fragments <- genemapper_table_to_fragments(gm_raw,
      dye_channel = "B",
      min_size_bp = 300
    )
  )

  lapply(test_fragments[1:2], function(x){
    x$peak_table_df <- x$peak_table_df[rep(FALSE, nrow(x$peak_table_df)), , drop = FALSE]
    invisible()
  })
  

  find_alleles_output <- find_alleles(
    test_fragments
  )

  testthat::expect_true(find_alleles_output$status == "warning")
})
