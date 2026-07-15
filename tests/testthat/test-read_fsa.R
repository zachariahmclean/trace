testthat::test_that("read_fsa builds unique_id from sample name + run id by default", {
  f <- system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr")
  fl <- read_fsa(f)

  expect_equal(fl[[1]]$unique_id, "FAC321_0000205983_Run_3130xl_2008-11-06_14-03_5134")
  # list name and the object's unique_id agree
  expect_equal(names(fl)[1], fl[[1]]$unique_id)
  # metadata is parsed and attached
  expect_equal(fl[[1]]$fsa_metadata$sample_name, "FAC321_0000205983")
})

testthat::test_that("read_fsa can fall back to the file name", {
  f <- system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr")
  fl <- read_fsa(f, unique_id_from_fsa = FALSE)

  expect_equal(fl[[1]]$unique_id, basename(f))
  expect_equal(names(fl)[1], basename(f))
})

testthat::test_that("read_fsa guarantees unique ids even with a repeated file", {
  f <- system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr")
  # the same file twice would otherwise collide
  fl <- read_fsa(c(f, f))

  ids <- vapply(fl, function(x) x$unique_id, character(1))
  expect_equal(length(unique(ids)), 2)
  expect_equal(unique(names(fl)), names(fl)) # names are unique
})
