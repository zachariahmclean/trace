testthat::test_that("parse_fsa_metadata extracts fields from a 3730xl fsa", {
  m <- parse_fsa_metadata(cell_line_fsa_list[[1]]$fsa)

  expect_equal(m$run_id, "Run_MRSPHILLIPS2_2023-04-14_11-30_0063")
  expect_equal(m$run_datetime, "2023-04-14 11:30:50")
  expect_equal(m$instrument, "MRSPHILLIPS2-15104-012")
  expect_equal(m$instrument_serial, "15104-012")
  expect_equal(m$well, "A7")
  expect_equal(m$capillary, 63)
  expect_equal(m$plate, "GT-2023-04-14-GT2")
  expect_equal(m$sample_name, "20230413_A07")
  expect_true(any(m$dye_names == "LIZ"))
  expect_true(m$n_offscale > 0)
})

testthat::test_that("parse_fsa_metadata generalises to another instrument (seqinr 3130xl example)", {
  ex <- seqinr::read.abif(system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr"))
  m <- parse_fsa_metadata(ex)

  expect_equal(m$run_id, "Run_3130xl_2008-11-06_14-03_5134")
  expect_equal(m$well, "B2")
  expect_equal(m$capillary, 4)
})

testthat::test_that("parse_fsa_metadata returns NA (not an error) when tags are missing", {
  m <- parse_fsa_metadata(list(Data = list()))
  expect_true(is.na(m$run_id))
  expect_true(is.na(m$instrument))
  expect_equal(m$n_offscale, 0L)
  expect_equal(length(m$dye_names), 0L)

  # falls back to zero-padded RUND_RUNT when RunN is absent
  m2 <- parse_fsa_metadata(list(Data = list(
    RUND.1 = list(2021L, 1L, 22L),
    RUNT.1 = list(9L, 5L, 3L, 0L)
  )))
  expect_equal(m2$run_id, "20210122_090503")
  expect_equal(m2$run_datetime, "2021-01-22 09:05:03")
})

testthat::test_that("parse_fsa_metadata tolerates a NULL / empty object", {
  expect_silent(m <- parse_fsa_metadata(NULL))
  expect_true(is.na(m$run_id))
})

testthat::test_that("extract_fsa_metadata returns one row per sample with expected columns", {
  out <- extract_fsa_metadata(cell_line_fsa_list[1:2])

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_true(all(c(
    "unique_id", "run_id", "run_datetime", "instrument", "instrument_serial",
    "plate", "well", "capillary", "sample_name", "dyes", "n_offscale", "n_saturated"
  ) %in% names(out)))
  expect_equal(out$unique_id, names(cell_line_fsa_list)[1:2])
})
