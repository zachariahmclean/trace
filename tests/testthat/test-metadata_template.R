testthat::test_that("generate_metadata_template returns the add_metadata columns with batch_run_id prefilled", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  template <- generate_metadata_template(fsa_list)

  expect_s3_class(template, "data.frame")
  expect_equal(nrow(template), length(fsa_list))

  # exactly the columns add_metadata() expects must be present
  metadata_params <- c(
    "unique_id", "metrics_group_id", "metrics_baseline_control",
    "batch_run_id", "batch_sample_id", "batch_sample_modal_repeat"
  )
  expect_true(all(metadata_params %in% names(template)))

  # provenance columns present
  expect_true(all(c("fsa_run_date", "fsa_well", "fsa_plate", "fsa_instrument",
                    "fsa_capillary", "fsa_sample_name") %in% names(template)))

  # batch_run_id is prefilled (not NA) and matches the parsed run id
  expect_false(any(is.na(template$batch_run_id)))
  m <- extract_fsa_metadata(fsa_list)
  expect_setequal(template$batch_run_id, m$run_id)

  # the to-be-filled columns are blank
  expect_true(all(is.na(template$metrics_group_id)))
  expect_true(all(is.na(template$batch_sample_id)))
  expect_true(all(is.na(template$batch_sample_modal_repeat)))
})

testthat::test_that("generate_metadata_template round-trips through add_metadata without a missing-column warning", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  template <- generate_metadata_template(fsa_list)

  status <- add_metadata(fsa_list, template)
  # add_metadata warns "...not supplied..." only when required columns are missing
  warn_text <- paste(unlist(status$warning_message), collapse = " ")
  expect_false(grepl("not supplied", warn_text))
})

testthat::test_that("generate_metadata_template accepts a file path", {
  f <- system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr")
  template <- generate_metadata_template(f)

  expect_equal(nrow(template), 1)
  expect_equal(template$batch_run_id, "Run_3130xl_2008-11-06_14-03_5134")
})

testthat::test_that("generate_metadata_template writes a csv when requested", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  out <- generate_metadata_template(fsa_list, output_csv = tmp)
  expect_true(file.exists(tmp))

  reread <- utils::read.csv(tmp, colClasses = "character")
  expect_equal(nrow(reread), nrow(out))
  expect_true("batch_run_id" %in% names(reread))
})

testthat::test_that("generate_metadata_template errors on invalid input", {
  expect_error(generate_metadata_template(42))
})
