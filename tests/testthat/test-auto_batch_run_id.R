testthat::test_that("set_batch_run_id_from_fsa fills batch_run_id from the fsa file", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  config <- load_config()

  status <- set_batch_run_id_from_fsa(fsa_list, config)

  # every fsa sample now has the RunN-derived run id
  expect_equal(fsa_list[[1]]$batch_run_id, "Run_MRSPHILLIPS2_2023-04-14_11-30_0063")
  expect_false(any(vapply(fsa_list, function(x) is.na(x$batch_run_id), logical(1))))
  # no mismatch -> okay status
  expect_equal(status$status, "okay")
})

testthat::test_that("set_batch_run_id_from_fsa overrides a conflicting value and warns", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  fsa_list[[1]]$batch_run_id <- "user_typed_wrong_run"
  config <- load_config()

  status <- set_batch_run_id_from_fsa(fsa_list, config)

  expect_equal(fsa_list[[1]]$batch_run_id, "Run_MRSPHILLIPS2_2023-04-14_11-30_0063")
  expect_equal(status$status, "warning")
  expect_true(any(grepl("user_typed_wrong_run", unlist(status$warning_message))))
})

testthat::test_that("batch_run_id_tag = 'RUND_RUNT' builds a date_time id", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  config <- load_config()
  config$batch_run_id_tag <- "RUND_RUNT"

  set_batch_run_id_from_fsa(fsa_list, config)
  expect_equal(fsa_list[[1]]$batch_run_id, "20230414_113050")
})

testthat::test_that("trace() errors early when correction is set but batch_sample_id is missing", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  expect_error(
    trace(fsa_list, correction = "batch"),
    "batch_sample_id"
  )
})

testthat::test_that("trace() errors early when grouped = TRUE but metrics_group_id is missing", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  expect_error(
    trace(fsa_list, grouped = TRUE),
    "metrics_group_id"
  )
})
