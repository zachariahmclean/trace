testthat::test_that("qc_report returns one row per sample with the expected columns", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  qc <- qc_report(processed)

  expect_s3_class(qc, "data.frame")
  expect_equal(nrow(qc), length(processed))
  expect_true(all(c(
    "unique_id", "batch_run_id", "run_date", "instrument", "capillary", "well",
    "plate", "ladder_avg_rsq", "ladder_min_rsq", "n_peaks", "modal_size",
    "modal_signal", "modal_prominence", "modal_saturated",
    "saturation_in_window", "qc_flags", "qc_pass"
  ) %in% names(qc)))
  expect_type(qc$qc_pass, "logical")
})

testthat::test_that("qc_report flags a sample with too few peaks", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  # force one sample to look empty
  processed[[1]]$peak_table_df <- processed[[1]]$peak_table_df[1:2, , drop = FALSE]

  qc <- qc_report(processed)
  row <- qc[qc$unique_id == processed[[1]]$unique_id, ]
  expect_true(grepl("few_peaks", row$qc_flags))
  expect_false(row$qc_pass)
})

testthat::test_that("qc_report flags a broken ladder via the worst-segment rsq", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  # a generous min_rsq threshold should flag the worst-fitting sample(s)
  qc_strict <- qc_report(processed, qc_min_rsq = 0.99999)
  expect_true(any(grepl("low_ladder_rsq", qc_strict$qc_flags)))

  # a permissive threshold should clear the ladder flag
  qc_loose <- qc_report(processed, qc_min_rsq = 0)
  expect_false(any(grepl("low_ladder_rsq", qc_loose$qc_flags)))
})

testthat::test_that("qc_report thresholds can be overridden via ...", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  # an impossibly high modal-signal requirement should flag every sample
  qc <- qc_report(processed, qc_min_modal_signal = 1e9)
  expect_true(all(grepl("low_signal", qc$qc_flags)))
})
