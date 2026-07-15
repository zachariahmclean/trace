testthat::test_that("generate_dashboard writes an html file with the expected structure", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  out <- tempfile(fileext = ".html")
  result <- generate_dashboard(processed, output_file = out, open_browser = FALSE)

  expect_equal(result, out)
  expect_true(file.exists(out))

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_true(grepl("dashboard-grid", html))
  expect_true(grepl("qc-table", html))
  expect_true(grepl("js-plotly-plot|plotly html-widget", html))

  # one tile per sample
  tile_count <- length(regmatches(html, gregexpr("class=\"dashboard-tile", html))[[1]])
  expect_equal(tile_count, length(processed))
})

testthat::test_that("generate_dashboard sorts failed-QC samples first and flags them", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  # force one sample to fail qc_report's few_peaks check
  processed[[1]]$peak_table_df <- processed[[1]]$peak_table_df[1:2, , drop = FALSE]
  failing_id <- processed[[1]]$unique_id

  out <- tempfile(fileext = ".html")
  generate_dashboard(processed, output_file = out, open_browser = FALSE)

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_true(grepl("qc-fail", html))
  expect_true(grepl("few_peaks", html))

  # the failing sample's tile should be the first one plotted
  headers <- regmatches(html, gregexpr("<h4>[^<]*</h4>", html))[[1]]
  expect_equal(headers[1], paste0("<h4>", failing_id, "</h4>"))
})

testthat::test_that("generate_dashboard respects sample_subset", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  keep_ids <- names(processed)[1:2]
  out <- tempfile(fileext = ".html")
  generate_dashboard(processed, output_file = out, sample_subset = keep_ids, open_browser = FALSE)

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  tile_count <- length(regmatches(html, gregexpr("class=\"dashboard-tile", html))[[1]])
  expect_equal(tile_count, 2)
})

testthat::test_that("generate_dashboard errors when sample_subset leaves nothing to plot", {
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
  processed <- trace(fsa_list)

  out <- tempfile(fileext = ".html")
  expect_error(
    generate_dashboard(processed, output_file = out, sample_subset = "not_a_real_id", open_browser = FALSE),
    "No samples to plot"
  )
})
