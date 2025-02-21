testthat::test_that("find ladder peaks", {

  test_processed <- data.frame(signal = cell_line_fsa_list[[1]]$fsa$Data$DATA.105, scan = 0:(length(cell_line_fsa_list[[1]]$fsa$Data$DATA.105) - 1))
  test_processed <- test_processed[which(test_processed$scan >= which.max(test_processed$signal)), ]
  test_processed$detrended_signal <- detrend_signal(test_processed$signal)
  test_processed$smoothed_signal <- pracma::savgol(
    test_processed$detrended_signal,
    21
  )


  ladder_sizes <- c(50, 75, 100, 139, 150, 160, 200, 300, 350, 400, 450, 490, 500)


  test_ladder_peaks <- find_ladder_peaks(
    test_processed,
    length(ladder_sizes),
    minimum_ladder_signal = NULL,
    sample_id = names(cell_line_fsa_list[1])
  )

  testthat::expect_true(length(test_ladder_peaks) >= length(ladder_sizes))


  test_ladder_peaks_20 <- find_ladder_peaks(
    test_processed,
    n_reference_sizes = 20,
    minimum_ladder_signal = NULL,
    sample_id = names(cell_line_fsa_list[1])
  )

  testthat::expect_true(length(test_ladder_peaks_20) == 20)
})


test_that("iterative ladder", {
  scans_162 <- c(
    1548,
    1624,
    1771,
    1913,
    2124,
    2142,
    2200,
    2259,
    2503,
    2802,
    3131,
    3376,
    3439,
    3758,
    4050,
    4287,
    4336
  )



  ladder_sizes <- c(50, 75, 100, 139, 150, 160, 200, 300, 350, 400, 450, 490, 500)


  iteration_result <- ladder_iteration(ladder_sizes, scans_162,
    choose = 4,
    max_combinations = 2500000
  )

  expect_true(round(mean(iteration_result$scan), 3) == 2877.923)
})


test_that("find ladders", {

  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
  suppressWarnings(
    find_ladders(
      fsa_list,
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )



  testthat::expect_true(all(fsa_list[[1]]$ladder_df$scan == c(1540, 1618, 1766, 1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)))
})


test_that("find ladders scan subset", {

  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
  suppressWarnings(
    find_ladders(fsa_list,
      ladder_sizes = c(200, 250, 300, 340, 350, 400, 450),
      scan_subset = c(2400, 4250),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )



  testthat::expect_true(all(fsa_list[[1]]$ladder_df$scan == c(2502, 2802, 3131, 3376, 3438, 3756, 4046)))
})





test_that("ladder minium signal", {

  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())


    test_ladders <- find_ladders(fsa_list,
                                 ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                                 max_combinations = 2500000,
                                 ladder_selection_window = 8,
                                 show_progress_bar = FALSE,
                                 minimum_ladder_signal = 100
    )




  testthat::expect_true(all(fsa_list[[1]]$ladder_df$scan == c(1540, 1618, 1766, 1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)))
})



# test_that("fix ladders", {
#   file_list <- trace::cell_line_fsa_list

#   suppressWarnings(
#     test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
#       ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#       max_combinations = 2500000,
#       ladder_selection_window = 8,
#       show_progress_bar = FALSE
#     )
#   )


#   example_list <- list(
#     "20230413_B03.fsa" = data.frame(
#       size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#       scan = c(1555, 1633, 1783, 1827, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)
#     )
#   )

#   suppressMessages(
#     suppressWarnings(
#       test_ladders_fixed_manual <- fix_ladders_manual(
#         test_ladders,
#         example_list
#       )
#     )
#   )

#   suppressWarnings(
#     test_ladders_fixed <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")
#   )



#   testthat::expect_true(all(test_ladders_fixed$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1633, 1783, 1927, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
# })



test_that("fix ladders manual", {
 example_list <- list(
  "20230413_A07.fsa" = data.frame(
    size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    scan = c(1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)
  )
 )

  fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())

  find_ladders(fsa_list,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    max_combinations = 2500000,
    ladder_selection_window = 8,
    show_progress_bar = FALSE
  )

  suppressMessages(
    fix_ladders_manual(
      fsa_list,
      example_list
    )
  )

  expect_true(nrow(fsa_list[[1]]$ladder_df) == 13)
})

