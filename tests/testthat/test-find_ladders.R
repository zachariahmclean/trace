testthat::test_that("find ladder peaks", {

  test_processed <- process_ladder_signal(cell_line_fsa_list[[1]]$fsa$Data$DATA.105,
    scans = 0:(length(cell_line_fsa_list[[1]]$fsa$Data$DATA.105) - 1),
    ladder_start_scan = 1000,
    smoothing_window = 21
  )


  ladder_sizes <- c(50, 75, 100, 139, 150, 160, 200, 300, 350, 400, 450, 490, 500)


  test_ladder_peaks <- find_ladder_peaks(
    test_processed,
    length(ladder_sizes),
    minimum_peak_signal = NULL,
    sample_id = names(file_list[1])
  )

  testthat::expect_true(length(test_ladder_peaks) >= length(ladder_sizes))


  test_ladder_peaks_32 <- find_ladder_peaks(
    test_processed,
    n_reference_sizes = 32,
    minimum_peak_signal = NULL,
    sample_id = names(file_list[1])
  )

  testthat::expect_true(length(test_ladder_peaks_32) == 32)
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



# test_that("fit ladder", {
#   file_list <- trace::cell_line_fsa_list

#   test_ladder_signal <- file_list[[1]]$Data$DATA.105
#   test_scans <- 0:(length(file_list[[1]]$Data$DATA.105) - 1)



#   test_fit <- fit_ladder(
#     ladder = test_ladder_signal,
#     scans = test_scans,
#     ladder_start_scan = NULL,
#     ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#     smoothing_window = 21,
#     minimum_peak_signal = NULL,
#     zero_floor = FALSE,
#     max_combinations = 2500000,
#     ladder_selection_window = 5
#   )


#   mod <- lm(scan ~ size, data = test_fit)

#   testthat::expect_true(round(mod$coefficients[[1]], 3) == 1347.69)
#   testthat::expect_true(round(mod$coefficients[[2]], 5) == 6.25118)
# })


# test_that("local southern", {
#   file_list <- trace::cell_line_fsa_list

#   test_ladder_signal <- file_list[[1]]$Data$DATA.105
#   test_scans <- 0:(length(file_list[[1]]$Data$DATA.105) - 1)



#   test_fit <- fit_ladder(
#     ladder = test_ladder_signal,
#     scans = test_scans,
#     ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#     ladder_start_scan = NULL,
#     smoothing_window = 21,
#     minimum_peak_signal = NULL,
#     zero_floor = FALSE,
#     max_combinations = 2500000,
#     ladder_selection_window = 5
#   )

#   # ladder all as one lm


#   mod <- lm(size ~ scan, data = test_fit)
#   single_rsqd <- summary(mod)$r.squared




#   mod_fit <- local_southern_fit(test_fit$scan, test_fit$size)

#   multi_rsqd <- sapply(mod_fit, function(x) summary(x$mod)$r.squared)

#   rsq_diff <- sum(multi_rsqd - single_rsqd)


#   predicted_sizes <- local_southern_predict(mod_fit, test_scans)

#   # plot(predicted_sizes, file_list[[1]]$Data$DATA.1)


#   testthat::expect_true(round(rsq_diff, 5) == -0.00343)
# })



test_that("find ladders", {

  fsa_list <- lapply(cell_line_fsa_list["20230413_B03.fsa"], function(x) x$clone())
  suppressWarnings(
    find_ladders(
      fsa_list,
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )



  testthat::expect_true(all(fsa_list$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1633, 1783, 1927, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
})


test_that("find ladders scan subset", {

  fsa_list <- lapply(cell_line_fsa_list["20230413_B03.fsa"], function(x) x$clone())
  suppressWarnings(
    find_ladders(fsa_list,
      ladder_sizes = c(200, 250, 300, 340, 350, 400, 450),
      scan_subset = c(2400, 4250),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )



  testthat::expect_true(all(fsa_list$`20230413_B03.fsa`$ladder_df$scan == c(2525, 2828, 3161, 3408, 3470, 3792, 4085)))
})





test_that("ladder minium height", {

  fsa_list <- lapply(cell_line_fsa_list["20230413_B03.fsa"], function(x) x$clone())


    test_ladders <- find_ladders(fsa_list,
                                 ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                                 max_combinations = 2500000,
                                 ladder_selection_window = 8,
                                 show_progress_bar = FALSE,
                                 minimum_peak_signal = 100
    )




  testthat::expect_true(all(test_ladders$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1633, 1783, 1927, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
})


test_that("ladder zero baseline", {
  fsa_list <- lapply(cell_line_fsa_list["20230413_B03.fsa"], function(x) x$clone())



  test_ladders <- find_ladders(fsa_list,
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               max_combinations = 2500000,
                               ladder_selection_window = 8,
                               show_progress_bar = FALSE,
                               zero_floor = TRUE
  )




  testthat::expect_true(all(test_ladders$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1633, 1783, 1927, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
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
    "20230413_A01.fsa" = data.frame(
      size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1593, 1671, 1825, 1971, 2208, 2269, 2329, 2581, 2888, 3228, 3479, 3543, 3872, 4170, 4412, 4460)
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

  expect_true(nrow(fsa_list[[1]]$ladder_df) == 16)
})

