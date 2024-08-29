
# remember that this test works on the installed version of the code, not devtools::load_all()!


test_that("{shinytest2} recording: fix_ladder-checkbox", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()


  file_list <- lapply(trace::cell_line_fsa_list[c("20230413_A01.fsa", "20230413_B03.fsa")], function(x) x$clone())


  test_ladders <- find_ladders(file_list,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    max_combinations = 2500000,
    ladder_selection_window = 8,
    show_progress_bar = FALSE
  )

  # sapply(test_ladders[[2]]$local_southern_mod, function(x)  summary(x$mod)$r.squared)

  example_list <- list(
    "20230413_B03.fsa" = data.frame(
      size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1555, 1633, 1783, 1827, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)
    )
  )

  suppressMessages(
    suppressWarnings(
      test_ladders_fixed_manual <- fix_ladders_manual(
        test_ladders,
        example_list
      )
    )
  )

  # sapply(test_ladders_fixed_manual[[2]]$local_southern_mod, function(x)  summary(x$mod)$r.squared)


  # to generate the values below I ran a test like below and copied the values
  # shinytest2::record_test(fix_ladders_interactive(test_ladders_fixed_manual))


  suppressMessages(
    shiny_app <- fix_ladders_interactive(test_ladders_fixed_manual)
  )


  app <- shinytest2::AppDriver$new(shiny_app, height = 945, width = 1619)
  # Update output value
  app$set_inputs(`plotly_afterplot-A` = "\"plot_module-plot\"", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_relayout-A` = "{\"width\":1037.984375,\"height\":400}", allow_no_input_binding_ = TRUE, priority_ = "event")
  # Update output value
  app$set_inputs(`sample_selection-warning_checkbox` = TRUE)
  # Update output value
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":1497,\"x\":1497,\"y\":4088}]", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":1506,\"x\":1506,\"y\":3913}]", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_afterplot-A` = "\"plot_module-plot\"", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_relayout-A` = "{\"shapes[3].x0\":1913.172314049587,\"shapes[3].x1\":1913.172314049587,\"shapes[3].y0\":0.050000000000000044,\"shapes[3].y1\":0.44999999999999996}", allow_no_input_binding_ = TRUE, priority_ = "event")

  rsq_table_html <- app$get_values()$output$`rsq_table-rsq_table`

  # Extract numbers with decimal points using regular expressions
  matches <- gregexpr("\\b\\d+(\\.\\d+)?\\b(?= </td> </tr>\\n)", rsq_table_html, perl = TRUE)

  # Extract matched numbers
  numbers <- regmatches(rsq_table_html, matches)[[1]]
  numbers <- round(as.numeric(numbers), 4)


  # I got these numbers by changing to sample 20230413_B03.fsa and moving 100 back to the right place, then round as above and remove the 1.0000
  expect_identical(numbers, c(0.9986, 0.9995, 0.9994, 0.9995, 0.9990, 1.0000, 1.0000, 0.9993, 0.9995, 1.0000, 1.000, 0.9993, 1.000, 0.9990))
})
