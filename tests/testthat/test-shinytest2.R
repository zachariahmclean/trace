
# remember that this test works on the installed version of the code, not devtools::load_all()!


test_that("{shinytest2} recording: fix_ladder-checkbox", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()
  config <- load_config()
  fsa_list <- lapply(trace::cell_line_fsa_list[1:2], function(x) x$clone())


  find_ladders(fsa_list,
    config,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    max_combinations = 2500000,
    ladder_selection_window = 8,
    show_progress_bar = FALSE
  )

  example_list <- list(
    "20230413_A08.fsa" = data.frame(
      size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1544, 1621, 1850, 1912, 2143, 2201, 2261, 2506, 2805, 3135, 3380, 3442, 3760, 4050, 4284, 4332)
    )
  )

  suppressMessages(
    suppressWarnings(
      fix_ladders_manual(
        fsa_list,
        example_list
      )
    )
  )


  # to generate the values below I ran a test like below and copied the values
  # shinytest2::record_test(fix_ladders_interactive(fsa_list))


  suppressMessages(
    shiny_app <- fix_ladders_interactive(fsa_list)
  )

  
  app <- shinytest2::AppDriver$new(shiny_app, width = 1235, height = 826)
  # Update output value

  app$set_inputs(`plotly_afterplot-A` = "\"plot_module-plot\"", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_relayout-A` = "{\"width\":782.2625122070312,\"height\":400}", allow_no_input_binding_ = TRUE, priority_ = "event")
  # Update output value
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":4946,\"x\":4946,\"y\":22}]", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`sample_selection-unique_id_selection` = "20230413_A08.fsa", wait_ = FALSE)


  rsq_table_html <- app$get_values()$output$`rsq_table-rsq_table`

  # Extract numbers with decimal points using regular expressions
  matches <- gregexpr("\\b\\d+(\\.\\d+)?\\b(?= </td> </tr>\\n)", rsq_table_html, perl = TRUE)

  # Extract matched numbers
  numbers <- regmatches(rsq_table_html, matches)[[1]]
  numbers <- round(as.numeric(numbers), 4)


  # I got these numbers by changing to second sample and moving 75 back to the right place, then round as above and remove the 1.0000
  expect_identical(numbers, c(0.9816, 0.9011, 0.9621, 0.9996, 0.9986, 1.0000, 1.0000, 0.9992, 0.9996, 1.0000, 1.0000, 0.9993, 1.0000, 0.9991))




  # Update output value
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":2143,\"x\":2143,\"y\":1877}]", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":1912,\"x\":1912,\"y\":1700}]", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_afterplot-A` = "\"plot_module-plot\"", allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_relayout-A` = "{\"shapes[2].x0\":1774.6904494382022,\"shapes[2].x1\":1774.6904494382022,\"shapes[2].y0\":0.050000000000000044,\"shapes[2].y1\":0.44999999999999996}", allow_no_input_binding_ = TRUE, priority_ = "event")

  rsq_table_html <- app$get_values()$output$`rsq_table-rsq_table`

  # Extract numbers with decimal points using regular expressions
  matches <- gregexpr("\\b\\d+(\\.\\d+)?\\b(?= </td> </tr>\\n)", rsq_table_html, perl = TRUE)

  # Extract matched numbers
  numbers <- regmatches(rsq_table_html, matches)[[1]]
  numbers <- round(as.numeric(numbers), 4)


  # I got these numbers by changing to second sample and moving 75 back to the right place, then round as above and remove the 1.0000
  expect_identical(numbers, c(0.9986, 0.9999, 0.9999, 0.9996, 0.9986, 1.0000, 1.0000, 0.9992, 0.9996, 1.0000, 1.0000, 0.9993, 1.0000, 0.9991))
})
