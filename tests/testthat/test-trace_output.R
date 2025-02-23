
test_that("trace_output error", {
  test_output <- trace_output$new("test")

  test_output$set_status("error", " error")
  print_test <- tryCatch(print(test_output), error = function(e) e )
 
  expect_true("error" %in% class(print_test))

})

 
 
test_that("trace_output warning", {
  test_output <- trace_output$new("test")

  test_output$set_status("warning", "hello world")

  expect_true(test_output$status == "warning")

})
