# init --------------------------------------------------------------------

testthat::test_that("fragments class initialization", {
  test_fragments <- fragments$new(
    unique_id = "test"
  )

  testthat::expect_equal(
    class(test_fragments)[1], "fragments"
  )
})
