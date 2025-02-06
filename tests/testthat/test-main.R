test_that("main trace", {

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

  frag_list <- trace_main(fsa_list)


})


test_that("main fragments", {

  fsa_list <- peak_table_to_fragments(example_data, data_format = "genemapper5")

  frag_list <- trace_main(fsa_list)

})
