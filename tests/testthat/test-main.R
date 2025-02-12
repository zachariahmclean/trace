test_that("main trace", {

  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

  frag_list <- trace_main(fsa_list)

  frag_list <- trace_main(fsa_list, grouped = TRUE, metadata_data.frame = metadata)


  # ladder fixing
  fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

  example_list <- list(
    "20230413_A07.fsa" = data.frame(
      size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)
    )
   )

   frag_list <- trace_main(fsa_list, ladder_df_list = example_list)

  
   

})


test_that("main fragments", {

  fsa_list <- peak_table_to_fragments(example_data, data_format = "genemapper5")

  frag_list <- trace_main(fsa_list)

})
