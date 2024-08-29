## code to prepare `cell_line_fsa_list` dataset goes here

cell_line_fsa_list <- read_fsa(list.files(full.names = TRUE, "data-raw/cell_line_fsa_list_files/"))

usethis::use_data(cell_line_fsa_list, overwrite = TRUE)
