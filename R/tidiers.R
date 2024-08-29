# genemapper5_tidier ------------------------------------------------------

clean_genemapper5 <- function(df,
                              peak_size_col,
                              peak_height_col,
                              unique_id,
                              dye_col,
                              dye_channel,
                              allele_col) {
  df <- as.data.frame(df)
  # rename cols based on used supplied name
  names(df)[names(df) == peak_size_col] <- "size"
  names(df)[names(df) == peak_height_col] <- "height"
  names(df)[names(df) == unique_id] <- "unique_id"
  names(df)[names(df) == dye_col] <- "Dye.Sample.Peak"
  names(df)[names(df) == allele_col] <- "allele"

  # extract channel
  df$dye <- substr(df$Dye.Sample.Peak, 1, 1)

  # tidy dataframe names to snake case
  names(df) <- tolower(names(df))
  names(df) <- gsub("\\.", "_", names(df))

  # filter for only the dye of interest
  df2 <- df[which(df$dye == dye_channel), , drop = FALSE]

  return(df2)
}

# for generic peak table

clean_generic <- function(df,
                          peak_size_col,
                          peak_height_col,
                          unique_id) {
  df <- as.data.frame(df)
  # rename cols based on used supplied name
  names(df)[names(df) == peak_size_col] <- "size"
  names(df)[names(df) == peak_height_col] <- "height"
  names(df)[names(df) == unique_id] <- "unique_id"

  # tidy dataframe names to snake case
  names(df) <- tolower(names(df))
  names(df) <- gsub("\\.", "_", names(df))

  return(df)
}
