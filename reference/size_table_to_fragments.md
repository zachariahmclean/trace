# Convert Size Table to Fragments

This function converts a size table data frame into a list of fragments.
class.

## Usage

``` r
size_table_to_fragments(df, min_size_bp = 200, max_size_bp = 1000)
```

## Arguments

- df:

  A data frame containing the size data with the columns "unique_id",
  "size", "signal".

- min_size_bp:

  Numeric value indicating the minimum size of the peak table to import.

- max_size_bp:

  Numeric value indicating the maximum size of the peak table to import.

## Value

A list of fragments objects.

## Details

This function takes a size table data frame and converts it into a list
of fragments objects. The column names must be "unique_id" (unique
sample id), "size" (base pair size), and "signal" (the signal associated
with fragment). The dataframe should be long (eg bind rows if the data
are separate).

## See also

[`repeat_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/repeat_table_to_fragments.md),
`size_table_to_fragments()`,
[`read_fsa()`](https://zachariahmclean.github.io/trace/reference/read_fsa.md)

## Examples

``` r
size_table <- trace::example_data
colnames(size_table)[c(2, 5, 6)] <- c("unique_id", "size", "signal")
fragments_list <- size_table_to_fragments(size_table)
```
