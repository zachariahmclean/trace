# Convert Repeat Table to Repeats Fragments

This function converts a repeat table data frame into a list of
fragments. class.

## Usage

``` r
repeat_table_to_fragments(df, min_repeat = 0, max_repeat = 1000)
```

## Arguments

- df:

  A data frame containing the repeat data with the columns "unique_id",
  "repeats", "signal".

- min_repeat:

  minimum repeat size

- max_repeat:

  maximum repeat size

## Value

A list of fragments objects.

## Details

This function takes a repeat table data frame and converts it into a
list of fragments objects. The column names must be "unique_id" (unique
sample id), "repeats" (number of repeats), and "signal" (the signal
associated with the number of repeats. for example a frequency count of
reads.). The dataframe should be long (eg bind rows if the data are
separate).

## Examples

``` r
repeat_table <- trace::example_data_repeat_table
test_fragments <- repeat_table_to_fragments(repeat_table)
```
