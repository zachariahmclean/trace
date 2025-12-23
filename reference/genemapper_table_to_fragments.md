# Convert Genemapper Peak Table to fragments class

This function converts a genemapper peak table data frame into a list of
fragments objects.

## Usage

``` r
genemapper_table_to_fragments(
  df,
  dye_channel,
  min_size_bp = 200,
  max_size_bp = 1000
)
```

## Arguments

- df:

  A data frame from genemapper containing the peak data.

- dye_channel:

  A character string indicating the Genemapper channel to extract data
  from. This is used to filter 'Dye.Sample.Peak' for the channel
  containing the data. For example, 6-FAM is often "B" while ladder is
  "O".

- min_size_bp:

  Numeric value indicating the minimum size of the peak table to import.

- max_size_bp:

  Numeric value indicating the maximum size of the peak table to import.

## Value

A list of fragments objects.

## Details

This function takes a peak table data frame (eg. Genemapper 5 output)
and converts it into a list of fragment objects. It uses the
"Sample.File.Name" as the unique id, so make sure that each file as a
unique name. Column names should contain: "Dye.Sample.Peak",
"Sample.File.Name", "Allele", "Size", "Height".

## See also

[`repeat_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/repeat_table_to_fragments.md),
[`size_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/size_table_to_fragments.md),
[`read_fsa()`](https://zachariahmclean.github.io/trace/reference/read_fsa.md)

## Examples

``` r
gm_raw <- trace::example_data

test_fragments <- genemapper_table_to_fragments(
  gm_raw,
  dye_channel = "B",
  min_size_bp = 400
)
```
