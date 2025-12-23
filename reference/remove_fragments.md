# Remove Samples from List

A convenient function to remove specific samples from a list of
fragments.

## Usage

``` r
remove_fragments(fragments_list, samples_to_remove)
```

## Arguments

- fragments_list:

  A list of fragments objects containing fragment data.

- samples_to_remove:

  A character vector containing the unique IDs of the samples to be
  removed.

## Value

A modified list of fragments with the specified samples removed.

## Examples

``` r
gm_raw <- trace::example_data

test_fragments <- genemapper_table_to_fragments(
  gm_raw,
  dye_channel = "B",
  min_size_bp = 300
)

all_fragment_names <- names(test_fragments)

# pull out unique ids of samples to remove
samples_to_remove <- all_fragment_names[c(1, 5, 10)]

samples_removed <- remove_fragments(test_fragments, samples_to_remove)
```
