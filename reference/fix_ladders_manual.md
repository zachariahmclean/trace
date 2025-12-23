# Fix ladders manually

Manually assign the ladder peaks for samples in a fragments_list

## Usage

``` r
fix_ladders_manual(
  fragments_list,
  ladder_df_list,
  warning_rsq_threshold = 0.998
)
```

## Arguments

- fragments_list:

  list of fragments objects

- ladder_df_list:

  a list of dataframes, with the names being the unique id and the value
  being a dataframe. The dataframe has two columns, size (indicating the
  bp of the standard) and scan (the scan value of the ladder peak). It's
  critical that the element name in the list is the unique id of the
  sample.

- warning_rsq_threshold:

  The value for which this function will warn you when parts of the
  ladder have R-squared values below the specified threshold.

## Value

This function modifies list of fragments objects in place with the
selected ladders fixed.

## Details

This function returns a fragments list the same length as was supplied.
It goes through each sample and either just returns the same fragments
if the unique id doesn't match the samples that need the ladder fixed,
or if it is one to fix, it will use the supplied dataframe in the
ladder_df_list as the ladder. It then reruns the bp sizing methods on
those samples.

This is best used with
[`fix_ladders_interactive()`](https://zachariahmclean.github.io/trace/reference/fix_ladders_interactive.md)
that can generate a `ladder_df_list`.

## Examples

``` r
config <- load_config()
fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())

trace:::find_ladders(fsa_list, config, show_progress_bar = FALSE)

# first manually determine the real ladder peaks using your judgment
# the raw ladder signal can be extracted
raw_ladder <- fsa_list[1]$raw_ladder

# or we can look at the "trace_bp_df" to see a dataframe that includes the scan and ladder signal
raw_ladder_df <- fsa_list[[1]]$trace_bp_df[, c("unique_id", "scan", "ladder_signal")]
plot(raw_ladder_df$scan, raw_ladder_df$ladder_signal)


# once you have figured what sizes align with which peak, make a dataframe. The
# fix_ladders_manual() function takes a list as an input so that multiple ladders
# can be fixed. Each sample would have the the list element name as it's unique id.

example_list <- list(
 "20230413_A07.fsa" = data.frame(
   size = c(100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
   scan = c(1909, 2139, 2198, 2257, 2502, 2802, 3131, 3376, 3438, 3756, 4046, 4280, 4328)
 )
)

trace:::fix_ladders_manual(
  fsa_list,
  example_list
)
#> Fixing ladder for 20230413_A07.fsa
```
