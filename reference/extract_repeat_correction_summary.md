# Extract repeat correction summary

Extracts a table summarizing the model used to correct repeat length

## Usage

``` r
extract_repeat_correction_summary(fragments_list)
```

## Arguments

- fragments_list:

  A list of fragments class objects obtained from the
  [`call_repeats()`](https://zachariahmclean.github.io/trace/reference/call_repeats.md)
  function when the `correction = "repeat"` parameter is used.

## Value

A data.frame

## Details

For each of the samples used for repeat correction, this table pulls out
the modal repeat length called by the model (`allele_repeat`), how far
that sample is on average from the linear model in repeat units by
finding the average residuals (`avg_residual`), and the absolute value
of the `avg_residual` (`abs_avg_residual`)

## Examples

``` r
fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
# import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
test_fragments <- trace(
   fsa_list, 
   grouped = TRUE, 
   metadata_data.frame = metadata, 
   correction = "repeat",
   show_progress_bar = FALSE
)
#> Finding ladders
#> Finding fragments
#> Finding alleles
#> Calling repeats
#> Repeat correction model: 4 samples used to build model
#> Repeat correction model: 2.87 bp increase per repeat
#> Assigning index peaks

# finally extract repeat correction summary
extract_repeat_correction_summary(test_fragments)
#>                                   unique_id batch_run_id batch_sample_id
#> 20230413_H07.fsa           20230413_H07.fsa     20230414        S-21-211
#> 20230413_H08.fsa           20230413_H08.fsa     20230414        S-21-212
#> S-21-211_20220630.fsa S-21-211_20220630.fsa     20220630        S-21-211
#> S-21-212_20220630.fsa S-21-212_20220630.fsa     20220630        S-21-212
#>                       batch_sample_modal_repeat allele_repeat avg_residual
#> 20230413_H07.fsa                            120      119.9629   0.01585004
#> 20230413_H08.fsa                            122      122.0190  -0.01585004
#> S-21-211_20220630.fsa                       120      119.9821   0.01776718
#> S-21-212_20220630.fsa                       122      122.0069  -0.01522901
#>                       abs_avg_residual
#> 20230413_H07.fsa            0.01585004
#> 20230413_H08.fsa            0.01585004
#> S-21-211_20220630.fsa       0.01776718
#> S-21-212_20220630.fsa       0.01522901

```
