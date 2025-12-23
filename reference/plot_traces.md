# Plot sample traces

Plot the raw trace data

## Usage

``` r
plot_traces(
  fragments_list,
  show_peaks = TRUE,
  n_facet_col = 1,
  sample_subset = NULL,
  xlim = NULL,
  ylim = NULL,
  x_axis = NULL,
  signal_color_threshold = 0.05
)
```

## Arguments

- fragments_list:

  A list of fragments or fragments objects containing fragment data.

- show_peaks:

  If peak data are available, TRUE will plot the peaks on top of the
  trace as dots.

- n_facet_col:

  A numeric value indicating the number of columns for faceting in the
  plot.

- sample_subset:

  A character vector of unique ids for a subset of samples to plot

- xlim:

  the x limits of the plot. A numeric vector of length two.

- ylim:

  the y limits of the plot. A numeric vector of length two.

- x_axis:

  A character indicating what should be plotted on the x-axis, chose
  between `size` or `repeats`. If neither is selected, an assumption is
  made based on if repeats have been called.

- signal_color_threshold:

  Threshold relative to tallest peak to color the dots (blue above,
  purple below).

## Value

plot traces from fragments object

## Details

A plot of the raw signal by bp size. Red vertical line indicates the
scan was flagged as off-scale. This is in any channel, so use your best
judgment to determine if it's from the sample or ladder channel.

If peaks are called, green is the tallest peak, blue is peaks above the
signal threshold (default 5%), purple is below the signal threshold. If
`force_whole_repeat_units` is used within
[`call_repeats()`](https://zachariahmclean.github.io/trace/reference/call_repeats.md),
the called repeat will be connected to the peak in the trace with a
horizontal dashed line.

The index peak will be plotted as a vertical dashed line when it has
been set using
[`assign_index_peaks()`](https://zachariahmclean.github.io/trace/reference/assign_index_peaks.md).

## Examples

``` r
fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
# import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
fragments_list <- trace(fsa_list)
#> Finding ladders
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding fragments
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding alleles
#> Calling repeats
#> Assigning index peaks

plot_traces(fragments_list, xlim = c(105, 150))
#> Warning: Error in plotting trace for fragment 20230413_A07.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_A08.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_A09.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_C01.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_C02.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_C03.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_D07.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_D08.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_D09.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_F01.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_F02.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_F03.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_G07.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_G08.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_G09.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_H07.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment 20230413_H08.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment S-21-211_20220630.fsa: figure margins too large
#> Warning: Error in plotting trace for fragment S-21-212_20220630.fsa: figure margins too large
```
