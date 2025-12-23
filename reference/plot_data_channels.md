# plot_data_channels

Plot the raw data from the fsa file

## Usage

``` r
plot_data_channels(fragments_list, sample_subset = NULL, n_facet_col = 1)
```

## Arguments

- fragments_list:

  A list of fragments objects.

- sample_subset:

  A character vector of unique ids for a subset of samples to plot

- n_facet_col:

  A numeric value indicating the number of columns for faceting in the
  plot.

## Value

a plot of the raw data channels

## Details

A plot of the raw data channels in the fsa file.

These plots are made using base R plotting. Sometimes these fail to
render in the viewing panes of IDEs (eg you get the error 'Error in
[`plot.new()`](https://rdrr.io/r/graphics/frame.html): figure margins
too large)'. If this happens, try saving the plot as a pdf using
traditional approaches (see grDevices::pdf). To get it to render in the
IDE pane, trying matching `n_facet_col` to the number of samples you're
attempting to plot, or using `sample_subset` to limit it to a single
sample.

## Examples

``` r
plot_data_channels(cell_line_fsa_list[1])

```
