# Calculate Repeat Instability Metrics

This function computes instability metrics from a list of fragments data
objects.

## Usage

``` r
calculate_instability_metrics(
  fragments_list,
  peak_threshold = 0.05,
  window_around_index_peak = c(NA_real_, NA_real_),
  percentile_range = c(0.5, 0.75, 0.9, 0.95),
  repeat_range = c(2, 5, 10, 20),
  index_modal_signal_threshold = NA_real_,
  index_signal_sum_threshold = NA_real_
)
```

## Arguments

- fragments_list:

  A list of "fragments" objects representing fragment data.

- peak_threshold:

  A single numeric value between 0 and 1 for the threshold of peak
  signals to be considered in the calculations, relative to the modal
  peak signal of the expanded allele.

- window_around_index_peak:

  A numeric vector (length 2) defining the range around the index peak.
  First number specifies repeats before the index peak, second after.
  For example, `c(-5, 40)` around an index peak of 100 would analyze
  repeats 95 to 140. The sign of the numbers does not matter (The
  absolute value is found).

- percentile_range:

  A numeric vector of percentiles to compute (e.g., c(0.5, 0.75, 0.9,
  0.95)).

- repeat_range:

  A numeric vector specifying ranges of repeats for the inverse quantile
  computation.

- index_modal_signal_threshold:

  A single numeric value for the minimum signal of the modal peak for
  the index samples (basically a quality control for the samples used to
  set the index peak or to calculate average_repeat_change or
  instability_index_change). This is only relevant when grouped = TRUE
  for the index peak assignment.

- index_signal_sum_threshold:

  A single numeric value for the minimum sum of all peaks for each index
  sample (basically a quality control for the samples used to set the
  index peak or to calculate average_repeat_change or
  instability_index_change). This is only relevant when grouped = TRUE
  for the index peak assignment.

## Value

A data.frame with calculated instability metrics for each sample.

## Details

Each of the columns in the supplied dataframe are explained below:

### General Information

- `unique_id`: A unique identifier for the sample (usually the fsa file
  name).

### Quality Control

- `QC_comments`: Quality control comments.

- `QC_modal_peak_signal`: Quality control status based on the modal peak
  signal (Low \< 500, very low \< 100).

- `QC_peak_number`: Quality control status based on the number of peaks
  (Low \< 20, very low \< 10).

- `QC_off_scale`: Quality control comments for off-scale peaks.
  Potential peaks that are off-scale are given. However, a caveat is
  that this could be from any of the channels (ie it could be from the
  ladder channel but is the same scan as the given repeat).

### settings used

- `peak_threshold`: THe peak_threshold parameter used.

- `lower_repeat_threshold`: The lower repeat limit based of the index
  repeat of each sample.

- `upper_repeat_threshold`: The upper repeat limit based of the index
  repeat of each sample.

- `index_modal_signal_threshold`: The index_modal_signal_threshold
  parameter used.

- `index_signal_sum_threshold`: The index_signal_sum_threshold parameter
  used.

### General sample metrics

- `modal_peak_repeat`: The repeat size of the modal peak.

- `modal_peak_signal`: The signal of the modal peak.

- `index_peak_repeat`: The repeat size of the index peak (the repeat
  value closest to the modal peak of the index sample).

- `index_weighted_mean_repeat`: The weighted mean repeat size (weighted
  on the signal of the peaks) of the index sample.

- `n_peaks_total`: The total number of peaks in the repeat table.

- `n_peaks_analysis_subset`: The number of peaks in the analysis subset.

- `n_peaks_analysis_subset_expansions`: The number of expansion peaks in
  the analysis subset.

- `min_repeat`: The minimum repeat size in the analysis subset.

- `max_repeat`: The maximum repeat size in the analysis subset.

- `mean_repeat`: The mean repeat size in the analysis subset.

- `weighted_mean_repeat`: The weighted mean repeat size (weight on peak
  signal) in the analysis subset.

- `median_repeat`: The median repeat size in the analysis subset.

- `max_signal`: The maximum peak signal in the analysis subset.

- `sum_signal`: The sum of the peak signal in the analysis subset.

- `max_delta_neg`: The maximum negative delta to the index peak.

- `max_delta_pos`: The maximum positive delta to the index peak.

- `skewness`: The skewness of the repeat size distribution.

- `kurtosis`: The kurtosis of the repeat size distribution.

### Repeat instability metrics

- `modal_repeat_change`: The difference between the modal repeat and the
  index repeat.

- `average_repeat_change`: The weighted mean of the sample (weighted by
  peak signal) subtracted by the weighted mean repeat of the index
  sample(s).

- `instability_index_change`: The instability index of the sample
  subtracted by the instability index of the index sample(s). This will
  be very similar to the average_repeat_change, with the key difference
  of instability_index_change being that it is an internally calculated
  metric for each sample, and therefore the random slight fluctuations
  of bp size (or systematic if across plates for example) will be
  removed. However, it requires the index peak to be correctly set for
  each sample, and if set incorrectly, can produce large arbitrary
  differences.

- `instability_index`: The instability index based on peak signal and
  distance to the index peak. (See Lee et al., 2010,
  [doi:10.1186/1752-0509-4-29](https://doi.org/10.1186/1752-0509-4-29)
  ).

- `instability_index_abs`: The absolute instability index. The absolute
  value is taken for the "Change from the main allele".

- `expansion_index`: The instability index for expansion peaks only.

- `contraction_index`: The instability index for contraction peaks only.

- `expansion_ratio`: The ratio of expansion peaks' signals to the main
  peak signal. Also known as "peak proportional sum" (See Genetic
  Modifiers of Huntington’s Disease (GeM-HD) Consortium, 2019,
  [doi:10.1016/j.cell.2019.06.036](https://doi.org/10.1016/j.cell.2019.06.036)
  ).

- `contraction_ratio`: The ratio of contraction peaks' signals to the
  main peak signal.

- `expansion_percentile_*`: The repeat size at specified percentiles of
  the cumulative distribution of expansion peaks.

- `expansion_percentile_for_repeat_*`: The percentile rank of specified
  repeat sizes in the distribution of expansion peaks.

## Examples

``` r
fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
# import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
test_fragments <- trace(fsa_list, grouped = TRUE, metadata_data.frame = metadata)
#> Finding ladders
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding fragments
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding alleles
#> Calling repeats
#> Assigning index peaks

test_metrics_grouped <- calculate_instability_metrics(
  fragments_list = test_fragments,
  peak_threshold = 0.05,
  window_around_index_peak = c(-40, 40)
)
```
