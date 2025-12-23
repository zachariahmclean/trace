# Main function for sample processing

The main function for the trace package that handles processing of
samples through the pipeline ready for the calculation of repeat
instability metrics.

## Usage

``` r
trace(
  fragments_list,
  metadata_data.frame = NULL,
  index_override_dataframe = NULL,
  ladder_df_list = NULL,
  config_file = NULL,
  ...
)
```

## Arguments

- fragments_list:

  A list of fragments objects containing fragment data, generated with
  either
  [`read_fsa()`](https://zachariahmclean.github.io/trace/reference/read_fsa.md),
  [`size_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/size_table_to_fragments.md),
  [`genemapper_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/genemapper_table_to_fragments.md),
  or
  [`repeat_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/repeat_table_to_fragments.md).

- metadata_data.frame:

  metadata passed to
  [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md)
  for grouping samples for metrics calculations or batch correction.

- index_override_dataframe:

  A data.frame to manually set index peaks. Column 1: unique sample IDs,
  Column 2: desired index peaks (the order of the columns is important
  since the information is pulled by column position rather than column
  name). Closest peak in each sample is selected so the number needs to
  just be approximate. Default: `NULL`. See
  [`assign_index_peaks()`](https://zachariahmclean.github.io/trace/reference/assign_index_peaks.md).

- ladder_df_list:

  A list of dataframes, with the names being the unique id and the value
  being a dataframe. The dataframe has two columns, size (indicating the
  bp of the standard) and scan (the scan value of the ladder peak). It's
  critical that the element name in the list is the unique id of the
  sample. Either manually figure out what scan the ladder peaks should
  be and generate the list, or use
  [`fix_ladders_interactive()`](https://zachariahmclean.github.io/trace/reference/fix_ladders_interactive.md)
  to interactively generate the ladder_df_list.

- config_file:

  The file path to a YAML file containing a list of parameters. This
  provides a central place to adjust parameters for the pipeline. Use
  the following command to make a copy of the default YAML file:
  `file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.
  The YAML file does not have to contain a complete list of parameters.
  Parameters directly passed to this function via ... will overwrite
  values in the `config_file`.

- ...:

  additional parameters from any of the functions in the pipeline
  detailed below may be passed to this function.These parameters are
  grouped by functionality:

  **Ladder-related parameters:**

  - `ladder_channel`: string, which channel in the fsa file contains the
    ladder signal. Default: `"DATA.105"`.

  - `signal_channel`: string, which channel in the fsa file contains the
    data signal. Default: `"DATA.1"`.

  - `ladder_sizes`: numeric vector, bp sizes of ladder used in fragment
    analysis. Default:
    `c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)`.

  - `ladder_start_scan`: single numeric indicating the scan number to
    start looking for ladder peaks (only required when ladder signal
    does not have large spike at start). Usually this can be
    automatically found (when set to NA) through the detection of the
    large spike at the start of the signal. Default: `NA`.

  - `minimum_ladder_signal`: single numeric for minimum signal of peak
    from the raw signal. If not set, the ladder will be fit to a set of
    the tallest peaks in the ladder channel. Default: `NA`.

  - `ladder_assign_left_to_right`: single logical for if the ladder
    should be assigned from the smallest base pair size to largest
    (TRUE), or if the order should be reversed and assigned from largest
    to smallest (FALSE), which can be helpful since the end often has
    cleaner signal than the start. Default: `TRUE`.

  - `ladder_selection_window`: single numeric for the ladder assigning
    algorithm. We iterate through the scans in blocks and test their
    linear fit (We can assume that the ladder is linear over a short
    distance). This value defines how large that block of peaks should
    be. Larger values should be better because the fit is tested in
    greater context, but larger numbers will make the fit increasingly
    slower. Default: `5`.

  - `ladder_top_n_branching`: single numeric. The ladder assigning
    algorithm branches as it tests the various combinations. This value
    defines how many branches should be created. If the correct
    combination is not found, you could try increasing this value, but
    it will make it increasingly slower. Default: `5`.

  - `ladder_branching_r_squared_threshold`: single numeric. The branches
    of the ladder assigning algorithm are pruned by R-squared values
    above this threshold to discard fits that are not promising. If the
    correct combination is not found, you could try decreasing this
    value, but it will make it increasingly slower. Default: `0.99`.

  - `min_scan`: single numeric indicating the lower scan limit to filter
    out scans below. Default: `NA`.

  - `max_scan`: single numeric indicating the upper scan limit to filter
    out scans above. Default: `NA`.

  - `max_combinations`: single numeric indicating what is the maximum
    number of ladder combinations that should be tested. Default:
    `2500000`.

  - `warning_rsq_threshold`: single numeric for the value for which this
    function will warn you when parts of the ladder have R-squared
    values below the specified threshold. Default: `0.998`.

  - `show_progress_bar`: single logical for showing progress bar.
    Default: `TRUE`.

  **Peak-finding parameters:**

  - `smoothing_window`: Single numeric value for signal smoothing window
    size passed to pracma::savgol(). Default: `21`.

  - `minimum_peak_signal`: Single numeric value for the minimum peak
    signal for a valid peak. To have no minimum signal set as "-Inf".
    Default: `20`

  - `min_bp_size`: Single numeric value for minimum bp size of peaks to
    consider. Default: `200`.

  - `max_bp_size`: Single numeric value for maximum bp size of peaks to
    consider. Default: `1000`.

  - `peak_scan_ramp`: Single numeric value to indicate how many scans
    (increasing in signal) should be either side of the peak maxima.
    Default: `5`.

  **Allele-calling parameters:**

  - `number_of_alleles`: Number of alleles to be returned for each
    fragment. Must either be 1 or 2. Default: `1`.

  - `peak_region_size_gap_threshold`: Gap threshold for identifying peak
    regions. Default: `6`.

  - `peak_region_signal_threshold_multiplier`: Multiplier for the peak
    signal threshold. Default: `1`.

  **Repeat-calling parameters:**

  - `assay_size_without_repeat`: An integer specifying the assay size
    without repeat for repeat calling. Default: `87`.

  - `repeat_size`: An integer specifying the repeat size for repeat
    calling. Default: `3`.

  - `correction`: A character vector of either "batch" to carry out a
    batch correction from common samples across runs (known repeat
    length not required), or "repeat" to use samples with validated
    modal repeat lengths to correct the repeat length. Default:
    `"none"`.

  - `force_whole_repeat_units`: A logical value specifying if the peaks
    should be forced to be whole repeat units apart. Default: `FALSE`.

  - `force_repeat_pattern`: A logical value specifying if the peaks
    should be re-called to fit the specific repeat unit pattern.
    Default: `FALSE`.

  - `force_repeat_pattern_size_period`: A numeric value to set the peak
    periodicity bp size. Default: `2.79`.

  - `force_repeat_pattern_size_window`: A numeric value for the size
    window when assigning the peak. Default: `0.5`.

  **Index peak assignment parameters:**

  - `grouped`: Logical value indicating whether samples should be
    grouped to share a common index peak. Default: `FALSE`.

## Value

A list of fragments objects ready for calculation of instability metrics
using
[`calculate_instability_metrics()`](https://zachariahmclean.github.io/trace/reference/calculate_instability_metrics.md)

## Details

This function processes samples through the full pipeline, applying the
library of functions within this package. Parameters can be adjusted
either by passing them directly to this function or by editing a
configuration file and supplying it to `config_file`. Below is a
breakdown of the pipeline stages and their key functionalities:

**Ladder-related parameters:** The ladder is used to calibrate the base
pair (bp) sizes of fragments. The ladder peaks are identified in the
ladder channel using
[`find_ladders`](https://zachariahmclean.github.io/trace/reference/find_ladders.md),
and a generalized additive model (GAM) with cubic regression splines is
used to fit the relationship between scan numbers and bp sizes. This
allows for accurate interpolation of bp sizes for all scans. Manual
inspection of ladder assignments is recommended to ensure correctness.
If ladders are broken, first try and adjust parameters. In a last
resort, ladders may be manually set by using the `ladder_df_list`
parameter. To generate this list, use the helper shiny app
[`fix_ladders_interactive`](https://zachariahmclean.github.io/trace/reference/fix_ladders_interactive.md).

**Peak-finding parameters:** Fragment peaks are identified in the
continuous trace data using
[`find_fragments`](https://zachariahmclean.github.io/trace/reference/find_fragments.md).
The signal is smoothed using a Savitzky-Golay filter, and peaks are
detected based on their signal intensity and bp size. Parameters such as
`min_bp_size` and `max_bp_size` allow filtering peaks outside the
desired range, while `peak_scan_ramp` controls the number of scans
around the peak maxima.

**Allele-calling parameters:** The main alleles within each fragment are
identified by clustering peaks into "peak regions" using
[`find_alleles`](https://zachariahmclean.github.io/trace/reference/find_alleles.md).
The tallest peak in each region is selected as the main allele. If
`number_of_alleles` is set to 2, the two tallest peaks in their
respective regions are selected, with the larger repeat size assigned as
the main allele.

**Repeat-calling parameters:** Repeat lengths are calculated using
[`call_repeats`](https://zachariahmclean.github.io/trace/reference/call_repeats.md),
based on the assay size without the repeat and the specified repeat
size. Options include batch correction to account for run-to-run
variability and repeat correction using samples with known modal repeat
lengths. The `force_whole_repeat_units` option ensures repeat lengths
are whole numbers, while `force_repeat_pattern` re-calls peaks to fit
the expected repeat pattern.

**Index peak assignment parameters:** The index peak is the reference
repeat length used for instability metrics calculations. It is typically
the inherited repeat length or the modal repeat length at a starting
time point. Samples can be grouped to share a common index peak using
[`assign_index_peaks`](https://zachariahmclean.github.io/trace/reference/assign_index_peaks.md),
or the index peak can be manually overridden using
`index_override_dataframe`.

**Metadata:** Metadata is used to group samples for metrics
calculations, batch correction, or repeat length validation. Use
[`add_metadata`](https://zachariahmclean.github.io/trace/reference/add_metadata.md)
to add metadata to the fragments list. Required (even if the column is
left blank) metadata fields include:

- `batch_run_id`: Groups samples by fragment analysis run.

- `batch_sample_id`: Links samples across batches for batch or repeat
  correction.

- `batch_sample_modal_repeat`: Specifies the validated repeat length for
  samples used in repeat correction.

- `metrics_group_id`: Groups samples for shared index peak assignment.

- `metrics_baseline_control`: Identifies samples used as baseline
  controls for index peak assignment.

## Examples

``` r
fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())

# import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
fragments_list <- trace(fsa_list, grouped = TRUE, metadata_data.frame = metadata)
#> Finding ladders
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding fragments
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding alleles
#> Calling repeats
#> Assigning index peaks

metrics <- calculate_instability_metrics(
  fragments_list,
  peak_threshold = 0.05,
  window_around_index_peak = c(-40, 40),
  percentile_range = c(0.5, 0.75, 0.9, 0.95),
  repeat_range = c(2, 5, 10, 20)
)
```
