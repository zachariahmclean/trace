# Changelog

## trace 1.0.0

CRAN release: 2025-12-22

- Major Changes
  - Consolidated Workflow with
    [`trace()`](https://zachariahmclean.github.io/trace/reference/trace.md)
    Function

    - Simplified three-step workflow:
      1.  Read in data (e.g.,
          [`read_fsa()`](https://zachariahmclean.github.io/trace/reference/read_fsa.md))
      2.  Process with main
          [`trace()`](https://zachariahmclean.github.io/trace/reference/trace.md)
          function
      3.  Analyze with extraction functions (e.g.,
          [`calculate_instability_metrics()`](https://zachariahmclean.github.io/trace/reference/calculate_instability_metrics.md))
    - Deprecated functions (now internal-only):
      - [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md),
        [`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md),
        [`find_fragments()`](https://zachariahmclean.github.io/trace/reference/find_fragments.md),
        [`find_alleles()`](https://zachariahmclean.github.io/trace/reference/find_alleles.md),
        [`call_repeats()`](https://zachariahmclean.github.io/trace/reference/call_repeats.md),
        [`assign_index_peaks()`](https://zachariahmclean.github.io/trace/reference/assign_index_peaks.md)
      - If used for legacy reasons access via the internal package with
        ‘:::’ (eg `trace:::find_ladders()`), these now require
        configuration via YAML file (see each functions help for more
        detail)
    - Configuration options:
      - Parameters can be passed via `...` or YAML config file
      - YAML support enables reproducible analysis configurations and
        easier sharing of processing parameters

  - Enhanced Ladder Assignment Algorithm

    - New branching approach selects best overall fit rather than greedy
      assignment
      - Parameters:
        - `ladder_top_n_branching`: Controls number of branches
          evaluated at each step
        - `ladder_branching_r_squared_threshold`: Aggressiveness of
          branch pruning
        - `ladder_assign_left_to_right`: Direction of assignment
          (small→large or large→small)
    - Parameter name updates:
      - For ladder processing, `minimum_peak_signal` to
        `minimum_ladder_signal`
      - For ladder processing, `scan_subset` split into `min_scan` and
        `max_scan`

  - Improved Metrics Calculations

    - Enhanced quality filtering in
      [`calculate_instability_metrics()`](https://zachariahmclean.github.io/trace/reference/calculate_instability_metrics.md):
      - New parameters: `index_modal_signal_threshold`,
        `index_signal_sum_threshold`
    - Fixed calculations:
      - Corrected weighting in skewness and kurtosis calculations
      - Fixed issue with average_repeat_change being non-zero for some
        of the index samples
- Data Handling Changes for Package Simplicity
  - Streamlined Data Import

    - Function updates:
      - `repeat_table_to_repeats()` to
        [`repeat_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/repeat_table_to_fragments.md)
        (now requires specific column names)
      - `peak_table_to_fragments()` split into:
        - [`genemapper_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/genemapper_table_to_fragments.md)
          (for GeneMapper files)
        - [`size_table_to_fragments()`](https://zachariahmclean.github.io/trace/reference/size_table_to_fragments.md)
          (generic import)

  - Strict Metadata Requirements

    - Require exact columns (if those column names not detected, those
      metadata will be skipped):
      - `unique_id`, `metrics_group_id`, `metrics_baseline_control` ,
        `batch_run_id`, `batch_sample_id`, `batch_sample_modal_repeat`
- Other Changes
  - Parameter Default Updates

    - Changed many defaults from `NULL` to `NA` to support configuration
      files

  - Fragment Detection Improvements

    - New `peak_scan_ramp` parameter in
      [`trace()`](https://zachariahmclean.github.io/trace/reference/trace.md)
      (default relaxed from 6 to 5 to improve peak detection at longer
      sizes)

  - Internal Changes

    - [`find_fragments()`](https://zachariahmclean.github.io/trace/reference/find_fragments.md)
      now modifies in place (advanced usage only)

## trace 0.5.0

CRAN release: 2025-01-16

- `expansion_ratio` metric updated so that it starts at 1 rather than 0

  - It now includes the relative signal of the index peak in the metric

- `find_fragments` filters `minimum_peak_signal` on the raw signal
  rather than the smoothed signal

- `extract_alleles` now will always return allele_2 so that the column
  names are always consistent

## trace 0.4.0

Along with a new name, this is a major update to the [previous version
of the package](https://github.com/zachariahmclean/instability). The
updates below are changes relative to the main branch of the
‘instability’ github branch.

### Breaking changes

- Most main functions in the pipeline modify in place rather than return
  a new object for efficient use of memory

  - This applies to the following functions:
    [`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md),
    [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md),
    [`find_alleles()`](https://zachariahmclean.github.io/trace/reference/find_alleles.md),
    [`call_repeats()`](https://zachariahmclean.github.io/trace/reference/call_repeats.md),
    and
    [`assign_index_peaks()`](https://zachariahmclean.github.io/trace/reference/assign_index_peaks.md).

- Repeat correction changes

  - The argument `repeat_length_correction` has been removed
    `from call_repeats()` and been replaced with `correction`
  - Now you chose between “repeat” or “batch” correction (see
    [`call_repeats()`](https://zachariahmclean.github.io/trace/reference/call_repeats.md)
    documentation for more info)

- Metadata categories renamed

  - Renamed for clarity and to make it clearer which functionality it’s
    associated with
  - plate_id -\> batch_id, group_id -\> metrics_group_id
  - size_standard_repeat_length -\> batch_sample_modal_repeat
  - “size_standard” removed due to unnecessary redundancy of category

- Renamed instability metrics for clarity:

  - `average_repeat_gain` to `average_repeat_change`
  - `modal_repeat_delta` to `modal_repeat_change`

- Renamed `number_of_peaks_to_return` to `number_of_alleles` in
  [`find_alleles()`](https://zachariahmclean.github.io/trace/reference/find_alleles.md)

- Allele and index peak information can no longer be directly interacted
  with in the class

  - Setters (`fragment_repeats$set_allele_peak()` &
    `fragment_repeats$set_index_peak()`) and getters
    (`fragment_repeats$get_allele_peak()` &
    `fragment_repeats$get_index_peak()`) have been introduced so users
    don’t accidentally break things if they chose to directly interact
    with objects.

- Removed backwards compatibility of index assignment in
  [`calculate_instability_metrics()`](https://zachariahmclean.github.io/trace/reference/calculate_instability_metrics.md)

  - `grouped` and `index_override_dataframe` parameters removed from
    [`calculate_instability_metrics()`](https://zachariahmclean.github.io/trace/reference/calculate_instability_metrics.md)

- Renamed `spike_location` parameter in
  [`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md)
  to `ladder_start_scan`

- Removed `fix_ladders_auto()` function

- The function `generate_instability_template()` was renamed to
  `generate_trace_template()`

- `fragment_trace` class is now instantiated after
  [`read_fsa()`](https://zachariahmclean.github.io/trace/reference/read_fsa.md)
  function rather than
  [`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md)

- Throughout the package renamed any mention of `height` to `signal` for
  consistency

### New features

- Batch correction

  - Correct systematic batch effects that occur across fragment analysis
    runs by using common samples across runs without having to know the
    repeat length

  - New associated function
    [`plot_batch_correction_samples()`](https://zachariahmclean.github.io/trace/reference/plot_batch_correction_samples.md)
    to help visualize this batch correction

  - See documentation for
    [`call_repeats()`](https://zachariahmclean.github.io/trace/reference/call_repeats.md)
    for more info

- New instability metric `instability_index_change`

- New function
  [`extract_ladder_summary()`](https://zachariahmclean.github.io/trace/reference/extract_ladder_summary.md)
  to get a summary of the ladder fits for all samples

- New function
  [`plot_data_channels()`](https://zachariahmclean.github.io/trace/reference/plot_data_channels.md)
  to help identify the channels containing the ladder and raw data

### Bug fixes

- If some samples do not having a metrics_group_id or
  metrics_baseline_control the calculate_instability_metrics() function
  will still run but those samples will have missing values
