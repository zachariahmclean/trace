# trace 1.0.0
-   Major Improvements
    -   Consolidated Workflow with `trace()` Function

        -   Simplified three-step workflow:
            1.  Read in data (e.g., `read_fsa()`)
            2.  Process with main `trace()` function
            3.  Analyze with extraction functions (e.g., `calculate_instability_metrics()`)
        -   Deprecated functions (now internal-only):
            -   `add_metadata()`, `find_ladders()`, `find_fragments()`, `find_alleles()`, `call_repeats()`, `assign_index_peaks()`
            -   If used for legacy reasons, these now require configuration via YAML file
        -   Configuration options:
            -   Parameters can be passed via `...` or YAML config file
            -   YAML support enables reproducible analysis configurations and easier sharing of processing parameters

    -   Enhanced Ladder Assignment Algorithm

        -   New branching approach selects best overall fit rather than greedy assignment
            -   Parameters:
                -   `ladder_top_n_branching`: Controls number of branches evaluated at each step
                -   `ladder_branching_r_squared_threshold`: Aggressiveness of branch pruning
                -   `ladder_assign_left_to_right`: Direction of assignment (small→large or large→small)
        -   Parameter updates:
            -   `minimum_peak_signal` → `minimum_ladder_signal`
            -   `scan_subset` split into `min_scan` and `max_scan`

    -   Improved Metrics Calculations

        -   Enhanced quality filtering in `calculate_instability_metrics()`:
            -   New parameters: `index_modal_signal_threshold`, `index_signal_sum_threshold`
        -   Fixed calculations:
            -   Corrected weighting in skewness and kurtosis calculations
-   Data Handling Improvements
    -   Streamlined Data Import

        -   Function updates:
            -   `repeat_table_to_repeats()` → `repeat_table_to_fragments()` (now requires specific column names)
            -   `peak_table_to_fragments()` split into:
                -   `genemapper_table_to_fragments()` (for GeneMapper files)
                -   `size_table_to_fragments()` (generic import)

    -   Strict Metadata Requirements

        -   Required columns (must exist even if empty):
            -   `unique_id`, `metrics_group_id`, `metrics_baseline_control` , `batch_run_id`, `batch_sample_id`, `batch_sample_modal_repeat`
-   Other Changes
    -   Parameter Default Updates

        -   Changed many defaults from `NULL` to `NA` to support configuration files

    -   Fragment Detection Improvements

        -   New `peak_scan_ramp` parameter in `trace()` (default relaxed from 6 to 5 to improve peak detection)

    -   Internal Changes

        -   `find_fragments()` now modifies in place (advanced usage only)

# trace 0.5.0

-   `expansion_ratio` metric updated so that it starts at 1 rather than 0

    -   It now includes the relative signal of the index peak in the metric

-   `find_fragments` filters `minimum_peak_signal` on the raw signal rather than the smoothed signal

-   `extract_alleles` now will always return allele_2 so that the column names are always consistent

# trace 0.4.0

Along with a new name, this is a major update to the [previous version of the package](https://github.com/zachariahmclean/instability). The updates below are changes relative to the main branch of the 'instability' github branch.

## Breaking changes

-   Most main functions in the pipeline modify in place rather than return a new object for efficient use of memory

    -   This applies to the following functions: `find_ladders()`, `add_metadata()`, `find_alleles()`, `call_repeats()`, and `assign_index_peaks()`.

-   Repeat correction changes

    -   The argument `repeat_length_correction` has been removed `from call_repeats()` and been replaced with `correction`
    -   Now you chose between "repeat" or "batch" correction (see `call_repeats()` documentation for more info)

-   Metadata categories renamed

    -   Renamed for clarity and to make it clearer which functionality it's associated with
    -   plate_id -\> batch_id, group_id -\> metrics_group_id
    -   size_standard_repeat_length -\> batch_sample_modal_repeat
    -   "size_standard" removed due to unnecessary redundancy of category

-   Renamed instability metrics for clarity:

    -   `average_repeat_gain` to `average_repeat_change`
    -   `modal_repeat_delta` to `modal_repeat_change`

-   Renamed `number_of_peaks_to_return` to `number_of_alleles` in `find_alleles()`

-   Allele and index peak information can no longer be directly interacted with in the class

    -   Setters (`fragment_repeats$set_allele_peak()` & `fragment_repeats$set_index_peak()`) and getters (`fragment_repeats$get_allele_peak()` & `fragment_repeats$get_index_peak()`) have been introduced so users don't accidentally break things if they chose to directly interact with objects.

-   Removed backwards compatibility of index assignment in `calculate_instability_metrics()`

    -   `grouped` and `index_override_dataframe` parameters removed from `calculate_instability_metrics()`

-   Renamed `spike_location` parameter in `find_ladders()` to `ladder_start_scan`

-   Removed `fix_ladders_auto()` function

-   The function `generate_instability_template()` was renamed to `generate_trace_template()`

-   `fragment_trace` class is now instantiated after `read_fsa()` function rather than `find_ladders()`

-   Throughout the package renamed any mention of `height` to `signal` for consistency

## New features

-   Batch correction

    -   Correct systematic batch effects that occur across fragment analysis runs by using common samples across runs without having to know the repeat length

    -   New associated function `plot_batch_correction_samples()` to help visualize this batch correction

    -   See documentation for `call_repeats()` for more info

-   New instability metric `instability_index_change`

-   New function `extract_ladder_summary()` to get a summary of the ladder fits for all samples

-   New function `plot_data_channels()` to help identify the channels containing the ladder and raw data

## Bug fixes

-   If some samples do not having a metrics_group_id or metrics_baseline_control the calculate_instability_metrics() function will still run but those samples will have missing values
