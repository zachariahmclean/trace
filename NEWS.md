# trace 0.2.0

# trace 0.1.0

Along with a new name, this is a major update to the [previous version of the package](https://github.com/zachariahmclean/instability). The updates below are changes relative to the main branch of the 'instability' github branch.

The nature of these changes are on focusing this tool for instability metrics calculation rather and accurate genotyping of repeat length.

## Breaking changes

-   Most main functions in the pipeline modify in place rather than return a new object for efficient use of memory.

    -   This applies to the following functions: `find_ladders()`, `add_metadata()`, `find_alleles()`, `call_repeats()`, and `assign_index_peaks()`.

-   Removed repeat calling using size standard of known repeat length

    -   The argument `repeat_length_correction` has been removed `from call_repeats()`.
    -   This functionality has been transferred to a newer and much simplier approach of batch correction (see documentation for `call_repeats()`.
    -   The associated plotting function `plot_size_standard_model()` was also removed.

-   Metadata categories renamed

    -   Renamed for clarity and to make it clearer which functionality it's associated with
    -   plate_id -\> batch_id, group_id -\> metrics_group_id, and as mentioned above size_standard & size_standard_repeat_length removed

-   Now only possible to find one allele

    -   Removed `number_of_peaks_to_return` in `find_alleles()`.

-   Allele and index peak information can no longer be directly interacted with in the class

    -   Setters (`fragment_repeats$set_allele_peak()` & f`ragment_repeats$set_index_peak()`) and getters (`fragment_repeats$get_allele_peak()` & `fragment_repeats$get_index_peak()`) have been introduced so users don't accidentally break things if they chose to directly interact with objects.

-   Removed backwards compatibility of index assignment in `calculate_instability_metrics()`

    -   `grouped` and `index_override_dataframe` parameters removed from `calculate_instability_metrics()`.

-   Renamed `spike_location` parameter in `find_ladders()` to `ladder_start_scan`

-   Removed `fix_ladders_auto()` function

-   The function `generate_instability_template()` was renamed to `generate_trace_template()`

-   `fragment_trace` class is now instantiated after `read_fsa()` function rather than `find_ladders()`

## New features

-   Batch correction

    -   Correct systematic batch effects that occur across fragment analysis runs by using common samples across runs without having to know the repeat length.

    -   New associated function `plot_batch_correction_samples()` to help visualize this batch correction.

    -   See documentation for `call_repeats()` for more info.

-   New function `extract_ladder_summary()` to get a summary of the ladder fits for all samples.

-   New function `plot_data_channels()` to help identify the channels containing the ladder and raw data.

## Bug fixes

-   If some samples do not having a metrics_group_id or metrics_baseline_control the calculate_instability_metrics() function will still run but those samples will have missing values.
