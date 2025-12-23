# Find Alleles

This function identifies main allele within each fragment object.

## Usage

``` r
find_alleles(fragments_list, config, ...)
```

## Arguments

- fragments_list:

  A list of fragment objects containing peak data.

- config:

  A trace_config object generated using
  [`load_config()`](https://zachariahmclean.github.io/trace/reference/load_config.md).

- ...:

  additional parameters from any of the functions in the pipeline
  detailed below may be passed to this function. This overwrites values
  in the `config`. These parameters include:

  - `number_of_alleles` Number of alleles to be returned for each
    fragment. Must either be 1 or 2. Being able to identify two alleles
    is for cases when you are analyze different human samples with a
    normal and expanded alleles and you can't do the preferred option of
    simply ignoring the normal allele in
    [`find_fragments()`](https://zachariahmclean.github.io/trace/reference/find_fragments.md)
    (eg setting the min_bp_size above the normal allele bp size).
    Default: `1`.

  - `peak_region_size_gap_threshold` Gap threshold for identifying peak
    regions. The peak_region_size_gap_threshold is a parameter used to
    determine the maximum allowed gap between peak sizes within a peak
    region. Adjusting this parameter affects the size range of peaks
    that can be grouped together in a region. A smaller value makes it
    more stringent, while a larger value groups peaks with greater size
    differences, leading to broader peak regions that may encompass
    wider size ranges. Default: `6`.

  - `peak_region_signal_threshold_multiplier` Multiplier for the peak
    signal threshold. The peak_region_signal_threshold_multiplier
    parameter allows adjusting the threshold for identifying peak
    regions based on peak signals. Increasing this multiplier value will
    result in higher thresholds, making it more stringent to consider
    peaks as part of a peak region. Conversely, reducing the multiplier
    value will make the criteria less strict, potentially leading to
    more peaks being grouped into peak regions. It's important to note
    that this parameter's optimal value depends on the characteristics
    of the data and the specific analysis goals. Choosing an appropriate
    value for this parameter can help in accurately identifying
    meaningful peak regions in the data. Default: `1`.

## Value

This function modifies list of fragments objects in place with alleles
added.

## Details

This function finds the main alleles for each fragment in the list by
identifying clusters of peaks ("peak regions") with the highest signal
intensities. This is based on the idea that PCR amplicons of repeats
have clusters of peaks (from somatic mosaicism and PCR artifacts) that
help differentiate the main allele of interest from capillary
electrophoresis noise/contamination.

If number_of_alleles = 1, the tallest of peaks will be selected as the
allele. This means that if your sample has multiple alleles, you have
two options i) make sure that your data is subsetted to only include the
allele of interest (using `min_bp_size` in
[`find_fragments()`](https://zachariahmclean.github.io/trace/reference/find_fragments.md)
to make sure that the smaller allele is excluded), or ii) setting
number_of_alleles = 2, which will pick the two tallest peaks in their
respective peak regions and set the main allele as the larger repeat
size, and allele_2 as the shorter repeat size. We recommend the
subsetting approach since that is far simpler and less likely to fail,
and the second option only if you're doing an experiment analysis a
large number of human samples where both the normal and expanded allele
repeat lengths vary, which makes it very difficult to find a common bp
size that excludes the normal allele.

The parameters `peak_region_signal_threshold_multiplier` and
`peak_region_size_gap_threshold` will only need to be adjusted in rare
cases if peaks are not being found for some reason. They influence the
criteria for identifying peak regions.
peak_region_signal_threshold_multiplier is multiplied to the mean height
of all the peaks to create a hight threshold for inclusion into the peak
region, so most of the time it's already a very low value and probably
only needs to be changed if you have very few peaks.
peak_region_size_gap_threshold is the distance between the peaks, either
bp size, or repeats if repeats have already been called.

## Examples

``` r
fsa_list <- lapply(cell_line_fsa_list[1], function(x) x$clone())
config <- load_config()

trace:::find_ladders(fsa_list, config, show_progress_bar = FALSE)

trace:::find_fragments(fsa_list,
  config,
  min_bp_size = 300,
  show_progress_bar = FALSE
)


trace:::find_alleles(
  fsa_list,
  config,
  peak_region_signal_threshold_multiplier = 1
)
```
