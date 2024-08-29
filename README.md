
# instability package

This package provides a pipeline for short tandem repeat instability
analysis from fragment analysis data. The inputs are fsa files or peak
tables (eg Genemapper 5 software output), and a user supplied metadata
data-frame. The functions identify ladders, calls peaks, and calculate
repeat instability metrics (ie expansion index or average repeat gain).

To report bugs or feature requests, please visit the Github issue
tracker [here](https://github.com/zachariahmclean/instability/issues).
For assistance or any other inquires, contact [Zach
McLean](mailto:zmclean@mgh.harvard.edu?subject=%5BGitHub%5D%20instability).

If you use this package, please cite
[this](https://www.nature.com/articles/s41467-024-47485-0) paper for
now.

# How to use the package

For an easy way to get started with your own data or to run an example,
use `trace::generate_instability_template()` to generate a
document with the pipeline pre-populated.

In this package, each sample is represented by an R6 ‘fragments’ object,
which are organised in lists. As a user, there are accessor functions
that iterate over these lists, so you shouldn’t need to interact with
the fragments object. However, if you do, the attributes of the objects
can be accessed with “\$”.

There are several important factors to a successful repeat instability
experiment and things to consider when using this package:

- (required) Each sample has a unique id, usually the file name

- (optional) Baseline control for your experiment. For example,
  specifying a sample where the modal allele is the inherited repeat
  length (eg a mouse tail sample) or a sample at the start of a
  time-course experiment. This is indicated with a `TRUE` in the
  `metrics_baseline_control` column of the metadata. For mice, if just a
  few samples have the inherited repeat height shorter than the expanded
  population, you could not worry about this and instead use the
  `index_override_dataframe` in `calculate_instability_metrics()`

- (optional) Using a common sample(s) of known repeat length to correct
  repeat size across fragment analysis runs. There are slight
  fluctuations of repeat length size across runs, so if samples are to
  be analyzed for different runs, you must correct the repeat length so
  they are comparable. This is usually achieved by running running
  positive control samples with a known validated repeat size of the
  modal peak in each fragment analysis run. These samples are then
  indicated in the metadata with `batch_run_id` (to group samples by
  fragment analysis run), `size_standard` (to indicate which sames are
  size standards), `size_standard_sample_id` (to enable checking that
  the modal peak is consistent across runs and plotting functions to
  check correction results), and `size_standard_repeat_length`
  (indicating the repeat length of the modal peak for each size standard
  sample).

- If starting from fsa files, the GeneScan™ 1200 LIZ™ dye Size Standard
  ladder assignment may not work very well. The ladder identification
  algorithm is optimized for GeneScan™ 500 LIZ™ or GeneScan™ 600 LIZ™.
  The 1200 LIZ™ ladder has a challenging pattern of ladder peaks to
  automatically assign. However, these ladders can be fixed by playing
  with the various parameters or manually with the built-in
  fix_ladders_interactive() app.

# Installation

You can install from
[GitHub](https://github.com/zachariahmclean/instability) with:

``` r
# install.packages("devtools")
devtools::install_github("zachariahmclean/instability")
```

Then load the package:

``` r
library(instability)
```

# Import data

First, we read in the raw data. In this case we will used example data
within this package, but usually this would be fsa files that are read
in using `read_fsa()`.

``` r
fsa_raw <- trace::cell_line_fsa_list
```

# Find ladders

The raw data are coerced into a list of ‘fragments’ class objects, which
is a fundamental data structure used in this pipeline. The ‘fragments’
class objects are R6 classes, so the individual elements can be accesses
with “\$”.

First we find the ladders and call bp size in the fsa file. The bp is
assigned using the local Southern method. Basically, for each data
point, linear models are made for the lower and upper 3 size standard
and the predicted sizes are averaged.

``` r
ladder_list <- find_ladders(cell_line_fsa_list,
  show_progress_bar = FALSE
)
```

visually inspect each ladder to make sure that the ladders were
correctly assigned

``` r
plot_ladders(ladder_list[1])
```

<img src="man/figures/README-plot_ladders-1.png" width="100%" />

If the ladders are are not assigned correctly, you can either try
fix_ladders_auto() (optimal for when just a single ladder peak is
wrong), or manually using the built-in fix_ladders_interactive() app.

![](man/figures/ladder_fixing.gif)

# Find fragments

The fragment peaks are identified in the raw continuous trace data.

``` r
peak_list <- find_fragments(ladder_list,
  min_bp_size = 300
)
```

Visually inspect the traces and called peaks to make sure they were
correctly assigned.

``` r
plot_traces(peak_list[1],
  xlim = c(400, 550),
  ylim = c(0, 1200)
)
```

<img src="man/figures/README-plot_traces-1.png" width="100%" />

Alternatively, if not starting from fsa files, this is where you would
use exported data from Genemapper if you would rather use the Genemapper
bp sizing and peak identification algorithms.

``` r
peak_list_genemapper <- peak_table_to_fragments(trace::example_data,
  data_format = "genemapper5",
  dye_channel = "B",
  min_size_bp = 300
)
```

# Add metadata

Metadata can be incorporated to allow additional functionality in the
`call_repeats()` (correcting repeat length across fragment analysis
runs) and `calculate_instability_metrics()` (calculating expansion index
or average repeat gain) functions. Some of the metadata fields are
optional but highly recommended. Prepare a file (eg spreadsheet saved as
.csv) with the following columns. If you use the specified column names,
it will be automatically parsed by `add_metadata()`, otherwise you will
need to match up which column name belongs to which metadata category
(as done below in `add_metadata()`):

| Metadata table column       | Functionality metadata is associated with                   | Description                                                                                                                                                                                                                                                                                                     |
|-----------------------------|-------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| unique_id                   | Required for adding metadata using `add_metdata()`          | The unique identifier for the fsa file. Usually the sample file name. This must be unique, including across runs.                                                                                                                                                                                               |
| group_id                    | `assign_index_peaks()`, allows setting `grouped = TRUE`     | This groups the samples for instability metric calculations. Provide a group id value for each sample. For example, in a mouse experiment and using the expansion index, you need to group the samples since they have the same metrics baseline control (eg inherited repeat length), so provide the mouse id. |
| metrics_baseline_control    | `assign_index_peaks()`, allows setting `grouped = TRUE`     | This is related to group_id. Indicate with ‘TRUE’ to specify which sample is the baseline control (eg mouse tail for inherited repeat length, or day-zero sample in cell line experiments)                                                                                                                      |
| batch_run_id                    | `call_repeats()`, allows setting `repeat_length_correction` | This groups the samples for correcting the repeat length. Provide a value for each fragment analysis run (eg date).                                                                                                                                                                                             |
| size_standard               | `call_repeats()`, allows setting `repeat_length_correction` | This is related to batch_run_id. Indicate with ‘TRUE’ to specify which sample is the size standard of the repeat length.                                                                                                                                                                                            |
| size_standard_sample_id     | `call_repeats()`, allows setting `repeat_length_correction` | Give an id for the size standard. This allows comparison of the size standard across fragment analysis runs and plotting to visualize corrections.                                                                                                                                                              |
| size_standard_repeat_length | `call_repeats()`, allows setting `repeat_length_correction` | This is related to size_standard. If the sample is a size standard, provide a numeric value of the modal repeat length.                                                                                                                                                                                         |

``` r

metadata_added_list <- add_metadata(
  fragments_list = peak_list,
  metadata_data.frame = trace::metadata,
  unique_id = "unique_id",
  group_id = "group_id",
  metrics_baseline_control = "metrics_baseline_control",
  batch_run_id = "batch_run_id",
  size_standard = "size_standard",
  size_standard_repeat_length = "size_standard_repeat_length"
)
```

# Identify modal peaks and call repeats

Next we identify the modal peaks with `find_alleles()` and convert the
base pair fragments to repeats with `call_repeats(`).

``` r
alleles_list <- find_alleles(
  fragments_list = metadata_added_list
)


repeats_list <- call_repeats(
  fragments_list = alleles_list,
  repeat_length_correction = "from_metadata"
)
```

We can view the distribution of repeat sizes and the identified modal
peak with a plotting function.

``` r
plot_traces(repeats_list[1], xlim = c(110, 150))
```

<img src="man/figures/README-plot_fragments-1.png" width="100%" />

We can also view the data used to generate the model for calling the
repeat size when we indicate size standard samples in the metadata and
have `repeat_length_correction = "from_metadata"` in `call_repeats()`.

``` r
plot_size_standard_model(repeats_list)
```

<img src="man/figures/README-plot_size_standard_model-1.png" width="100%" />

In this case the dots are basically overlapping and in the middle of the
linear model, indicating that we have correctly identified the known
tallest peak used for the repeat size standards. If the wrong peak was
selected for one of the samples, the dots would be shifted across 3 bp
and no longer overlapping.

# Assign index peaks

A key part of several instability metrics is the index peak. This is the
repeat length used as the reference for relative instability metrics
calculations, like expansion index or average repeat gain. In the
metadata, samples are grouped by a `group_id` and a subset of the
samples are set as `metrics_baseline_control`, meaning they are the
samples taken at day 0 in this experiment. This allows us to set
`grouped = TRUE` and set the index peak for the expansion index and
other metrics. For mice, if just a few samples have the inherited repeat
height shorter than the expanded population, you could not worry about
this and instead use the `index_override_dataframe` in
`calculate_instability_metrics()`.

``` r

index_list <- assign_index_peaks(
  repeats_list,
  grouped = TRUE
)
```

We can validate that the index peaks were assigned correctly with a
dotted vertical line added to the trace. This is perhaps more useful in
the context of mice where you can visually see when the inherited repeat
length should be in the bimodal distribution.

``` r
plot_traces(index_list[1], xlim = c(110, 150))
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

# Calculate instability metrics

All of the information we need to calculate the repeat instability
metrics has now been identified. We can finally use
`calculate_instability_metrics` to generate a dataframe of per-sample
metrics.

``` r
metrics_grouped_df <- calculate_instability_metrics(
  fragments_list = index_list,
  peak_threshold = 0.05
)
```

These metrics can then be used to quantify repeat instability. For
example, this reproduces Figure 7e of [our
manuscript](https://www.nature.com/articles/s41467-024-47485-0).

First, prepare the data for plotting by removing poor quality samples
and finding the average repeat gain relative to the DMSO group for each
cell line

``` r
library(dplyr)


plot_data <- metrics_grouped_df |>
  dplyr::left_join(metadata, by = dplyr::join_by(unique_id)) |>
  dplyr::filter(
    day > 0,
    modal_peak_height > 500
  ) |>
  dplyr::group_by(group_id) |>
  dplyr::mutate(
    rel_gain = average_repeat_gain / median(average_repeat_gain[which(treatment == 0)]),
    genotype = forcats::fct_rev(genotype)
  )
```

Then we can plot the instability metrics

``` r
library(ggplot2)

ggplot(
  plot_data,
  aes(genotype, rel_gain, colour = genotype)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(vars(as.factor(treatment)), nrow = 1) +
  labs(
    y = "Average repeat gain\n(relative to DMSO)",
    x = "PMS1 pseudoexon status"
  ) +
  theme(legend.position = "none")
```

<img src="man/figures/README-ggplot-1.png" width="100%" />
