# fragments object

An R6 Class representing a fragments object.

## Details

This is the parent class of both fragments and fragments object. The
idea is that shared fields and methods are both inherited from this
object, but it is not itself directly used.

## Public fields

- `unique_id`:

  unique id of the sample usually the file name

- `input_method`:

  pathway used to import data (either 'fsa', 'size', or 'repeats')

- `metrics_group_id`:

  sample grouping for metrics calculations. Associated with
  [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md).

- `metrics_baseline_control`:

  logical to indicate if sample is the baseline control. Associated with
  [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md).

- `batch_run_id`:

  fragment analysis run. Associated with
  [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md).

- `batch_sample_id`:

  An id for the sample used as size standard for repeat calculation.
  Associated with
  [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md).

- `batch_sample_modal_repeat`:

  Validated repeat length for the modal repeat repeat in that sample.
  Associated with
  [`add_metadata()`](https://zachariahmclean.github.io/trace/reference/add_metadata.md).

- `fsa`:

  The whole fsa file, output from seqinr::read.abif()

- `raw_ladder`:

  The raw data from the ladder channel

- `raw_data`:

  The raw data from the sample channel

- `scan`:

  The scan number

- `off_scale_scans`:

  vector indicating which scales were too big and off scale. Note can be
  in any channel

- `ladder_df`:

  A dataframe of the identified ladder from
  [`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md).
  Scan is the scan number of peak and size is the associated bp size.

- `ladder_total_combinations_tested`:

  A numeric value indicating how many total combinations were tested
  during ladder fit

- `trace_bp_df`:

  A dataframe of bp size for every scan from
  [`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md).

- `peak_table_df`:

  A dataframe containing the fragment peak level information.

- `repeat_table_df`:

  A dataframe containing the fragment peak level information with the
  repeat size added. May or may not be the same as peak_table_df
  depending on what options are chosen in `call_repeats`.

## Methods

### Public methods

- [`fragments$new()`](#method-fragments-new)

- [`fragments$print()`](#method-fragments-print)

- [`fragments$plot_trace()`](#method-fragments-plot_trace)

- [`fragments$plot_ladder()`](#method-fragments-plot_ladder)

- [`fragments$plot_data_channels()`](#method-fragments-plot_data_channels)

- [`fragments$get_allele_peak()`](#method-fragments-get_allele_peak)

- [`fragments$set_allele_peak()`](#method-fragments-set_allele_peak)

- [`fragments$get_index_peak()`](#method-fragments-get_index_peak)

- [`fragments$set_index_peak()`](#method-fragments-set_index_peak)

- [`fragments$plot_fragments()`](#method-fragments-plot_fragments)

- [`fragments$clone()`](#method-fragments-clone)

------------------------------------------------------------------------

### Method `new()`

initialization function that is not used since the child classes are the
main object of this package.

#### Usage

    fragments$new(unique_id, input_method, object)

#### Arguments

- `unique_id`:

  unique_id

- `input_method`:

  pathway used to import data (either 'fsa', 'size', or 'repeats')

- `object`:

  object to be inserted into either 'fsa', 'peak_table_df', or
  'repeat_table_df'

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

A function to print informative information to the console

#### Usage

    fragments$print()

------------------------------------------------------------------------

### Method `plot_trace()`

plot the trace data

#### Usage

    fragments$plot_trace(
      show_peaks = TRUE,
      x_axis = NULL,
      ylim = NULL,
      xlim = NULL,
      signal_color_threshold = 0.05,
      plot_title = NULL
    )

#### Arguments

- `show_peaks`:

  A logical to say if the called peaks should be plotted on top of the
  trace. Only valid for fragments objects.

- `x_axis`:

  Either "size" or "repeats" to indicate what should be plotted on the
  x-axis.

- `ylim`:

  numeric vector length two specifying the y axis limits

- `xlim`:

  numeric vector length two specifying the x axis limits

- `signal_color_threshold`:

  A threshold value to colour the peaks relative to the tallest peak.

- `plot_title`:

  A character string for setting the plot title. Defaults to the unique
  id of the object

#### Returns

A base R plot

------------------------------------------------------------------------

### Method `plot_ladder()`

plot the ladder data

#### Usage

    fragments$plot_ladder(xlim = NULL, ylim = NULL, plot_title = NULL)

#### Arguments

- `xlim`:

  numeric vector length two specifying the x axis limits

- `ylim`:

  numeric vector length two specifying the y axis limits

- `plot_title`:

  A character string for setting the plot title. Defaults to the unique
  id of the object

#### Returns

A base R plot

------------------------------------------------------------------------

### Method [`plot_data_channels()`](https://zachariahmclean.github.io/trace/reference/plot_data_channels.md)

plot the raw data channels in the fsa file. It identifies every channel
that has "DATA" in its name.

#### Usage

    fragments$plot_data_channels()

#### Returns

A base R plot

------------------------------------------------------------------------

### Method `get_allele_peak()`

This returns a list with the allele information for this object.

#### Usage

    fragments$get_allele_peak()

------------------------------------------------------------------------

### Method `set_allele_peak()`

This sets a single allele size/repeat. It searches through the
appropriate peak table and finds the closest peak to the value that's
provided.

#### Usage

    fragments$set_allele_peak(allele, unit, value)

#### Arguments

- `allele`:

  Either `1` or `2`, indicating which allele information should be set.
  Allele 1 is the only one used for repeat instability metrics
  calculations.

- `unit`:

  Either "size" or "repeats" to indicate if the value you're providing
  is bp size or repeat length.

- `value`:

  Numeric vector (length one) of the size/repeat length to set.

------------------------------------------------------------------------

### Method `get_index_peak()`

This returns a list with the index peak information for this object.

#### Usage

    fragments$get_index_peak()

------------------------------------------------------------------------

### Method `set_index_peak()`

This sets the index repeat length. It searches through the repeat table
and finds the closest peak to the value that's provided.

#### Usage

    fragments$set_index_peak(value)

#### Arguments

- `value`:

  Numeric vector (length one) of the repeat length to set as index peak.

------------------------------------------------------------------------

### Method [`plot_fragments()`](https://zachariahmclean.github.io/trace/reference/plot_fragments.md)

This plots the peak/repeat table as a histogram

#### Usage

    fragments$plot_fragments(ylim = NULL, xlim = NULL, plot_title = NULL)

#### Arguments

- `ylim`:

  numeric vector length two specifying the y axis limits

- `xlim`:

  numeric vector length two specifying the x axis limits

- `plot_title`:

  A character string for setting the plot title. Defaults to the unique
  id of the object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    fragments$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
