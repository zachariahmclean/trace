# Read fsa file

Read fsa file into memory and create fragments object

## Usage

``` r
read_fsa(files)
```

## Arguments

- files:

  a chr vector of fsa file names. For example, return all the fsa files
  in a directory with 'list.files("example_directory/", full.names =
  TRUE, pattern = ".fsa")'.

## Value

A list of fragments objects

## Details

read_fsa is just a wrapper around
[`seqinr::read.abif()`](https://rdrr.io/pkg/seqinr/man/read.abif.html)
that reads the fsa file into memory and stores it inside a fragments
object. That enables you to use the next function
[`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md).

## See also

[`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md),
[`plot_data_channels()`](https://zachariahmclean.github.io/trace/reference/plot_data_channels.md)

## Examples

``` r
fsa_file <- read_fsa(system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr"))
plot_data_channels(fsa_file)

```
