# Load Trace Config File

Load configuration file that can be passed to individual trace modules.

## Usage

``` r
load_config(config_file = NULL)
```

## Arguments

- config_file:

  The file path to a YAML file containing a full list of parameters. If
  NULL, it will load a default YAML.

## Value

A trace_config object.

## Details

Use this function to load in a config file for passing to individual
trace modules. This is only required if creating a custom pipeline
rather than the main function. This provides a central place to adjust
parameters for the pipeline. Either edit the YAML directly, or modify
the object (which is basically just a list).

Use the following command to make a copy of the YAML file:
`file.copy(system.file("extdata/trace_config.yaml", package = "trace"), ".")`.

## Examples

``` r
config <- load_config()
```
