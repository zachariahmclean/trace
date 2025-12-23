# Fix ladders interactively

An app for fixing ladders

## Usage

``` r
fix_ladders_interactive(fragment_trace_list)
```

## Arguments

- fragment_trace_list:

  A list of fragments objects containing fragment data

## Value

interactive shiny app for fixing ladders

## Details

This function helps you fix ladders that are incorrectly assigned. Run
`fix_ladders_interactive()` and provide output from `find_ladders`. In
the app, for each sample, click on line for the incorrect ladder size
and drag it to the correct peak.

Once you are satisfied with the ladders for all the broken samples,
click the download button to generate a file that has the ladder
correction data. Read this file back into R using readRDS, then use
[`fix_ladders_manual()`](https://zachariahmclean.github.io/trace/reference/fix_ladders_manual.md)
and supply the ladder correction data as `ladder_df_list`. This allows
the manually corrected data to be saved and used within a script so that
the correct does not need to be done every time. An example of what you
would need to do:

ladder_df_list \<- readRDS('path/to/exported/data.rds')
test_ladders_fixed \<- fix_ladders_manual(test_ladders_broken,
ladder_df_list)

## See also

[`fix_ladders_manual()`](https://zachariahmclean.github.io/trace/reference/fix_ladders_manual.md),
[`find_ladders()`](https://zachariahmclean.github.io/trace/reference/find_ladders.md)

## Examples

``` r
# to create an example, lets brake one of the ladders
brake_ladder_list <-  list(
   "20230413_A08.fsa" = data.frame(
     size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
     scan = c(1544, 1621, 1850, 1912, 2143, 2201, 2261, 2506, 2805, 3135, 3380, 3442, 3760, 
              4050, 4284, 4332)
   )
 )

fsa_list <- lapply(cell_line_fsa_list, function(x) x$clone())
# import data with read_fsa() to generate an equivalent list to cell_line_fsa_list
test_fragments <- trace(fsa_list, ladder_df_list = brake_ladder_list)
#> Finding ladders
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Fixing ladder for 20230413_A08.fsa
#> < fix_ladders_manual warning >
#> sample 20230413_A08.fsa has badly fitting ladder for bp sizes: 35-75, 50-100, 75-139
#> Finding fragments
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  11%  |                                                                              |===========                                                           |  16%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  26%  |                                                                              |======================                                                |  32%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  42%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  63%  |                                                                              |================================================                      |  68%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  79%  |                                                                              |===========================================================           |  84%  |                                                                              |===============================================================       |  89%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Finding alleles
#> Calling repeats
#> Assigning index peaks

if (interactive()) {
  fix_ladders_interactive(test_fragments)
}

# once you have corrected your ladders in the app, 
# export the data for incorporation into the script.
# You can then re-import the data and fix ladders as described in the help details.
```
