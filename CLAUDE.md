# CLAUDE.md

Guidance for working on **trace** — an R package for Tandem Repeat Analysis by Capillary Electrophoresis (STR instability metrics from fragment-analysis `.fsa` files). See the paper: Jiang et al. 2026, *J Huntington's Dis*.

> Machine-specific setup (R executable path, which shell to run R through, Pandoc availability, absolute paths) lives in **`CLAUDE.local.md`** (untracked). Read it for the exact command invocations on this machine.

## Dev workflow

- **Tests** (use `load_all`, don't install):
  ```r
  devtools::load_all(".", quiet = TRUE)
  testthat::test_file("tests/testthat/test-<name>.R")
  ```
  Full sweep: loop `test_file()` over `list.files("tests/testthat", "^test-.*\\.R$")`, skipping `shinytest2`/`plotting` (slow/interactive).
- **Docs:** run `devtools::document(".")` after any roxygen change — both `@export` (updates NAMESPACE) and `@field`/`@param` (updates `man/*.Rd`) require it.
- **Check:** `devtools::check(".")` (see `CLAUDE.local.md` for the flags this environment needs).

## Pipeline & architecture

User-facing flow: `read_fsa()` (or `*_to_fragments()`) → `trace()` → `calculate_instability_metrics()` / `extract_*()`.

`trace()` ([R/main.R](R/main.R)) orchestrates internal steps (each in its own file, callable as `trace:::fn()`):
`add_metadata` → `set_batch_run_id_from_fsa` → `find_ladders` → `find_fragments` → `find_alleles` → `call_repeats` → `assign_index_peaks`.

- **`fragments` is an R6 class** ([R/classes.R](R/classes.R)), mutated **by reference** — internal steps and `add_metadata()` change objects in place and often return only a `trace_output` status object (warnings in `$warning_message` (list); status in `$status`). `trace()` clones its inputs first.
- **Config** is a strict S3 list ([R/config.R](R/config.R), `inst/extdata/trace_config.yaml`). `validate_inputs()` errors on any unknown or missing key, so **a new parameter must be added to BOTH the YAML and the `expected` list**. Params reach functions via `...` or a YAML `config_file`.
- FSA files are read with `seqinr::read.abif()`; the whole ABIF object is stored on `$fsa`, parsed header fields on `$fsa_metadata` ([R/fsa_metadata.R](R/fsa_metadata.R)).

## Conventions

- Base R style throughout (no tidyverse in package code; dplyr/ggplot2 are Suggests, used only in examples/tests). There is no `%||%` operator — write explicit NA guards.
- Internal helpers are unexported, dot-prefixed (`.fsa_chr`) or marked `@keywords internal`.
- Exported functions get roxygen with `@export` and `@examples` (examples typically `lapply(cell_line_fsa_list, \(x) x$clone())` then `trace(...)`).
- Tests: testthat 3e. Fixtures use bundled data (`cell_line_fsa_list`, `metadata`, `example_data`) and the always-available `system.file("abif/...", package = "seqinr")` (a different-instrument `.fsa`, good for generality). No `.fsa` files ship in the package (`data-raw/` is `.Rbuildignore`d, not run in checks).

## Gotchas

- Changing the bundled `metadata` / `cell_line_fsa_list` data objects breaks tests that hard-code values (e.g. `batch_run_id` `"20230414"`/`"20220630"`). Prefer leaving bundled data alone.
- `read_fsa()` builds `unique_id` = `SpNm` + `_` + `RunN` by default; `trace()` derives `batch_run_id` from the fsa by default (`auto_batch_run_id`). Both are opt-out. Examples relying on hand-set metadata `batch_run_id` must pass `auto_batch_run_id = FALSE`.
- Dev docs `REFACTOR_SPEC.md` and `CLAUDE.md`/`CLAUDE.local.md`, plus `.claude/`, are `.Rbuildignore`d so they don't ship in the package.
