# TRACE refactor — implementation spec (Phase 1)

Covers three additive, mostly backward-compatible features:
1. FSA metadata parsing (foundation)
2. Auto-derived `batch_run_id` with cross-check (default on, opt-out)
3. Metadata template generator
4. `qc_report()` sample-QC panel

Grounded in: R6 `fragments` class ([R/classes.R](R/classes.R)), `read_fsa()` ([R/constructors.R](R/constructors.R)), `add_metadata()` ([R/add_metadata.R](R/add_metadata.R)), `trace()` ([R/main.R](R/main.R)), `find_ladders()` ([R/find_ladders.R](R/find_ladders.R)), `extract_ladder_summary()` ([R/extract_functions.R](R/extract_functions.R)).

Empirical basis (project 2 Q111, 386 samples/35 runs + curated good/bad): telemetry does **not** explain batch effects; modal prominence ↓ with repeat size (ρ=−0.50); bad samples fail heterogeneously so QC must be a multi-metric panel using `min_rsq` (not `avg_rsq`), `n_peaks`, modal signal, saturation.

---

## 0. ABIF tags (confirmed present in all sampled files; parse defensively — other instruments may omit)

| Tag | Meaning | Used for |
|-----|---------|----------|
| `RunN.1` | Run name e.g. `Run_MRSPHILLIPS2_2021-01-22_12-20_3553` | **primary `batch_run_id`** |
| `RUND.1` / `RUNT.1` | run start date (list y,m,d) / time (list h,m,s,cs) | fallback `batch_run_id`, run timestamp |
| `MCHN.1` | instrument name | provenance |
| `HCFG.4` | `...SerialNumber=15104-012;` | instrument serial (parse) |
| `LANE.1` | capillary number | spatial QC |
| `TUBE.1` | well (e.g. `H12`) | provenance / plate layout |
| `CTNM.1` / `CTID.1` | container/plate name | provenance / alt grouping |
| `SpNm.1` | sample name typed at machine | cross-check vs `unique_id` |
| `DyeN.1..5` | dye names (6-FAM,VIC,NED,PET,LIZ) | channel validation |
| `OfSc.1` | off-scale scan list | saturation QC (already read) |
| `Satd.1` | camera-saturated scan list | saturation QC (currently unused) |

---

## 1. FSA metadata parsing

### 1a. New internal helper — `R/fsa_metadata.R`
```r
parse_fsa_metadata <- function(abif) {
  d <- abif$Data
  list(
    run_id           = .fsa_run_id(d),          # RunN.1 or fallback
    run_name         = .chr(d$RunN.1),
    run_datetime     = .fsa_datetime(d),        # ISO "2021-01-22 12:20:06" or NA
    instrument       = .chr(d$MCHN.1),
    instrument_serial= .parse_serial(d$HCFG.4), # regex SerialNumber=([^;]+)
    plate            = .chr(d$CTNM.1),
    well             = .chr(d$TUBE.1),
    capillary        = .num(d$LANE.1),
    sample_name      = .chr(d$SpNm.1),
    dye_names        = .dye_names(d),            # named chr vector DyeN.1..n
    n_offscale       = length(d$OfSc.1 %||% integer()),
    n_saturated      = length(d$Satd.1 %||% integer())
  )
}
```
Helpers (same file, internal):
- `.fsa_run_id(d)`: `if (!is.null(d$RunN.1)) trimws(d$RunN.1[1])` else `.fsa_run_id_fallback(d)`.
- `.fsa_run_id_fallback(d)`: build from `RUND.1` (list y,m,d) + `RUNT.1` (list h,m,s) via `sprintf("%04d%02d%02d_%02d%02d%02d", ...)` — **zero-padded** (fixes the lossy `paste(collapse="")` users do today). Return `NA_character_` if either missing.
- `.fsa_datetime(d)`: `sprintf("%04d-%02d-%02d %02d:%02d:%02d", ...)` from RUND.1/RUNT.1, else NA.
- `.chr/.num`: coerce first element or NA; `%||%` null-coalesce.
- All wrapped so a missing tag yields NA, never an error.

### 1b. Class field — [R/classes.R](R/classes.R)
Add to public fields (near `fsa = NULL`):
```r
#' @field fsa_metadata Parsed ABIF header fields (list); populated by read_fsa()
fsa_metadata = NULL,
```

### 1c. Populate at read — [R/constructors.R](R/constructors.R) `read_fsa()`
After line 36 (`fragments_list[[i]] <- fragments$new(...)`):
```r
fragments_list[[i]]$fsa_metadata <- parse_fsa_metadata(fragments_list[[i]]$fsa)
```

### 1d. New exported accessor — [R/extract_functions.R](R/extract_functions.R)
```r
extract_fsa_metadata <- function(fragments_list)  # one row per sample, flattens dye_names to dye_1..n
```
`@export`, `@examples`, used by users + shiny. Returns `data.frame` with `unique_id` + scalar columns (drop the dye vector or expand to `dye_ladder`, `dye_signal`).

---

## 2. Auto-derived `batch_run_id` (default on, opt-out)

### 2a. Config — [R/config.R](R/config.R) + `inst/extdata/trace_config.yaml`
Add:
```yaml
auto_batch_run_id: true        # FSA input: derive batch_run_id from FSA run_id
batch_run_id_tag: "RunN"       # "RunN" (default) or "RUND_RUNT"
```
Add validation in config.R (logical / character checks alongside existing).

### 2b. New internal step — `R/fsa_metadata.R`
```r
set_batch_run_id_from_fsa <- function(fragments_list, config) {
  output <- trace_output$new("set_batch_run_id_from_fsa")
  mismatches <- character()
  for (x in fragments_list) {
    if (!identical(x$input_method, "fsa") || is.null(x$fsa_metadata)) next
    fsa_run <- x$fsa_metadata$run_id
    if (is.na(fsa_run)) next
    if (!is.na(x$batch_run_id) && x$batch_run_id != fsa_run)
      mismatches <- c(mismatches, sprintf("%s (metadata='%s' fsa='%s')",
                                           x$unique_id, x$batch_run_id, fsa_run))
    x$batch_run_id <- fsa_run          # FSA is authoritative
  }
  if (length(mismatches))
    output$set_status("warning", paste0(
      "User-supplied batch_run_id disagreed with the FSA run id (FSA used). ",
      "Check for run mix-ups:\n  ", paste(mismatches, collapse = "\n  ")))
  output
}
```
Optional `sample_name` cross-check (lower priority): warn when `fsa_metadata$sample_name` is set and clearly disjoint from `unique_id`.

### 2c. Wire into `trace()` — [R/main.R](R/main.R)
After the `add_metadata` block (after line 144), before `trace_fsa`:
```r
if (input_type == "fsa" && isTRUE(config$auto_batch_run_id)) {
  s <- set_batch_run_id_from_fsa(fragments_list, config); s$print_status()
}
```
Runs whether or not metadata was supplied (so plotting-only users benefit).

### 2d. Resolve TODO at [main.R:135](R/main.R#L135)
Before pipeline, if `config$correction %in% c("batch","repeat")` or `isTRUE(config$grouped)` and the required metadata fields are all NA across samples **and** can't be FSA-derived (non-fsa input), `stop()` with a specific message naming the missing field(s).

### 2e. Backward-compat note (NEWS.md)
Default flips run grouping to FSA-derived. Set `auto_batch_run_id = FALSE` for the old behavior. `batch_sample_id` / `batch_sample_modal_repeat` are unchanged (not FSA-derivable).

---

## 3. Metadata template generator

### New exported function — `R/metadata_template.R`
```r
generate_metadata_template <- function(input, output_csv = NULL)
```
- `input`: a `read_fsa()` list, a character vector of `.fsa` paths, or a directory.
  - If paths/dir → `read_fsa()` internally (only headers needed; full read is fine).
- Returns `data.frame`, one row per sample, ordered by `batch_run_id` then `well`:

| column | value |
|--------|-------|
| `unique_id` | filename |
| `metrics_group_id` | `NA` (fill) |
| `metrics_baseline_control` | `NA` (fill) |
| `batch_run_id` | **prefilled** from `fsa_metadata$run_id` |
| `batch_sample_id` | `NA` (fill) |
| `batch_sample_modal_repeat` | `NA` (fill) |
| `fsa_run_date` `fsa_well` `fsa_plate` `fsa_instrument` `fsa_capillary` `fsa_sample_name` | read-only provenance |

- The six canonical columns match `add_metadata()` ([add_metadata.R:46-49](R/add_metadata.R#L46)); the `fsa_*` columns are ignored by `add_metadata()` (it reads only the 6), so they're safe passthrough.
- If `output_csv` given, `write.csv(..., row.names = FALSE, na = "")`.
- Docs: "recommended starting point — fill the blank columns, then pass to `trace(metadata_data.frame = ...)`. `batch_run_id` is pre-filled from the file; do not edit unless you know a run was mislabeled."

---

## 4. `qc_report()` — sample-QC panel

### New exported function — `R/qc_report.R`
```r
qc_report <- function(fragments_list, config_file = NULL, ...)
```
Run **after** the pipeline (needs `ladder_df`, `peak_table_df`, allele). Returns one row per sample:

**Provenance** (from `fsa_metadata`): `unique_id, batch_run_id, fsa_run_date, instrument, capillary, well, plate`.

**Metrics + per-metric flags** (thresholds from config, overridable via `...`):

| metric | source | flag condition | catches |
|--------|--------|----------------|---------|
| `ladder_avg_rsq` | `ladder_fit_cor()` mean | — | (report) |
| `ladder_min_rsq` | `ladder_fit_cor()` min | `< qc_min_rsq` (0.998) | **broken ladder** (P-21-207 had avg ok, min=0.90) |
| `n_peaks` | `nrow(peak_table_df)` | `< qc_min_peaks` (5) | **failed/empty** sample |
| `modal_signal` | `get_allele_peak()$allele_signal` | `< qc_min_modal_signal` (e.g. 500) | low-signal sample |
| `modal_saturated` | modal `off_scale` or `>= qc_saturation_ceiling` (32000) | TRUE | over-loaded; peak ratios unreliable |
| `saturation_in_window` | count `off_scale` peaks within `qc_window` bp of modal | `> 0` | saturation hitting analysis region |
| `modal_prominence` | modal signal ÷ mean(±1 neighbor) | `< qc_prominence_min` (1.0) | suspect modal call (neighbor taller) |

**Output columns**: the metrics above, plus
- `qc_flags`: `;`-joined names of tripped checks (`""` if none)
- `qc_pass`: `length(flags)==0`

**Size-awareness note (documented, v1 behavior):** `modal_prominence` declines with repeat size (ρ=−0.50). v1 uses an absolute floor (`< 1.0` = neighbor taller than the called modal → genuinely suspect). Document that long-repeat cohorts have intrinsically low prominence and should rely on visual checks / a lower threshold; full residual-vs-size normalization is a v2 enhancement.

### Config additions ([R/config.R](R/config.R) + yaml)
```yaml
qc_min_rsq: 0.998
qc_min_peaks: 5
qc_min_modal_signal: 500
qc_saturation_ceiling: 32000
qc_prominence_min: 1.0
qc_window: 25
```

### Refactor opportunity
Factor the prominence + in-window-saturation calc into a small internal helper reused by `qc_report()` and any plotting. Reuse `ladder_fit_cor()` (already backing `extract_ladder_summary()`).

---

## 5. Smaller fixes (fold in)
- `add_metadata()` ([add_metadata.R:95-110](R/add_metadata.R#L95)): the `lapply` mutates R6 by reference but its result is discarded and `output` returned. Add a comment documenting the by-reference contract (or return the list) so it's not mistaken for a bug.
- Optional dye/channel validation in `find_ladders()`: before using `config$ladder_channel`/`signal_channel`, compare the channel's `DyeN.*` to expected (LIZ for ladder) and warn on mismatch instead of silently using whatever is at `DATA.105`/`DATA.1`.

---

## 6. Tests (`tests/testthat/`)
- `test-fsa_metadata.R`: parse a bundled `.fsa`; assert `run_id`, `well`, `capillary`, datetime format; assert NA (not error) when tags absent (synthesize an abif with tags removed). Test `RunN` vs `RUND_RUNT` fallback.
- `test-add_metadata.R` (extend): two files from different runs → `batch_run_id` auto-set; supplying a wrong `batch_run_id` emits the mismatch warning; `auto_batch_run_id = FALSE` preserves user value.
- `test-metadata_template.R`: returns the 6 canonical cols + `fsa_*`; `batch_run_id` prefilled; round-trips through `add_metadata()` without warning about missing columns.
- `test-qc_report.R`: a known-good sample → `qc_pass == TRUE`; a few-peaks / broken-ladder fixture trips `n_peaks` / `ladder_min_rsq` respectively.

**Fixtures:** need a real `.fsa` with headers. Options: reuse `seqinr` example (`system.file("abif/...", package="seqinr")`) for single-file parse (verify it carries `RUND`/`RunN`; if not, only the fallback/NA paths are testable with it), and add 1–2 small real `.fsa` (one good, one bad) to `tests/testthat/fixtures/` for the panel + auto-batch tests.

---

## 7. Build order & rationale
1. **§1 parsing** — foundation, zero behavior change, ship + test first.
2. **§3 template** — fast win, depends only on §1, immediately useful in scripts + shiny.
3. **§2 auto-batch** — the safety change; riskiest (default flip) so it lands after §1/§3 are solid, with the opt-out + NEWS entry.
4. **§4 qc_report** — largest new surface; depends on §1 and the pipeline.
5. **§5 fixes** — alongside.

## 8. Docs / packaging
- `devtools::document()` after adding `@export`/`@field` (updates NAMESPACE + man/).
- NEWS.md: entries for the three new exported functions + the `auto_batch_run_id` default-on note.
- pkgdown: add the new functions to `_pkgdown.yml` reference; consider a "QC & batch metadata" vignette.
- Add `^REFACTOR_SPEC\.md$` to `.Rbuildignore` (this file is a dev doc, not package content).
