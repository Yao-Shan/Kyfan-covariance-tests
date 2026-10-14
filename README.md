# Ky–Fan covariance test simulations

Reproducible R workflows for archived tables, null comparisons with LC, CLX,
CZZW and WZZ, mean-shift sensitivity, and Appendix C plots. Keep the directory
structure intact: the four reports load the ten modules in
`simulation_reproducibility/R/` and the supplied CSV inputs.

## Requirements

R, `rmarkdown`, `knitr`, and Pandoc (usually bundled with RStudio). Install the R
packages with `install.packages(c("rmarkdown", "knitr"))`. Calculations and PDF
plots use base R; `RhpcBLASctl` is optional. Simulations use at most seven workers,
with one BLAS/OpenMP thread per worker. Report 01 falls back to serial execution
on Windows; reports 02 and 03 use PSOCK workers and support Windows parallelism.

## Run the reports

Run these commands in R from the repository root:

```r
# Rebuild saved Tables 1–3; no new simulation by default.
rmarkdown::render("simulation_reproducibility/01_reproduce_tables.Rmd")

# Small execution checks for the comparison and mean-shift workflows.
rmarkdown::render("simulation_reproducibility/02_compare_methods_H0.Rmd",
  params = list(profile = "smoke", R = 3, B = 39, workers = 1))
rmarkdown::render("simulation_reproducibility/03_mean_shift_sensitivity.Rmd",
  params = list(profile = "smoke", R = 3, B = 39, workers = 1))

# Validate the saved comparison summary and write eight Appendix C PDFs.
rmarkdown::render("simulation_reproducibility/04_plot_appendix_C.Rmd")
```

Smoke runs check execution; they do not assess statistical calibration. To run
reports 02 and 03 at their full defaults, omit `params`. To simulate Tables 1–3
with report 01, pass `params = list(mode = "rerun")`. Other parameters can be
changed through the YAML header or `params`.

| Report | Default behavior |
| --- | --- |
| `01_reproduce_tables.Rmd` | Archived tables; rerun defaults: R = B = 1000, one worker |
| `02_compare_methods_H0.Rmd` | 180 scenarios, six methods; R = B = 1000, seven workers |
| `03_mean_shift_sensitivity.Rmd` | 180 base scenarios, six mean settings, k = 1 and 3; R = B = 500, seven workers |
| `04_plot_appendix_C.Rmd` | Read the supplied summary and write eight PDFs to `manuscript/figures/`; no simulations |

The mean-shift report compares zero means, a common shift, and unequal means
with difference norms 1 and 2. Unequal means are sensitivity experiments outside
the paper's common-mean null model.

## Results, seeds and resuming

Simulation outputs go to `simulation_results/`, including configurations,
replicate records, summaries, checkpoints, session information and
`source_snapshot/`. The master seed defaults to `20260925`; scenario, replicate
and purpose determine random streams. Changing worker count preserves seeds and
run identity. Checkpoints are saved every ten replications by default.

Run identifiers include scientific parameters and the full contents of loaded
source/reference inputs; even comment edits create a new run identity. Matching
inputs, configuration and R environment allow resuming with `resume = TRUE`.
Use frozen input snapshots to resume older runs. Avoid concurrent identical
configurations. In reports 02 and 03, `COMPLETED.txt` indicates all requested
replications finished; also inspect `failed` counts and replicate `status`.
Those reports record failed computations without counting them as non-rejections.
Report 01 stops on computation errors and does not write a completion marker.

## Saved inputs and provenance

- `simulation_reproducibility/data/table_archive.csv` has 1,020 archive cells:
  360 main-table results from the completed R = B = 1000 rerun
  (`rerun_main_20261006`) and 660 historical extension results that were not rerun.
  The historical random streams were not retained, so exact recovery of those
  experiments is not guaranteed.
- `simulation_reproducibility/data/manuscript_reference.csv` is the shared
  main-table reference for reports 01 and 03. Report 03 recovers reference counts
  and unrounded MCSEs from saved rates using R = 1000.
- `manuscript/supplementary/comparison_summary.csv` contains the saved comparison
  summary used by report 04. Rendering that report validates the selected 720
  Appendix C cells and their Wilson intervals before plotting.

Saved reference values are not new simulations. Method sources and assumptions
are documented in report 02. WZZ follows the authors' public 2022 manuscript,
[arXiv:2212.12338](https://arxiv.org/abs/2212.12338), corresponding to the 2024
journal article; the final journal formulas have not been independently
re-audited here.

## Completed checks

The structural refactor was checked against the original implementation:
small-run statistics, p-values, decisions and seeds matched; serial, two-worker
and resumed runs agreed. All 1,020 archive cells, 17 generated LaTeX tables and
eight appendix plots matched, and all four reports rendered successfully.
These code edits did not repeat the full Monte Carlo experiments or independently
re-audit the original papers' formulas.
