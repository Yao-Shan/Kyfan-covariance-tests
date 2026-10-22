# Ky–Fan covariance tests

Simulation code for two-sample covariance equality using fixed-order Ky–Fan
permutation tests. Ordinary R scripts run the experiments; R Markdown reports
are optional.

## Files

```text
simulation/   R scripts for experiments and figures
R/            Data generators, tests, simulation runners, and plotting functions
data/         Saved results for reference comparisons and Appendix C
reports/      Optional R Markdown reports
results/      Generated outputs (created when needed; excluded from Git)
```

| Script in `simulation/` | Purpose |
|:--|:--|
| `01_reproduce_tables.R` | Simulate the Ky–Fan results for Tables 1–3 |
| `02_compare_methods_H0.R` | Compare Ky–Fan orders 1 and 3 with LC, CLX, CZZW, and WZZ under covariance equality |
| `04_plot_appendix_C.R` | Validate saved comparison results and draw eight Wilson interval figures |

The corresponding files in `reports/` explain the designs and display results.
`R/simulation_core.R` contains the shared data generators and Ky–Fan tests;
`R/entrywise_tests.R` contains LC, CLX, and CZZW; `R/wzz_test.R` contains WZZ.
The remaining modules handle experiment scheduling, output, and plotting.

## Run

The scripts require only base R. They were checked with R 4.5.1. Open
`Kyfan-covariance-tests.Rproj` in RStudio and source a script, or run these
commands from the repository directory:

```sh
Rscript simulation/01_reproduce_tables.R --smoke
Rscript simulation/02_compare_methods_H0.R --smoke
Rscript simulation/04_plot_appendix_C.R
```

`--smoke` uses a small design with three repetitions, 39 permutations or multiplier
draws, and one worker. It checks execution, not test size. The figure script reads
saved results and performs no simulation. Absolute script paths also work when
launching from another directory.

Edit the `params` list at the top of a script to change settings. Without
`--smoke`, the first two scripts run their full default designs:

| Experiment | Repetitions `R` | Draws `B` | Workers | Summary rows |
|:--|--:|--:|--:|--:|
| Tables 1–3 | 1,000 | 1,000 | 1 | 360 |
| Six-method comparison | 1,000 | 1,000 | 7 | 1,080 |

The design uses `(p,m) = (10,50), (10,150), (40,150)`, four covariance structures,
and all 15 pairs of the N, T, E, G, and L distributions. The nominal level is 0.05.
Appendix C uses only `m=150`, giving 720 cells.
Worker counts are capped at seven, with one BLAS thread per worker. The table
runner uses serial execution on Windows.

New outputs go to `results/tables/`, `results/comparison/`, and
`results/figures/`. Simulation runs save summaries,
replicate results, checkpoints, configuration, and a runnable source snapshot.
Matching configurations can resume interrupted runs. Source changes, including
comments or switching between a script and a report, produce different run IDs.
Use a run's original source snapshot to resume that version.

## Saved data

| File in `data/` | Contents |
|:--|:--|
| `comparison_summary.csv` | 1,080 six-method results, with counts, MCSEs, and 95% Wilson intervals; `R=B=1000` |
| `manuscript_reference.csv` | 360 saved Table 1–3 rates from the completed October 6, 2026 rerun |
| `table_archive.csv` | 1,020 saved table cells: the 360 main results plus 660 historical extension results |

These files contain estimates, not individual observations. Fresh simulation in
scripts 01 and 02 does not require saved CSVs. Script 01 accepts an optional fixed
reference and an explicit `mode="archived"` for viewing saved tables.
`manuscript_reference.csv` is an optional fixed reference for Tables 1–3;
script 04 reads `comparison_summary.csv` to draw Appendix C figures.

## Reports

HTML reports additionally require `rmarkdown`, `knitr`, and Pandoc (included with
RStudio). For example, from the repository directory:

```r
rmarkdown::render(
  "reports/02_compare_methods_H0.Rmd",
  params = list(profile = "smoke", R = 3, B = 39, workers = 1)
)
rmarkdown::render("reports/04_plot_appendix_C.Rmd")
```

Reports and scripts use the same workflow functions. Relative input and output
paths resolve from the repository directory. PDF reports additionally require
a LaTeX installation; the figure script itself uses R's built-in PDF device.

## Interpretation

LC denotes Li–Chen, CLX Cai–Liu–Xia, CZZW Chang–Zhou–Zhou–Wang, and WZZ
Wang–Zhu–Zhang. Each test retains its own calibration. WZZ follows the authors'
public 2022 manuscript corresponding to their 2024 article; this code does not
claim a formula-by-formula audit of the final journal version.

The experiments measure rejection rates under covariance equality. Some settings
fall outside competing methods' sufficient conditions. See the reports for
these qualifications.
Failures are recorded separately and must not be treated as non-rejections.
