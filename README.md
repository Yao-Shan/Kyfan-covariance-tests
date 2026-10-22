# Ky–Fan covariance tests

Code for two-sample covariance equality using fixed-order Ky–Fan permutation
tests. R scripts run the simulations; R Markdown reports display experiments
or test user-supplied data.

## Files

```text
simulation/   R scripts for experiments and figures
R/            Data generators, tests, simulation runners, and plotting functions
data/         Saved results for reference comparisons and Appendix C
reports/      Optional R Markdown reports
results/      Generated outputs (created when needed; excluded from Git)
```

| Purpose | Script in `simulation/` | Report in `reports/` |
|:--|:--|:--|
| Simulate the Ky–Fan results for Tables 1–3 | `kyfan_tables.R` | `kyfan_tables.Rmd` |
| Compare Ky–Fan orders 1 and 3 with LC, CLX, CZZW, and WZZ under covariance equality | `h0_comparison.R` | `h0_comparison.Rmd` |
| Validate saved comparison results and draw eight Wilson interval figures | `appendix_plots.R` | `appendix_plots.Rmd` |
| Test supplied data and report the observed statistic and p-value, with 1,000 permutations by default | — | `kyfan_data_test.Rmd` |

Each simulation script has a matching report that explains the design and
displays results. The supplied-data report accepts two CSV files or R data objects.
`R/simulation_core.R` contains the shared data generators and Ky–Fan tests.
Each comparison test has its own file and function:

| File | Method | Function |
|:--|:--|:--|
| `R/lc_test.R` | Li–Chen (LC) | `method_lc(X, Y)` |
| `R/clx_test.R` | Cai–Liu–Xia (CLX) | `method_clx(X, Y)` |
| `R/czzw_test.R` | Chang–Zhou–Zhou–Wang (CZZW) | `method_czzw(X, Y, B)` |
| `R/wzz_test.R` | Wang–Zhu–Zhang (WZZ) | `method_wzz(X, Y)` |

`R/test_helpers.R` provides shared input checks and CLX/CZZW covariance-entry
products. Load it before LC, CLX, or CZZW when using those files directly;
the comparison script loads all dependencies automatically. Test functions
return a statistic, p-value, status, and diagnostics. The remaining modules
handle experiment scheduling, output, and plotting.

## Run

The scripts require only base R. They were checked with R 4.5.1. Open
`Kyfan-covariance-tests.Rproj` in RStudio and source a script, or run these
commands from the repository directory:

```sh
Rscript simulation/kyfan_tables.R --smoke
Rscript simulation/h0_comparison.R --smoke
Rscript simulation/appendix_plots.R
```

`--smoke` uses a small design with three repetitions, 39 permutations or multiplier
draws, and one worker. It checks execution, not test size. The figure script reads
saved results and performs no simulation. Absolute script paths also work when
launching from another directory.

Edit the `params` list at the top of a script to change settings. Without
`--smoke`, `kyfan_tables.R` and `h0_comparison.R` run their full default designs:

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
`kyfan_tables.R` and `h0_comparison.R` does not require saved CSVs.
`kyfan_tables.R` accepts an optional fixed
reference and an explicit `mode="archived"` for viewing saved tables.
`manuscript_reference.csv` is an optional fixed reference for Tables 1–3;
`appendix_plots.R` reads `comparison_summary.csv` to draw Appendix C figures.

## Reports

HTML reports additionally require `rmarkdown`, `knitr`, and Pandoc (included with
RStudio). For example, from the repository directory:

```r
rmarkdown::render(
  "reports/h0_comparison.Rmd",
  params = list(profile = "smoke", R = 3, B = 39, workers = 1)
)
rmarkdown::render("reports/appendix_plots.Rmd")
```

Simulation reports and scripts use the same workflow functions. Relative input and output
paths resolve from the repository directory. PDF reports additionally require
a LaTeX installation; the figure script itself uses R's built-in PDF device.

## Test supplied data

Use two numeric CSV files with matching measurement columns, or two R matrices
or data frames. Groups must have equal sample sizes. Choose `k` before inspecting
the data; its default is 1. From the repository directory:

```r
rmarkdown::render(
  "reports/kyfan_data_test.Rmd",
  params = list(x_csv = "/path/group_x.csv", y_csv = "/path/group_y.csv",
                k = 1, B = 1000)
)
# Alternatively, replace the CSV parameters with X = X_data, Y = Y_data.
```

The report displays the Ky–Fan norm, its scaled statistic, and the plus-one
p-value. It saves results and permutation statistics under `results/data_test/`.
With no input, it renders an input guide without running a test. Numerical ties
use a relative tolerance of `100 * .Machine$double.eps` times the observed
statistic; the existing simulation code is unchanged.

For calculations without a report, source `R/simulation_core.R` followed by
`R/kyfan_test.R` and call `method_kyfan(X, Y, k = 1, B = 1000)`.
See the report for input requirements and the scope of permutation validity.

## Interpretation

LC denotes Li–Chen, CLX Cai–Liu–Xia, CZZW Chang–Zhou–Zhou–Wang, and WZZ
Wang–Zhu–Zhang. Each test retains its own calibration. WZZ follows the authors'
public 2022 manuscript corresponding to their 2024 article; this code does not
claim a formula-by-formula audit of the final journal version.

The experiments measure rejection rates under covariance equality. Some settings
fall outside competing methods' sufficient conditions. See the reports for
these qualifications.
Failures are recorded separately and must not be treated as non-rejections.
