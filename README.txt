Ky-Fan covariance permutation tests

This repository contains self-contained R Markdown programs for a two-sample
covariance testing study. The simulations examine null rejection probabilities.
The procedures and their calibration are described in the corresponding Rmd.

Files

simulation_reproducibility/01_reproduce_tables.Rmd
  Rebuild saved manuscript tables or run the original simulation design.
  The default mode is "archived": it reconstructs embedded saved results
  without running new Monte Carlo experiments. Select mode="rerun" to simulate.

simulation_reproducibility/02_compare_methods_H0.Rmd
  Compare fixed Ky-Fan orders k=1 and k=3 with LC, CLX, CZZW and WZZ under H0.
  Defaults: 1,000 repetitions and 1,000 permutation/multiplier draws.
  The six methods share each generated sample. CZZW and WZZ retain their
  historical random-seed tokens, so changing the labels does not change draws.

simulation_reproducibility/04_plot_appendix_C.Rmd
  Reproduce eight Appendix C figures from the saved comparison CSV.
  The file recomputes and checks pointwise 95% Wilson intervals before drawing.
  It uses m=150, p=10 and 40, four covariance structures, and 15 distribution
  pairs: 720 results across the six methods. It does not run simulations.

manuscript/supplementary/comparison_summary.csv
  Saved six-method summary: 180 simulation settings x 6 methods = 1,080 rows.
  This file is required for the default figure-reproduction workflow.

Requirements and execution

Install R, rmarkdown and knitr. Rendering an HTML report also needs Pandoc,
which is available in RStudio. Statistical calculations and PDF plotting use
base R. RhpcBLASctl is optional; when available, it limits BLAS threads.
The simulation files allow at most seven workers.

Open a selected Rmd in RStudio and choose Knit, or run from the repository root:

  rmarkdown::render("simulation_reproducibility/01_reproduce_tables.Rmd")
  rmarkdown::render("simulation_reproducibility/04_plot_appendix_C.Rmd")

The first command rebuilds saved tables in the default archived mode.
The second redraws the appendix figures. Outputs go to simulation_results/
and manuscript/figures/, respectively. These generated files are ignored by Git.

For a short comparison check before a full run:

  rmarkdown::render("simulation_reproducibility/02_compare_methods_H0.Rmd",
    params=list(profile="smoke", R=3, B=39, workers=1))

The smoke profile checks that the program runs; its output is not an estimate
at the precision used in the paper. Knit with default parameters performs the
full comparison. Check each Rmd's parameters before starting a simulation.

Method abbreviations

LC: Li and Chen (2012).
CLX: Cai, Liu and Xia (2013).
CZZW: Chang, Zhou, Zhou and Wang (2017).
WZZ: Wang, Zhu and Zhang (2024).

Provenance

The comparison CSV is the six-method export used in the manuscript. The source
run also included SY, which has been removed from this export and the current
comparison code. Only method names were changed to CZZW/WZZ in the saved data;
the numerical results were retained. The Rmd files explain model applicability
and distinguish saved results from new Monte Carlo output.

New simulation runs save configurations, random-seed information, and checkpoints.
The comparison program also preserves an input Rmd snapshot.
Code or parameter changes can create a new run identifier. To resume a run,
keep the matching input file and configuration. Archived-mode table reconstruction
only writes the rebuilt table outputs.
The large historical runs and unpublished manuscript are kept separately;
they are not needed to run these self-contained Rmd files.
