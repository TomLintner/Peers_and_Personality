# Analysis code: peers and personality in two classroom network panels

Everything needed to reproduce the reported results, in run order. Exploratory model-building steps,
timing and scaling tests, inspection scripts and manuscript-building utilities are not included.

## Running it

```r
setwd("<path>/code_bundle")
source("00_setup.R")          # installs any missing packages, creates work/ and output/
source("R/03_build_sample1.R")
```

Every script in `R/` sources `00_setup.R` itself, so each can be run on its own once the scripts it
depends on have run. Paths come from `00_setup.R`: `RAW` for input data, `WORK` for intermediates,
`OUT` for fits, chains, tables and figures. Seed 2012 throughout.

`00_setup.R` installs RSiena, lavaan, ggplot2, ggridges, svglite, openxlsx, officer, flextable, sna
and png from CRAN. **RSienaTest** supplies `sienaBayes`, which RSiena does not; it is on R-Forge,
source only, so it needs a compiler toolchain (Rtools on Windows). Without it, scripts 20 to 33 will
not run and the rest will.

## Data and what is included

```
work/objects/   1.2 MB, 10 files  - built siena objects, sienaBayes setups, siena08 results
work/chains/    385 MB, 29 files  - bayes12h_stage1-9.rds, ozbayes_stage1-20.rds
work/mom/       empty             - the converged method-of-moments fits, see below
data_raw/       empty             - raw survey files and sociometric matrices, not distributed
output/         created on first run
```

**Included: the posterior chains.** `work/chains/` holds the 29 sienaBayes stage files, 2,650 draws
for Sample 1 across nine stages and 6,000 for Sample 2 across twenty. These are parameter draws, not
individual-level data, and they are enough to regenerate every reported table and figure without
re-estimating. Scripts 30, 31, 50, 60, 61, 62, 63, 64 and 65 need only these plus `work/objects/`,
and run in seconds to minutes.

**Not included: the raw data.** Individual survey responses and sociometric matrices are withheld
because classroom matrices identify individuals within classes of 13 to 32 students even without
names, and consent did not cover public release. Available under a data-use agreement; see the
article's Open Science Statement. Scripts 01 to 11 need them.

**Not included: the two method-of-moments fits.** `gacr_final_ans.rds` and `oz_final_ans.rds` are
485 MB and 339 MB, because `n3 = 40000` stores the Phase 3 statistics by iteration. Script 12
(goodness of fit) expects them in `work/mom/`; scripts 10 and 11 regenerate them from the raw data.
They hold no individual-level data and can be added to `work/mom/` if wanted.

**Two intermediates were never persisted:** `gacr_cfa2.rds` (written by 02, read by 03) and
`oz_bins.rds` (written by 04, read by 06) existed only in a temporary directory. Running 01 to 06 in
order recreates them; the saved copies are gone.

## Run order

### Data preparation
| Script | Output |
|---|---|
| `01_bins_sample1.R` | `gacr_bfi_bins.csv`, `gacr_bins.rds` |
| `02_cfa_sample1.R` | `gacr_cfa2.rds` — longitudinal CFA per domain, lavaan |
| `03_build_sample1.R` | `gacr_siena.rds` — RSiena group object, 27 classes |
| `04_bins_sample2.R` | `oz_bfi_bins.csv`, `oz_bins.rds` |
| `05_bins_sample2_refine.R` | revised `oz_bfi_bins.csv` after the reliability screen |
| `06_build_sample2.R` | `oz_siena.rds` — RSiena group object, 18 classes |

### Method of moments
| Script | Output |
|---|---|
| `10_mom_sample1.R` | `gacr_final_ans.rds`, 205 parameters, tconv.max 0.2456 |
| `11_mom_sample2.R` | `oz_final_ans.rds`, 151 parameters, tconv.max 0.2210 |

`n3` is raised in stages to 3,000, 10,000, 20,000 and 40,000 because `tconv.max` has a floor near
`sqrt(p / n3)`: at p = 205 and n3 = 3,000 that floor is 0.26, so 0.25 is unreachable. `thetaBound = 50`
keeps the class-specific rate parameters bounded.

### Goodness of fit
| Script | Output |
|---|---|
| `12_gof_phase3.R` | `mom_s1_deps.rds`, `mom_s2_deps.rds`, `gof_s{1,2}_table.csv` |
| `13_gof_summarise.R` | `gof_summary.csv`, `gof_classroom_table.csv` |

`sienaGOF` needs the simulated networks, which the stored fits do not carry, so 12 re-runs Phase 3
only (`nsub = 0`, `prevAns` = the saved fit, `returnDeps = TRUE`) rather than re-estimating.

### Bayesian multilevel estimation
| Script | Output |
|---|---|
| `20_bayes_sample1.R` | `gacr_bayes_full.rds` — first long run |
| `21_bayes_sample1_stages.R` | `bayes12h_stage1-9.rds` chained via `prevBayes` |
| `22_bayes_sample1_continue.R` | further stages, `bayes12h_accumulated.rds` |
| `23_bayes_sample1_tail.R` | retained-draw summaries after burn-in |
| `24_bayes_setup_sample2.R` | `oz_bayes_init.rds` |
| `25_bayes_rebuild_sample2.R` | `oz_bayes_group_fixed.rds`, `oz_bayes_setup.rds` |
| `26_bayes_sample2.R` | `ozbayes_stage1-20.rds`, per-stage diagnostics |
| `27_bayes_sample2_continue.R` | further stages |
| `28_bayes_sample2_tail.R` | retained-draw summaries after burn-in |

Sample 1 reaches 2,650 draws (900 discarded), Sample 2 6,000 (1,200 discarded). 19 parameters vary
across classes, 30 are constant. Script 25 exists because the Sample 2 data objects carried class
`"sienadata"` where Sample 1 carried `"siena"`, which `checkImpossibleChanges` rejects, and because
31 forbidden changes had to be set to missing.

### Diagnostics
| Script | Output |
|---|---|
| `30_diag_cumulative.R` | `bayes_cumulative_diag.csv` — R-hat and ESS by stage |
| `31_diag_drift.R` | first-half vs second-half drift in posterior SD units |
| `32_diag_mds_sample1.R` | `defending_victims_mds.png` |
| `33_diag_mds_sample2.R` | `moment_of_change_mds.png` |

The drift statistic, not R-hat, is what identified the burn-in: R-hat rose while ESS rose for ten
stages of Sample 2 because the chain was still in a transient.

### Wave-1 assortativity
| Script | Output |
|---|---|
| `40_assortativity.R` | `homophily_classroom.csv`, `homophily_dataset.csv` |
| `41_assortativity_tables.R` | `fig_homophily.png/.svg` (Figure 1), Table 3 |

Newman's coefficient centred on a within-class permutation null, 1,000 permutations, pooled with
tie-count weights, bootstrapped across classes, class p values combined by Stouffer's method.

### Between-class heterogeneity
| Script | Output |
|---|---|
| `50_heterogeneity_27classes.R` | `gacr_meta27_table.csv` — `siena08` across the 27 classes |
| `51_meta_two_study.R` | `meta_two_study.csv` — fixed-effect pooling of the two samples |

### Tables and figures
| Script | Output |
|---|---|
| `60_posterior_table_sample1.R` | `bayes_full_table.csv`, `bayes_full_table_mu.csv` |
| `61_posterior_table_sample2.R` | `oz_full_table.csv`, `oz_full_table_mu.csv` |
| `62_manuscript_tables.R` | `paper_table1_samples.csv` … `paper_table3_estimates.csv` |
| `63_figures_ridge.R` | Figures 2 to 6, PNG at 400 dpi and SVG |
| `64_figures_sample1.R` | Sample 1 trace and posterior figures |
| `65_figures_sample2.R` | Sample 2 trace and posterior figures |
| `70_supplemental_descriptives.R` | Tables S1 and S2 |
| `71_supplemental_document.R` | `Supplemental_Material.docx`, Tables S1–S6 and Figures S1–S6 |