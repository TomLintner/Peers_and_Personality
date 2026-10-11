## Packages and paths. Every script in R/ sources this first.
## Run from the bundle root: setwd("<...>/code_bundle"); source("R/03_build_sample1.R")

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
RAW  <- file.path(ROOT, "data_raw")   # raw survey files and sociometric matrices; not distributed
WORK <- file.path(ROOT, "work")       # intermediates
OUT  <- file.path(ROOT, "output")     # fits, chains, tables, figures
SC   <- WORK                          # alias the build scripts use
OBJ  <- file.path(WORK, "objects")    # built siena objects, setups, meta results
CHN  <- file.path(WORK, "chains")     # sienaBayes stage files
MOM  <- file.path(WORK, "mom")        # converged method-of-moments fits

for (d in c(WORK, OBJ, CHN, MOM, OUT, file.path(OUT, c("gof",
                                        "tables", "figures", "supplemental"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("installing ", p)
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE) }
  if (!requireNamespace(p, quietly = TRUE)) stop("could not install ", p) }

invisible(lapply(c("RSiena", "lavaan", "ggplot2", "ggridges", "svglite",
                   "openxlsx", "officer", "flextable", "sna", "png"), need))

## RSienaTest provides sienaBayes, which RSiena does not. It is on R-Forge, source only, so a
## compiler toolchain is needed (Rtools on Windows). Without it, scripts 20 to 33 cannot run.
if (!requireNamespace("RSienaTest", quietly = TRUE)) {
  message("installing RSienaTest from R-Forge (source build)")
  try(install.packages("RSienaTest", repos = "http://R-Forge.R-project.org", type = "source"),
      silent = TRUE) }
if (!requireNamespace("RSienaTest", quietly = TRUE))
  message("RSienaTest unavailable: the Bayesian scripts (20-33) will not run; the rest will")

set.seed(2012)
