# Source the chip_lod R port. Usage: source("r/load.R") from the repo root, or with chdir = TRUE.
local({
  here <- if (!is.null(sys.frame(1)$ofile)) dirname(sys.frame(1)$ofile) else "r"
  for (f in c("rng.R", "model.R", "lod.R", "io.R", "simulate.R", "analyses.R", "freq_analyses.R")) {
    sys.source(file.path(here, "R", f), envir = globalenv())
  }
})
