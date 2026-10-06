# Usage: Rscript r/run_all.R [params.json] [outdir]
args <- commandArgs(trailingOnly = TRUE)
params_path <- if (length(args) >= 1) args[1] else "config/params.json"
outdir <- if (length(args) >= 2) args[2] else "results/r"
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
for (f in c("rng.R", "model.R", "lod.R", "io.R", "simulate.R", "analyses.R", "freq_analyses.R")) {
  sys.source(file.path(script_dir, "R", f), envir = globalenv())
}
t0 <- Sys.time()
run_all(read_params(params_path), outdir)
cat(sprintf("wrote results to %s in %.1f s\n", outdir, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
