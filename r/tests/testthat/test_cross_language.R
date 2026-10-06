# Byte-for-byte comparison of results/python and results/r (sha256 via the system tool, plus raw bytes).
sha256_file <- function(path) {
  tool <- Sys.which(c("sha256sum", "shasum"))
  tool <- tool[nzchar(tool)][1]
  args <- if (basename(tool) == "shasum") c("-a", "256", shQuote(path)) else shQuote(path)
  strsplit(system2(tool, args, stdout = TRUE), " ")[[1]][1]
}

test_that("Python and R result files are byte-identical", {
  py <- file.path(ROOT, "results", "python")
  rr <- file.path(ROOT, "results", "r")
  files <- sort(union(list.files(py), list.files(rr)))
  skip_if(length(files) == 0, "no results yet: run the Python and R pipelines first")
  for (f in files) {
    expect_true(file.exists(file.path(py, f)), info = paste("missing in python:", f))
    expect_true(file.exists(file.path(rr, f)), info = paste("missing in r:", f))
    if (file.exists(file.path(py, f)) && file.exists(file.path(rr, f))) {
      expect_identical(sha256_file(file.path(rr, f)), sha256_file(file.path(py, f)), info = f)
      a <- readBin(file.path(py, f), "raw", file.size(file.path(py, f)))
      b <- readBin(file.path(rr, f), "raw", file.size(file.path(rr, f)))
      expect_identical(b, a, info = f)
    }
  }
})
