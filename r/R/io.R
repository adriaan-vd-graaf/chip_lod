# Hand-written, byte-deterministic writers (mirror python/src/chip_lod/io.py).
# integer -> "%d", double -> "%.12e", NA -> "NA" (CSV) / null (JSON). Non-finite doubles are an error.

format_value <- function(x) {
  if (length(x) != 1) stop("format_value expects a scalar")
  if (is.na(x) && !is.nan(x)) return("NA")
  if (is.integer(x)) return(sprintf("%d", x))
  if (is.double(x)) {
    if (!is.finite(x)) stop(sprintf("refusing to write non-finite value %s", format(x)))
    return(sprintf("%.12e", x))
  }
  stop(sprintf("unsupported value type %s", typeof(x)))
}

format_column <- function(col) {
  out <- character(length(col))
  for (i in seq_along(col)) out[i] <- format_value(col[[i]])
  out
}

# df: a data.frame (or named list of equal-length columns) whose column types decide the format.
write_csv <- function(path, df) {
  cols <- lapply(df, format_column)
  lines <- c(paste(names(df), collapse = ","), do.call(paste, c(cols, sep = ",")))
  data <- paste0(paste(lines, collapse = "\n"), "\n")
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(data)), con)
}

.json_scalar <- function(x) {
  if (is.null(x) || (length(x) == 1 && is.na(x) && !is.nan(x))) return("null")
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.character(x)) return(paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x)), '"'))
  format_value(x)
}

.json_lines <- function(obj, indent) {
  pad <- strrep("  ", indent)
  pad_in <- strrep("  ", indent + 1)
  if (is.list(obj)) {
    nm <- names(obj)
    if (!is.null(nm)) {
      if (length(obj) == 0) return("{}")
      items <- character(length(obj))
      for (i in seq_along(obj)) {
        items[i] <- paste0(pad_in, .json_scalar(nm[i]), ": ", .json_lines(obj[[i]], indent + 1))
      }
      return(paste0("{\n", paste(items, collapse = ",\n"), "\n", pad, "}"))
    }
    if (length(obj) == 0) return("[]")
    items <- character(length(obj))
    for (i in seq_along(obj)) items[i] <- paste0(pad_in, .json_lines(obj[[i]], indent + 1))
    return(paste0("[\n", paste(items, collapse = ",\n"), "\n", pad, "]"))
  }
  .json_scalar(obj)
}

to_json <- function(obj) paste0(.json_lines(obj, 0), "\n")

write_json_file <- function(path, obj) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(to_json(obj))), con)
}

read_params <- function(path) jsonlite::fromJSON(path, simplifyVector = TRUE)

sum_loop <- function(xs) {
  s <- 0
  for (x in xs) s <- s + x
  s
}

mean_loop <- function(xs) {
  if (length(xs) == 0) stop("mean of empty sequence")
  sum_loop(xs) / length(xs)
}

# Inverse empirical CDF: element at 1-based index max(1, ceil(q N)) of the sorted values.
quantile_ecdf <- function(xs, q) {
  n <- length(xs)
  if (n == 0) stop("quantile of empty sequence")
  s <- sort(xs)
  idx <- ceiling(q * n)
  if (idx < 1) idx <- 1
  if (idx > n) idx <- n
  s[idx]
}
