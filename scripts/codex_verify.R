dir.create(".tmp", showWarnings = FALSE)
dir.create(".r-lib", showWarnings = FALSE)

tmp_dir <- normalizePath(".tmp", winslash = "/", mustWork = TRUE)
lib_dir <- normalizePath(".r-lib", winslash = "/", mustWork = TRUE)

Sys.setenv(TMP = tmp_dir, TEMP = tmp_dir, TMPDIR = tmp_dir)
.libPaths(c(lib_dir, .libPaths()))

files <- c(
  "R/solver_utils.R",
  "R/get_B.R",
  "R/blockfit_solve.R",
  "tests/testthat/test_correlation_structure.R",
  "tests/codex_active_set_smoke.R"
)

for (f in files) {
  cat("Parsing:", f, "\n")
  parse(file = f)
}

cat("Checking installed package\n")
pkg_path <- find.package("lgspline", lib.loc = lib_dir, quiet = TRUE)

if (length(pkg_path) == 0) {
  stop("Install lgspline into .r-lib first; scripts/codex_setup.R prints the command.")
}

cat("Running focused smoke\n")
source("tests/codex_active_set_smoke.R")

cat("Focused verification passed\n")
cat("Next: testthat::test_file('tests/testthat/test_correlation_structure.R')\n")
