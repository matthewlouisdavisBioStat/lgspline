dir.create(".tmp", showWarnings = FALSE)
dir.create(".r-lib", showWarnings = FALSE)

tmp_dir <- normalizePath(".tmp", winslash = "/", mustWork = TRUE)
lib_dir <- normalizePath(".r-lib", winslash = "/", mustWork = TRUE)

Sys.setenv(TMP = tmp_dir, TEMP = tmp_dir, TMPDIR = tmp_dir)
.libPaths(c(lib_dir, .libPaths()))

required <- c("testthat", "quadprog")
installed <- rownames(installed.packages(lib.loc = .libPaths()))
missing <- setdiff(required, installed)

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org", lib = lib_dir)
}

cat("tmp:", tmp_dir, "\n")
cat("lib:", lib_dir, "\n")
cat("install package:\n")
cat(
  "$env:TMP=(Resolve-Path '.\\.tmp').Path; ",
  "$env:TEMP=$env:TMP; $env:TMPDIR=$env:TMP; ",
  "& '", file.path(R.home(), "bin", "Rscript.exe"), "' ",
  "-e \".libPaths(c('.r-lib', .libPaths())); ",
  "install.packages('.', repos=NULL, type='source', lib='.r-lib')\"\n",
  sep = ""
)
cat("then: Rscript scripts/codex_verify.R\n")
