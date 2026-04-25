# Suppress R CMD check notes for symbols exported to cluster workers or used
# in package-internal control flow where static analysis cannot resolve them.
utils::globalVariables(c("parallel_qr", "pkg_lib_paths"))
