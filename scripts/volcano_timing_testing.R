## Volcano timing benchmark for lgspline
#
#  Goal:
#   revisit the volcano timing benchmark now that parallel_qr is available.
#   The script benchmarks the same volcano surface fit with and without
#   parallel_qr across increasing partition counts, up to K + 1 = 250.

.libPaths(c(".r-lib", .libPaths()))

suppressPackageStartupMessages({
  library(lgspline)
  library(parallel)
})


## Data

make_volcano_long <- function() {
  data("volcano")
  out <- Reduce("rbind", lapply(seq_len(nrow(volcano)), function(i) {
    t(sapply(seq_len(ncol(volcano)), function(j) {
      c(i, j, volcano[i, j])
    }))
  }))
  colnames(out) <- c("Length", "Width", "Height")
  out
}


## Benchmark helper

run_one_volcano_fit <- function(volcano_long,
                                K,
                                parallel_qr,
                                cl,
                                opt = FALSE,
                                include_quadratic_interactions = TRUE,
                                include_quartic_terms = TRUE,
                                include_constrain_second_deriv = FALSE,
                                reps = 2L) {

  times <- numeric(reps)

  for (r in seq_len(reps)) {
    gc()
    set.seed(1234 + r)

    fit_time <- system.time({
      fit <- lgspline(
        volcano_long[, c("Length", "Width")],
        volcano_long[, "Height"],
        include_quadratic_interactions = include_quadratic_interactions,
        include_quartic_terms = include_quartic_terms,
        include_constrain_second_deriv = include_constrain_second_deriv,
        unique_penalty_per_partition = FALSE,
        opt = opt,
        return_G = FALSE,
        return_Ghalf = FALSE,
        return_U = FALSE,
        return_varcovmat = FALSE,
        estimate_dispersion = FALSE,
        K = K,
        cl = cl,
        parallel_eigen = FALSE,
        parallel_trace = FALSE,
        parallel_aga = FALSE,
        parallel_matmult = FALSE,
        parallel_qr = parallel_qr,
        parallel_unconstrained = FALSE,
        parallel_find_neighbors = FALSE,
        parallel_penalty = FALSE,
        parallel_make_constraint = FALSE,
        include_warnings = TRUE,
        verbose = FALSE
      )
    })

    stopifnot(inherits(fit, "lgspline"))
    times[r] <- unname(fit_time["elapsed"])
  }

  data.frame(
    K = K,
    partitions = K + 1L,
    parallel_qr = parallel_qr,
    reps = reps,
    elapsed_mean = mean(times),
    elapsed_median = median(times),
    elapsed_min = min(times),
    elapsed_max = max(times),
    stringsAsFactors = FALSE
  )
}


## Main benchmark

benchmark_volcano_parallel_qr <- function(K_grid = c(24L, 49L, 99L, 249L),
                                          reps = 2L,
                                          opt = FALSE,
                                          n_workers = NULL) {

  volcano_long <- make_volcano_long()

  if (is.null(n_workers)) {
    det <- parallel::detectCores(logical = FALSE)
    if (is.na(det) || det < 2L) {
      det <- parallel::detectCores()
    }
    if (is.na(det) || det < 2L) {
      n_workers <- 2L
    } else {
      n_workers <- min(4L, det)
    }
  }

  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  ## Warm start to avoid first-call overhead dominating the smallest benchmark.
  invisible(run_one_volcano_fit(
    volcano_long = volcano_long,
    K = K_grid[[1]],
    parallel_qr = FALSE,
    cl = cl,
    opt = opt,
    reps = 1L
  ))

  results <- do.call(
    rbind,
    lapply(K_grid, function(K) {
      rbind(
        run_one_volcano_fit(
          volcano_long = volcano_long,
          K = K,
          parallel_qr = FALSE,
          cl = cl,
          opt = opt,
          reps = reps
        ),
        run_one_volcano_fit(
          volcano_long = volcano_long,
          K = K,
          parallel_qr = TRUE,
          cl = cl,
          opt = opt,
          reps = reps
        )
      )
    })
  )

  results$workers <- n_workers
  results$opt <- opt

  speed_summary <- do.call(
    rbind,
    lapply(split(results, results$partitions), function(df_part) {
      serial_row <- df_part[!df_part$parallel_qr, , drop = FALSE]
      qr_row <- df_part[df_part$parallel_qr, , drop = FALSE]

      data.frame(
        partitions = serial_row$partitions,
        K = serial_row$K,
        workers = serial_row$workers,
        opt = serial_row$opt,
        serial_median = serial_row$elapsed_median,
        parallel_qr_median = qr_row$elapsed_median,
        speedup = serial_row$elapsed_median / qr_row$elapsed_median,
        percent_change = 100 * (
          qr_row$elapsed_median - serial_row$elapsed_median
        ) / serial_row$elapsed_median,
        stringsAsFactors = FALSE
      )
    })
  )

  write.csv(results,
            file = "scripts/volcano_timing_results.csv",
            row.names = FALSE)
  write.csv(speed_summary,
            file = "scripts/volcano_timing_summary.csv",
            row.names = FALSE)

  list(
    results = results,
    summary = speed_summary
  )
}


## Run

if (identical(Sys.getenv("RUN_VOLCANO_TIMING_DEFAULT"), "true")) {
  bench <- benchmark_volcano_parallel_qr()

  cat("\nVolcano Timing Results\n")
  print(bench$results)

  cat("\nVolcano Timing Summary\n")
  print(bench$summary)
}
