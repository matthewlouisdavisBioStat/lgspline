#' lgspline: Lagrangian Multiplier Smoothing Splines
#'
#' @description
#' Allows for common S3 methods including print, summary, coef, plot, and
#' predict, with additional inference methods provided.
#'
#' @details
#' This package implements various methods for working with lgspline models,
#' including printing, summarizing, plotting, and conducting statistical inference.
#'
#' @docType methods
#' @keywords internal
#' @name lgspline-methods
#' @rdname lgspline-methods
#' @aliases lgspline-methods
NULL

#' Print Method for lgspline Objects
#'
#' @description
#' Provides a standard print method for lgspline model objects to display
#' key model characteristics.
#'
#' @param x An lgspline model object
#' @param ... Additional arguments (not used)
#'
#' @export
#' @method print lgspline
print.lgspline <- function(x, ...) {
  cat("Lagrangian Multiplier Smoothing Spline Model\n")
  cat("============================================\n")
  cat("Model Family (Link Function):",
      paste0(paste0(x$family)[1],
             " (",
             paste0(x$family)[2],
             ")", collapse = ""),
      "\n")
  cat("Number of Observations:", x$N, "\n")
  cat("Number of Predictors:", x$q, "\n")
  cat("Number of Partitions:", x$K + 1, "\n")
  cat("Basis Functions per Partition:", x$p, "\n")
  invisible(x)
}

#' Print Method for lgspline Object Summaries
#'
#' @param x An lgspline model object
#' @export
print.summary.lgspline <- function(x, ...) {
  cat("Lagrangian Multiplier Smoothing Spline Model Summary\n")
  cat("====================================================\n")
  cat("Model Family:", x$model_family[[1]], "\n")
  cat("Model Family (Link Function):",
      paste0(paste0(x$model_family)[1],
             " (",
             paste0(x$model_family)[2],
             ")", collapse = ""),
      "\n")
  cat("Observations:", x$observations, "\n")
  cat("Predictors:", x$predictors, "\n")
  cat("Partitions:", x$knots + 1, "\n")
  cat("Basis Functions per Partition:", x$basis_functions, "\n")
  if(length(unlist(x$coefficients)) > 1){
    cat("----------------------------------------------------\n")
    cat("Coefficients and Wald Inference: \n")
    print(x$coefficients)
    cat('\n')
    cat("Dispersion:", x$sigmasq_tilde, "\n")
    cat("Effective degrees of freedom:", x$N - x$trace_XUGX, "\n")
    cat("Critical value for confidence intervals: ", x$cv, "\n")
    cat("----------------------------------------------------\n")
  }
  invisible(x)
}

#' Summary method for lgspline Objects
#'
#' @param object An lgspline model object
#' @export
summary.lgspline <- function(object, ...) {
  ## Create a brief summary
  summary_list <- list(
    model_family = object$family,
    observations = object$N,
    predictors = object$q,
    knots = object$K,
    basis_functions = object$p,
    estimate_dispersion = ifelse(object$estimate_dispersion &
                                 object$sigmasq_tilde != 1,
                                 'Yes',
                                 'No'),
    cv = object$critical_value
  )

  ## Typical summaries for Wald inference, like lm() or glm()
  if(object$return_varcovmat){
    tr <- try({
      wald_res <- Reduce('cbind',
                         object$wald_univariate())
      if(object$family$family == 'gaussian' &
         object$family$link == 'identity'){
        stat <- 't value'
      } else {
        stat <- 'z value'
      }
      colnames(wald_res) <- c('Estimate',
                              'Std. Error',
                              stat,
                              'CI LB',
                              'CI UB',
                              'Pr(>|z|)')
      wald_res <- wald_res[,c('Estimate',
                              'Std. Error',
                              stat,
                              'Pr(>|z|)',
                              'CI LB',
                              'CI UB')]
    }, silent = TRUE)
    if(any(class(tr) != 'try-error')){
      summary_list$coefficients <- tr
    } else {
      summary_list$coefficients <- cbind(unlist(object$B))
    }
  } else {
    summary_list$coefficients <- cbind(unlist(object$B))
  }
  summary_list$sigmasq_tilde <- object$sigmasq_tilde
  summary_list$effective_df <- object$N - object$trace_XUGX

  ## Custom print method for the summary
  print_summary <- function(x) {
    cat("Lagrangian Multiplier Smoothing Spline Model Summary\n")
    cat("====================================================\n")
    cat("Model Family (Link Function):",
        paste0(paste0(x$model_family)[1],
               " (",
               paste0(x$model_family)[2],
               ")", collapse = ""),
        "\n")
    cat("Observations:", x$observations, "\n")
    cat("Predictors:", x$predictors, "\n")
    cat("Partitions:", x$knots + 1, "\n")
    cat("Basis Functions per Partition:", x$basis_functions, "\n")
    cat("----------------------------------------------------\n")
    cat("Wald Inference: \n")
    print(x$coefficients)
    cat('\n')
    cat("Dispersion:", x$sigmasq_tilde, "\n")
    cat("Effective degrees of freedom:", x$effective_df, "\n")
    cat("Critical value for confidence intervals: ", x$cv, "\n")
    cat("----------------------------------------------------\n")
  }

  ## Set custom print method
  class(summary_list) <- "summary.lgspline"
  attr(summary_list, "print") <- print_summary
  return(summary_list)
}

#' Print method for lgspline summary objects
#'
#' @param x A summary.lgspline object
#' @export
print.summary.lgspline <- function(x, ...) {
  # Use the custom print function stored in the attribute
  print_method <- attr(x, "print")
  print_method(x)
  invisible(x)
}

#' Plot method for lgspline Objects
#'
#' @param object An lgspline model object
#' @param ... Additional plotting arguments
#' @noRd
plot.lgspline <- function(object, ...) {
  # Use the model's internal plotting function
  object$plot(...)
}

#' Predict method for lgspline Objects
#'
#' @param object An lgspline model object
#' @param ... Additional prediction arguments
#' @noRd
predict.lgspline <- function(object, ...) {
  # Use the model's internal predict function
  object$predict(...)
}

#' Extract model coefficients
#'
#' @param object An lgspline model object
#' @noRd
coef.lgspline <- function(object, ...) {
  object$B  # Return coefficients for each partition
}

#' Wald Univariate Inference for Lagrangian Multiplier Smoothing Spline Model
#'
#' Performs univariate Wald inference on model coefficients
#'
#' @param object An lgspline model object
#' @param ... Additional arguments passed to wald univariate
#' @return A matrix with columns for estimate, standard error, standardized test statistic, two-sided p-value,
#'         and confidence interval bounds for all coefficients
#' @noRd
wald_univariate <- function(object, ...) {
  ## Ensure variance-covariance matrix is available
  if(is.null(object$varcovmat)) {
    stop("Wald tests require return_varcovmat = TRUE during model fitting")
  }

  ## Compute Wald test results
  wald_tests <- object$wald_univariate(...)

  ## Construct results matrix
  res_matrix <- cbind(
    estimate = wald_tests$est,
    std_error = wald_tests$se,
    statistic = wald_tests$stat,
    p_value = wald_tests$pval,
    lower_ci = wald_tests$interval_lb,
    upper_ci = wald_tests$interval_ub
  )

  ## Row names based on model coefficients
  rownames(res_matrix) <- unlist(lapply(seq_along(object$B), function(k) {
    paste0("partition", k, "_", names(object$B[[k]]))
  }))

  return(res_matrix)
}

#' Posterior sampling method
#'
#' @param object An lgspline model object
#' @param ... Additional arguments for posterior generation
#' @noRd
generate_posterior <- function(object, ...) {
  object$generate_posterior(...)
}

#' Find function extrema
#'
#' @param object An lgspline model object
#' @param ... Additional arguments for extrema finding
#' @noRd
find_extremum <- function(object, ...) {
  object$find_extremum(...)
}

# Register S3 methods
#' @export
print.lgspline
#' @export
summary.lgspline
#' @export
print.summary.lgspline
#' @export
find_extremum
#' @export
generate_posterior
#' @export
plot.lgspline
#' @export
predict.lgspline
#' @export
coef.lgspline
#' @export
wald_univariate
