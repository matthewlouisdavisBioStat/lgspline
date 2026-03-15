#' Get Centers for Partitioning
#'
#' @param data Matrix of predictor data (already processed: binary/excluded
#'   columns zeroed out by the caller when \code{data_already_processed = TRUE})
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param cluster_args List with custom centers and kmeans args
#' @param cluster_on_indicators Include binary predictors in clustering
#' @param data_already_processed Logical; when TRUE the caller has already
#'   zeroed out binary / excluded columns, so \code{get_centers} skips its
#'   own binary-zeroing step.  Defaults to FALSE for backward compatibility.
#'
#' @details
#' Returns partition centers via:
#'
#' 1. Custom supplied centers if provided as a valid \eqn{(K+1) \times q} matrix
#'
#' 2. kmeans clustering on all non-spline variables if \code{cluster_on_indicators=TRUE}
#'    or if \code{data_already_processed=TRUE}
#'
#' 3. kmeans clustering excluding binary variables if \code{cluster_on_indicators=FALSE}
#'    and \code{data_already_processed=FALSE}
#'
#' @return Matrix of cluster centers
#'
#' @keywords internal
#' @export
get_centers <- function(data, K, cluster_args, cluster_on_indicators,
                        data_already_processed = FALSE) {

  ## If custom centers isn't null, return them
  #  These will be supplied to the first position
  #  Defaults to NA
  if(any(!is.na(cluster_args[[1]]))){
    return(cluster_args[[1]])

    ## Partition clusters including 0/1 predictors,
    #  or when the caller has pre-processed the data
    #  suppressWarnings prevents irrelevant "did not converge in 10 iterations"
  } else if(cluster_on_indicators || data_already_processed){
    km <- suppressWarnings(do.call(stats::kmeans,
                  c(list(x = data, centers = K + 1),
                    as.list(cluster_args[-1]))))

  } else {
    bin <- which(apply(data, 2, is_binary))
    data[,bin] <- 0
    km <- suppressWarnings(do.call(stats::kmeans,
                  c(list(x = data, centers = K + 1),
                    as.list(cluster_args[-1]))))
  }
  return(km$centers)
}

#' Find Neighboring Cluster Partitions Using Midpoint Distance Criterion
#'
#' @description
#' Identifies neighboring partitions by evaluating whether the midpoint between
#' cluster centers is closer to those centers than to any other center.
#'
#' @param centers Matrix; rows are cluster center coordinates
#' @param parallel Logical; use parallel processing
#' @param cl Cluster object for parallel execution
#' @param neighbor_tolerance Numeric; scaling factor for distance comparisons
#'
#' @return List where element i contains indices of centers neighboring center i
#'
#' @keywords internal
#' @export
find_neighbors <- function(centers, parallel, cl, neighbor_tolerance) {
  num_centers <- nrow(centers)
  neighbors <- vector("list", num_centers)

  if (!is.null(cl) & parallel) {
    ## Create index pairs for all comparisons
    pairs <- expand.grid(i = 1:(num_centers-1), j = 2:num_centers)
    pairs <- pairs[pairs$j > pairs$i, ]

    ## Split into chunks
    chunk_size <- max(1, ceiling(nrow(pairs) / length(cl)))
    chunks <- split(1:nrow(pairs), ceiling(seq_along(1:nrow(pairs))/chunk_size))

    ## Export data to cluster
    parallel::clusterExport(cl, "centers", envir = environment())

    ## Process chunks in parallel
    results <- parallel::parLapply(cl, chunks, function(chunk_indices) {
      chunk_results <- vector("list", length(chunk_indices))

      for (idx in seq_along(chunk_indices)) {
        pair_idx <- chunk_indices[idx]
        i <- pairs$i[pair_idx]
        j <- pairs$j[pair_idx]

        ## Check if midpoint is neighbor
        mid <- (centers[i,] + centers[j,])/2
        dist_ij <- sum((mid - centers[i,])^2) / neighbor_tolerance

        distances_to_others <- apply(centers[-c(i,j),,drop=FALSE], 1,
                                     function(k)sum((mid - k)^2))

        if(all(dist_ij < distances_to_others)) {
          chunk_results[[idx]] <- c(i, j)
        }
      }

      return(do.call(rbind, chunk_results[!sapply(chunk_results, is.null)]))
    })

    ## Combine results
    neighbor_pairs <- do.call(rbind, results)

    ## Convert to neighbor list
    if (!is.null(neighbor_pairs)) {
      for (k in 1:nrow(neighbor_pairs)) {
        i <- neighbor_pairs[k, 1]
        j <- neighbor_pairs[k, 2]
        neighbors[[i]] <- c(neighbors[[i]], j)
        neighbors[[j]] <- c(neighbors[[j]], i)
      }
    }

  } else {
    ## Sequential neighbor computation
    for (i in 1:(num_centers-1)) {
      for (j in (i+1):num_centers) {
        ## Get midpoint between centers i and j
        mid <- (centers[i,] + centers[j,])/2

        ## Compute scaled distance from midpoint to center i
        dist_ij <- sum((mid - centers[i,])^2) / neighbor_tolerance

        ## Get distances from midpoint to all other centers
        distances_to_others <- apply(centers[-c(i,j),,drop=FALSE], 1,
                                     function(k)sum((mid - k)^2))

        ## If midpoint closer to i,j than others, they are neighbors
        if(all(dist_ij < distances_to_others)) {
          neighbors[[i]] <- c(neighbors[[i]], j)
          neighbors[[j]] <- c(neighbors[[j]], i)
        }
      }
    }
  }

  return(neighbors)
}


#' Create Data Partitions Using Clustering
#'
#' @description
#' Partitions data support into clusters using Voronoi-like diagrams.
#' Standardization for clustering is handled internally; centers, knots, and
#' the returned \code{assign_partition} function are all on the **raw
#' (natural) scale** of the original predictors.
#'
#' @param data Numeric matrix of predictor variables (\strong{raw scale})
#' @param cluster_args Parameters for clustering
#' @param cluster_on_indicators Logical to include binary predictors
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#' @param do_not_cluster_on_these Columns to exclude from clustering
#' @param neighbor_tolerance Scaling factor for neighbor detection
#' @param standardize Logical; whether to standardize data internally before
#'   clustering.  Should equal \code{standardize_predictors_for_knots} from
#'   \code{lgspline.fit}.  Default \code{TRUE}.
#' @param standardize_mode Character; \code{"auto"} (default) selects
#'   \code{"minmax"} for a single effective clustering dimension and
#'   \code{"normal"} for multiple dimensions.  Can be forced to either value.
#' @param dummy_adder Small constant added to numerator during standardization
#'   to avoid exact-zero values (matches \code{lgspline.fit}'s
#'   \code{dummy_adder}).  Default \code{0}.
#' @param dummy_dividor Small constant added to denominator during
#'   standardization to avoid division by zero (matches
#'   \code{lgspline.fit}'s \code{dummy_dividor}).  Default \code{0}.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{centers}{Cluster center coordinates on the \strong{raw} scale.}
#'   \item{knots}{Knot points between centers on the \strong{raw} scale.}
#'   \item{assign_partition}{Function that accepts \strong{raw}-scale new data
#'     and returns integer-like partition assignments (0.5, 1.5, \ldots).}
#'   \item{neighbors}{List of neighboring partition indices.}
#'   \item{standardize_transf}{The forward standardization function used
#'     internally (for diagnostic use only).}
#'   \item{standardize_inv_transf}{The inverse standardization function
#'     (for diagnostic use only).}
#'   \item{centers_std}{Cluster centers on the standardized scale (for
#'     diagnostic use only).}
#' }
#'
#' @keywords internal
#' @export
make_partitions <- function(data,
                            cluster_args,
                            cluster_on_indicators,
                            K,
                            parallel,
                            cl,
                            do_not_cluster_on_these,
                            neighbor_tolerance,
                            standardize      = TRUE,
                            standardize_mode = "auto",
                            dummy_adder      = 0,
                            dummy_dividor    = 0) {

  q <- ncol(data)

  ## If K is 0, return as is
  if(K == 0){
    return(list(
      centers                = matrix()[-1,,drop=FALSE],
      knots                  = matrix()[-1,,drop=FALSE],
      assign_partition       = function(x) 1,
      neighbors              = list(),
      standardize_transf     = function(X) X,
      standardize_inv_transf = function(X) X,
      centers_std            = matrix()[-1,,drop=FALSE]
    ))
  }

  ## Internal standardization for clustering
  #  Cluster on a standardized copy of the data so that variables measured on
  #  very different scales do not dominate the kmeans objective.  Centers and
  #  knots are back-transformed to raw scale before returning so that
  #  downstream consumers (constraint matrix, return_list, etc.) see raw-scale
  #  coordinates without any additional transformation.
  ##
  #  standardize_mode "auto": choose "minmax" when only one effective clustering
  #  dimension remains after zeroing excluded / binary columns, otherwise
  #  "normal".  Can be forced to "minmax" or "normal" by the caller.

  if(standardize_mode == "auto"){
    ## Count columns that will actually vary during clustering
    effective_dims <- q - length(do_not_cluster_on_these)
    if(!cluster_on_indicators){
      bin_check <- which(apply(data, 2, is_binary))
      active    <- setdiff(seq_len(q), do_not_cluster_on_these)
      effective_dims <- length(setdiff(active, bin_check))
    }
    standardize_mode <- if(effective_dims <= 1L) "minmax" else "normal"
  }

  if(standardize){
    if(standardize_mode == "minmax"){
      std_mins <- apply(data, 2, min)
      std_maxs <- apply(data, 2, max)
      std_transf <- function(X){
        for(j in seq_len(ncol(X))){
          X[,j] <- (X[,j] - std_mins[j] + dummy_adder) /
            (std_maxs[j] - std_mins[j] + dummy_dividor)
        }
        X
      }
      std_inv_transf <- function(X){
        for(j in seq_len(ncol(X))){
          X[,j] <- X[,j] * (std_maxs[j] - std_mins[j] + dummy_dividor) +
            std_mins[j] - dummy_adder
        }
        X
      }
    } else {
      ## Normal (z-score) standardization
      std_means <- apply(data, 2, mean)
      std_sds   <- apply(data, 2, function(x) tryCatch(sd(x), error = function(e) 1))
      std_transf <- function(X){
        for(j in seq_len(ncol(X))){
          X[,j] <- (X[,j] - std_means[j] + dummy_adder) /
            (std_sds[j] + dummy_dividor)
        }
        X
      }
      std_inv_transf <- function(X){
        for(j in seq_len(ncol(X))){
          X[,j] <- X[,j] * (std_sds[j] + dummy_dividor) +
            std_means[j] - dummy_adder
        }
        X
      }
    }
    data_std <- std_transf(data)
  } else {
    data_std       <- data
    std_transf     <- function(X) X
    std_inv_transf <- function(X) X
  }

  ## Zero out excluded / binary columns in the standardized copy   -
  #  We operate on data_std from here so that find_neighbors and kmeans work
  #  on the same scale.  The originals (data, data_std) are not modified.
  data_cluster <- data_std

  if(!cluster_on_indicators){
    bin <- which(apply(data_cluster, 2, is_binary))
    data_cluster[,bin] <- 0
  }

  if(length(do_not_cluster_on_these) > 0){
    data_cluster[,do_not_cluster_on_these] <- 0
  }

  ## Compute CVT centers on the standardized scale       -
  #  Pass data_already_processed = TRUE so get_centers does not re-zero binary
  #  columns (we already did it above on data_cluster).
  #  Custom centers from cluster_args[[1]] are returned as-is; they are
  #  assumed to be on the STANDARDIZED scale when supplied via cluster_args
  #  (same assumption as before).  If the caller supplies raw-scale custom
  # centers they should pass standardize = FALSE.
  initial_points_std <- get_centers(data_cluster, K, cluster_args,
                                    cluster_on_indicators = TRUE,
                                    data_already_processed = TRUE)

  ## Re-order cluster indices by ascending L1-norm of their standardized
  #  center. This makes cluster 1 the "leftmost/smallest" and cluster K+1
  #  the "rightmost/largest", so that knot rownames like "1_2", "2_3" reflect
  #  actual spatial adjacency.  For 1-D this gives strict left-to-right order.
  #  For multi-D it gives a consistent, reproducible numbering that is
  #  insensitive to kmeans initialisation order.
  ## Re-index partitions so numbering follows the geometry more naturally
  #  and is reproducible across kmeans initialisations.
  l1_order <- order(apply(initial_points_std,
                          1,
                          function(r) sum(abs(r))))
  ## overwrite the initial_points_std above
  initial_points_std <- initial_points_std[l1_order, , drop = FALSE]

  ## Find neighboring partitions (on standardized scale)     -
  #  find_neighbors uses Euclidean distances; working on the standardized scale
  #  keeps comparisons numerically stable and scale-invariant.
  neighbors <- find_neighbors(initial_points_std,
                              parallel,
                              cl,
                              neighbor_tolerance)

  ## Compute knots as midpoints between neighboring centers
  find_knots <- function(pts, nbrs) {
    additional_points <- list()
    for (i in seq_along(nbrs)) {
      for (j in nbrs[[i]]) {
        if(j > i) next       ## avoid double-counting
        c1  <- pts[i,]
        c2  <- pts[j,]
        mid <- (c1 + c2) / 2
        additional_points[[length(additional_points) + 1]] <- mid
        names(additional_points)[length(additional_points)] <-
          paste0(i, '_', j)
      }
    }
    if(length(additional_points) > 0){
      return(do.call(rbind, additional_points))
    } else {
      return(matrix(nrow = 0, ncol = ncol(pts)))
    }
  }

  knots_std <- find_knots(initial_points_std, neighbors)

  ## Back-transform centers and knots to raw scale       -
  initial_points_raw <- std_inv_transf(initial_points_std)
  knots_raw <- if(nrow(knots_std) > 0) std_inv_transf(knots_std) else knots_std

  ## Partition-assignment closure
  #  The returned function accepts RAW-scale new data.  It standardizes
  #  internally using the same transform fitted above, then uses nearest-
  #  neighbor lookup against the standardized centers.
  #  Capture everything needed by value (inside the closure environment).
  assign_partition <- function(new_data) {

    ## Standardize new data with the same transform used for training
    new_data_std <- std_transf(new_data)

    ## Apply same column zeroing as during clustering
    if(!cluster_on_indicators){
      bin_new <- which(apply(new_data_std, 2, is_binary))
      if(length(bin_new) > 0) new_data_std[, bin_new] <- 0
    }
    if(length(do_not_cluster_on_these) > 0){
      new_data_std[, do_not_cluster_on_these] <- 0
    }

    ## Batch the nearest-neighbor lookup so large prediction sets do not
    #  allocate one full distance matrix at once.
    bytes_per_row <- 8 * nrow(initial_points_std)
    batch_size    <- floor((0.5 * 1024^3) / bytes_per_row)
    total_rows    <- nrow(new_data_std)
    assignments   <- numeric(total_rows)

    for(i in seq(1, total_rows, batch_size)) {
      nd_idx <- min(i + batch_size - 1, total_rows)
      batch  <- new_data_std[i:nd_idx, , drop = FALSE]
      nn     <- FNN::get.knnx(initial_points_std, batch, k = 1)
      assignments[i:nd_idx] <- nn$nn.index
    }

    return(assignments - 0.5)
  }

  ## Return raw-scale centers/knots together with the assignment closure
  #  and the internal standardization objects used to build them.
  rownames(initial_points_raw) <- paste0("center_", seq_len(nrow(initial_points_raw)))
  return(list(
    centers                = initial_points_raw,   ## RAW scale
    knots                  = knots_raw,            ## RAW scale
    assign_partition       = assign_partition,     ## Accepts RAW scale input
    neighbors              = neighbors,
    standardize_transf     = std_transf,           ## For diagnostics only
    standardize_inv_transf = std_inv_transf,       ## For diagnostics only
    centers_std            = initial_points_std    ## Standardized, for reference
  ))
}

