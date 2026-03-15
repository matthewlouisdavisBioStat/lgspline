#' Process Input Arguments for lgspline
#'
#' @description
#' Parses formula and data arguments, performs factor encoding,
#' resolves variable roles (spline, linear-with-interactions,
#' linear-without-interactions), constructs exclusion patterns,
#' and validates inputs. Called internally by \code{\link{lgspline}}.
#'
#' Users may call this function directly to inspect how their formula
#' and data are interpreted before fitting.
#'
#' @param predictors Default: NULL. Formula or numeric matrix/data frame
#'   of predictor variables.
#' @param y Default: NULL. Numeric response vector.
#' @param formula Default: NULL. Optional formula; alias for predictors
#'   when a formula object.
#' @param response Default: NULL. Alternative name for \code{y}.
#' @param data Default: NULL. Data frame for formula interface.
#' @param weights Default: NULL. Alias for \code{observation_weights}.
#' @param observation_weights Default: NULL. Numeric observation weight
#'   vector.
#' @param family Default: \code{gaussian()}. GLM family object.
#' @param K Default: NULL. Number of interior knots.
#' @param custom_knots Default: NULL. Custom knot matrix.
#' @param auto_encode_factors Default: TRUE. Logical; auto one-hot encode
#'   factor and character columns when using the formula interface.
#' @param include_2way_interactions Default: TRUE. Logical.
#' @param include_3way_interactions Default: TRUE. Logical.
#' @param just_linear_with_interactions Default: NULL. Integer vector or
#'   character vector of column names.
#' @param just_linear_without_interactions Default: NULL. Integer vector or
#'   character vector of column names.
#' @param exclude_interactions_for Default: NULL. Integer vector or character
#'   vector of column names.
#' @param exclude_these_expansions Default: NULL. Character vector of
#'   expansion names to exclude.
#' @param offset Default: \code{c()}. Vector of column indices or names
#'   to include as offsets.
#' @param no_intercept Default: FALSE. Logical; remove intercept.
#' @param do_not_cluster_on_these Default: \code{c()}. Vector of column
#'   indices or names to exclude from clustering.
#' @param include_quartic_terms Default: NULL. Logical or NULL.
#' @param cluster_args Default: \code{c(custom_centers = NA, nstart = 10)}.
#'   Named vector of clustering arguments.
#' @param include_warnings Default: TRUE. Logical.
#' @param dummy_fit Default: FALSE. Logical; run the full preprocessing path but stop short of fitting nonzero coefficients.
#' @param include_constrain_second_deriv Default: TRUE. Logical.
#' @param standardize_response Default: TRUE. Logical.
#' @param ... Additional arguments (checked for depreciated names).
#'
#' @return A named list containing:
#' \describe{
#'   \item{predictors}{Numeric matrix of predictor variables with column
#'     names stripped for positional indexing.}
#'   \item{y}{Numeric response vector.}
#'   \item{og_cols}{Character vector of original predictor column names,
#'     or NULL if none were available.}
#'   \item{replace_colnames}{Logical; TRUE if og_cols is available and
#'     column renaming should be applied post-fit.}
#'   \item{just_linear_with_interactions}{Integer vector of column indices,
#'     or NULL.}
#'   \item{just_linear_without_interactions}{Integer vector of column indices,
#'     or NULL.}
#'   \item{exclude_interactions_for}{Integer vector of column indices,
#'     or NULL.}
#'   \item{exclude_these_expansions}{Character vector of positional-notation
#'     expansion names, or NULL.}
#'   \item{offset}{Integer vector of column indices, or \code{c()}.}
#'   \item{no_intercept}{Logical.}
#'   \item{do_not_cluster_on_these}{Numeric vector of column indices,
#'     or \code{c()}.}
#'   \item{observation_weights}{Numeric vector or NULL.}
#'   \item{K}{Integer or NULL, possibly updated by cluster_args or
#'     all-linear detection.}
#'   \item{include_3way_interactions}{Logical, possibly updated by
#'     formula parsing.}
#'   \item{include_quartic_terms}{Logical or NULL, possibly updated
#'     by number of predictors.}
#'   \item{data}{Data frame, possibly with factor columns one-hot encoded.}
#'   \item{include_constrain_second_deriv}{Logical, possibly set FALSE
#'     when no numeric predictors remain.}
#'   \item{factor_groups}{Named list mapping original factor column names to
#'     integer vectors of their one-hot indicator column positions within
#'     the predictor matrix. Used by lgspline.fit to impose sum-to-zero
#'     constraints on encoded factor levels. NULL when no factors were
#'     encoded, or when the formula interface was not used.}
#' }
#'
#' @seealso \code{\link{lgspline}} for the main fitting interface.
#'
#' @examples
#' \dontrun{
#' data("Theoph")
#' df <- Theoph[, c("Time", "Dose", "conc", "Subject")]
#' processed <- process_input(
#'   predictors = conc ~ spl(Time) + Time*Dose,
#'   data = df,
#'   auto_encode_factors = TRUE,
#'   include_warnings = TRUE
#' )
#' str(processed$predictors)
#' processed$og_cols
#' processed$just_linear_without_interactions
#' processed$factor_groups
#' }
#'
#' @keywords internal
#' @export
process_input <- function(
    predictors = NULL,
    y = NULL,
    formula = NULL,
    response = NULL,
    data = NULL,
    weights = NULL,
    observation_weights = NULL,
    family = gaussian(),
    K = NULL,
    custom_knots = NULL,
    auto_encode_factors = TRUE,
    include_2way_interactions = TRUE,
    include_3way_interactions = TRUE,
    just_linear_with_interactions = NULL,
    just_linear_without_interactions = NULL,
    exclude_interactions_for = NULL,
    exclude_these_expansions = NULL,
    offset = c(),
    no_intercept = FALSE,
    do_not_cluster_on_these = c(),
    include_quartic_terms = NULL,
    cluster_args = c(custom_centers = NA, nstart = 10),
    include_warnings = TRUE,
    dummy_fit = FALSE,
    include_constrain_second_deriv = TRUE,
    standardize_response = TRUE,
    ...
){

  ## Backward compat: expansions_only in ... maps
  #  to dummy_fit with deprecation warning
  #  shur -> schur for correct spelling
  dots <- list(...)
  if(!is.null(dots$expansions_only) && dots$expansions_only){
    if(include_warnings){
      warning(
        "'expansions_only' is depreciated; use 'dummy_fit = TRUE' instead.",
        " The full pipeline now executes with coefficients set to 0.")
    }
    dummy_fit <- TRUE
  }
  if(!is.null(dots$shur_correction_function)){
    if(include_warnings){
      warning("'shur_correction_function' option is depreciated, and won't be",
              "applied; repeat with 'schur_correction_function' instead.")
    }
  }

  ## Track one-hot factor groups for later sum-to-zero constraints
  #  Stores original factor name -> indicator column names until the
  #  final predictor matrix has been resolved.
  factor_groups <- list()

  ## Update naming conventions, if first argument is a formula and second is a
  # data frame, assumed by R-like interfaces for user convenience.
  if(any(!is.null(predictors)) & any(!is.null(y))){
    if(any(inherits(y,'data.frame') & inherits(predictors, "formula"))){
      data <- y
    }
  }

  ## Update naming conventions, if response supplied in place of "y" for user
  # convenience.
  if(any(is.null(predictors)) & any(!is.null(formula))){
    predictors <- formula
  } else if(any(is.null(predictors))){
    stop(
      '\n \t Predictors argument is NULL without formula supplied.',
      ' Either supply a formula to predictors OR formula argument, or a',
      ' data frame to predictors argument, or a matrix of numeric predictors',
      ' to the predictor argument. \n')
  }

  ## Update naming conventions, if response supplied in place of "y" for user
  # convenience.
  if(any(is.null(y)) & any(!is.null(response))){
    y <- response
    response <- NULL
  }

  ## Weights is just an R-friendly argument to be passed to observation_weights
  # actually used by the function.
  if(any(is.null(observation_weights)) &
     any(!is.null(weights))){
    observation_weights <- weights
    weights <- NULL
  }

  ## Check cluster args for compatibility
  if(any(!is.na(cluster_args[[1]]))){
    ncluster <- try({nrow(cluster_args[[1]])},
                    silent = TRUE)
    if(any(inherits(ncluster, 'try-error'))){
      stop('\n \t custom_centers should be a matrix; do not include any other',
           ' arguments within cluster_args if you include custom_centers. \n')
    }
    if(!is.null(K)){
      if(ncluster != (K+1) & include_warnings){
        warning('\n \t K must be equal to number of custom_centers minus 1. ',
                'Updating K for compatibility. \n')
        K <- ncluster - 1
      }
    } else {
      K <- ncluster - 1
    }
  }

  ## Check data and formula argument
  if(!is.null(data) & !inherits(predictors, "formula")) {
    stop("\n \t If submitting data argument, formula must be supplied and',
         ' variables must match. Otherwise, use predictors and y (or response)",
         " arguments directly.\n",
         "\t Example: lgspline(y ~ spl(x1, x2) + x3 + x4*x5, data = my_data)\n",
         "\t Example: lgspline(y ~ ., data = my_data)\n")
  }

  ## Handle formula interface
  if(inherits(predictors, "formula")) {

    ## Allow s() as an alias for spl() for mgcv-style formulas.
    #  Replace s(...) with spl(...) in formula before any parsing.
    #  Uses pattern matching to avoid replacing 's' inside other function
    #  names like abs(), cos(), is(), etc. The regex anchors on a
    #  non-alphanumeric/non-underscore/non-dot character (or start of string)
    #  immediately before the 's(', so only bare 's(' is matched.
    form_str_alias <- paste(deparse(predictors, width.cutoff = 500L),
                            collapse = "")
    form_str_alias <- gsub("(^|[^a-zA-Z0-9_.])s\\(", "\\1spl(", form_str_alias)
    predictors <- as.formula(form_str_alias)

    ## Check data argument
    if(is.null(data)) {
      stop(
        "\n \t When using formula interface, data argument must be ",
        "provided.\n",
        "\t Example: lgspline(y ~ spl(x1, x2) + x3 + x4*x5, data = my_data)\n",
        "\t Example: lgspline(y ~ ., data = my_data)\n")
    }

    ## Try to coerce data to data.frame
    tryCatch({
      data <- as.data.frame(data)
    }, error = function(e) {
      stop("\n \t Could not coerce data argument to data.frame. ",
           "Please provide data in a format coercible to data.frame. ",
           "Examples: data.frame, tibble, or matrix.")
    })

    ## Auto one-hot encode factor/character columns
    # For reviewers: the decision to include one-hot rather than
    # dummy-intercept is for interpretation as conditional random effects
    # models. We would not like this to drop a random effect level,
    # The flat_ridge_penalty should address lack of identifiability to an extent
    # Users can manually supply the variable names after one-hot encoding
    # if dummy-intercept is truly desired.
    if(auto_encode_factors && !is.null(data)){
      non_numeric_cols <- which(!sapply(data, is.numeric))
      ## Exclude response column from encoding
      resp_name <- NULL
      if(inherits(predictors, "formula")){
        resp_name <- as.character(predictors[[2]])
      }
      non_numeric_cols <- non_numeric_cols[
        !names(non_numeric_cols) %in% resp_name
      ]
      if(length(non_numeric_cols) > 0){
        for(col_idx in rev(non_numeric_cols)){
          col_name <- colnames(data)[col_idx]
          ## One-hot encode this column
          encoded <- create_onehot(data[[col_name]])
          colnames(encoded) <- paste0(col_name, "_", colnames(encoded))
          ## Remove original, append indicators
          data <- cbind(data[, -col_idx, drop = FALSE], encoded)
          ## Track these as linear-without-interactions only if the
          #  original column was not involved in any interaction term
          new_names <- colnames(encoded)
          ## Track these indicator columns for the later sum-to-zero constraints.
          #  Stored by original factor name, then resolved to indices below.
          factor_groups[[col_name]] <- new_names
          form_check <- paste(deparse(predictors, width.cutoff = 500L),
                              collapse = "")
          col_in_interaction <- grepl(
            paste0("(\\*|:)\\s*", col_name, "(\\s|\\)|\\+|$)"), form_check
          ) || grepl(
            paste0("(\\+|\\(|~|\\s)", col_name, "\\s*(\\*|:)"), form_check
          )
          if(!col_in_interaction){
            if(is.null(just_linear_without_interactions)){
              just_linear_without_interactions <- new_names
            } else {
              just_linear_without_interactions <- c(
                just_linear_without_interactions, new_names
              )
            }
          }
          ## Update formula if needed: replace col_name with encoded names
          if(inherits(predictors, "formula")){
            form_str <- deparse(predictors)
            replacement <- paste0("(", paste(new_names, collapse = " + "), ")")
            form_str <- gsub(col_name, replacement, form_str, fixed = TRUE)
            predictors <- as.formula(form_str)
          }
        }
      }
    }

    ## Check column names are present
    if(any(is.null(colnames(data)))){
      stop('\n \t Column names of data must be supplied if formula is',
           ' supplied, and column names must match what is provided in formula',
           '.\n')
    }

    ## Stringed formula
    form_paste0 <- paste0(predictors)

    ## If formula is y ~ ., replace with y ~ spl(x1, x2, ...)
    if(form_paste0[1] == '~' &
       (gsub(' ', '', form_paste0[3]) %in% c('.', '0+.','.+0')) &
       length(form_paste0) == 3){
      if(form_paste0[2] == '.'){
        stop('\n \t . ~ . is not valid for this function, specify the',
             ' response directly, like "y ~ . \n"')
      } else if(gsub(' ', '', form_paste0[3]) == '0+.' |
                gsub(' ', '', form_paste0[3]) == '.+0'){
        ## No intercept
        predictors <- as.formula(
          paste0(form_paste0[2],
                 ' ~ 0 + spl(',
                 paste(colnames(data)[-which(colnames(data) ==
                                               form_paste0[2])],
                       collapse = ', '),
                 ')')
        )
      } else {
        ## With intercept
        predictors <- as.formula(
          paste0(form_paste0[2],
                 ' ~ spl(',
                 paste(colnames(data)[-which(colnames(data) ==
                                               form_paste0[2])],
                       collapse = ', '),
                 ')')
        )
      }
    }

    ## [Change PATCH] Handle spl(...) + . formula pattern.
    #  When the formula contains spl() terms AND a trailing dot,
    #  replace the dot with explicit column names of non-spl,
    #  non-response columns and treat them as
    #  just_linear_without_interactions.
    #  For example:
    #    conc ~ spl(Time) + Time*Dose + .
    #  becomes (with data having columns Time, Dose, Subject):
    #    conc ~ spl(Time) + Time*Dose + Subject
    #  and Subject is added to just_linear_without_interactions.
    #  Variables already named elsewhere in the formula (inside spl(),
    #  in interactions, etc.) are excluded from the dot expansion so
    #  they do not appear twice. The downstream formula parser then
    #  correctly classifies variables that appear in interactions as
    #  linear_with_int rather than linear_without_int.
    form_str_full <- paste(deparse(predictors, width.cutoff = 500L),
                           collapse = "")
    has_spl <- grepl("spl\\(", form_str_full)
    has_dot <- grepl("\\+\\s*\\.", form_str_full) |
      grepl("\\.\\s*\\+", form_str_full)
    if(has_spl & has_dot){
      ## Extract spl variable names
      spl_matches <- gregexpr("spl\\(([^)]+)\\)", form_str_full)
      spl_contents <- regmatches(form_str_full, spl_matches)[[1]]
      spl_var_names <- unique(trimws(unlist(strsplit(
        gsub("spl\\((.*)\\)", "\\1", spl_contents), ","
      ))))

      ## Response name
      resp_nm <- as.character(predictors[[2]])

      ## All variables already explicitly named in the formula
      #  (outside spl() and the dot itself). Parse the terms to find them.
      form_terms_tmp <- try(terms(predictors), silent = TRUE)
      if(!inherits(form_terms_tmp, 'try-error')){
        explicit_vars <- unique(c(
          spl_var_names,
          rownames(attr(form_terms_tmp, "factors"))[-1]
        ))
      } else {
        explicit_vars <- spl_var_names
      }

      ## Remaining columns = all data columns minus response minus explicit
      remaining_cols <- setdiff(colnames(data),
                                c(resp_nm, explicit_vars))

      if(length(remaining_cols) > 0){
        ## Replace the dot with explicit column names
        dot_replacement <- paste(remaining_cols, collapse = " + ")
        ## Handle both "+ ." and ". +" patterns
        form_str_new <- sub("\\+\\s*\\.", paste0("+ ", dot_replacement),
                            form_str_full)
        form_str_new <- sub("\\.\\s*\\+", paste0(dot_replacement, " + "),
                            form_str_new)
        predictors <- as.formula(form_str_new)

        ## Append to just_linear_without_interactions
        if(is.null(just_linear_without_interactions)){
          just_linear_without_interactions <- remaining_cols
        } else {
          just_linear_without_interactions <- unique(c(
            just_linear_without_interactions, remaining_cols
          ))
        }
      }
    }

    ## If offset() appears in the formula, remove offset from the formula
    # and set offset = names of terms inside the offset() operator
    if(any(grepl('offset', paste0(predictors)))){
      if(any(grepl("offset\\(", deparse(predictors)))){
        ## Extract the original formula as string
        form_str <- deparse(predictors)

        ## Find all offset terms using regex
        offset_matches <- gregexpr("offset\\(([^)]+)\\)", form_str)
        offset_vars <- regmatches(form_str, offset_matches)[[1]]

        ## Extract variable names from offset terms
        offset_var_names <- gsub("offset\\((.*)\\)", "\\1", offset_vars)

        ## Remove whitespace
        offset_var_names <- trimws(offset_var_names)

        ## Add to offset vector
        offset <- c(offset, offset_var_names)

        ## Replace offset(var) with var in formula
        for(i in seq_along(offset_vars)) {
          form_str <- gsub(offset_vars[i], offset_var_names[i],
                           form_str,
                           fixed=TRUE)
        }

        ## Convert back to formula
        predictors <- as.formula(form_str)
      }
    }

    ## Check that all formula predictors are numeric in data
    terms <- terms(predictors)
    term_labels <- attr(terms, "term.labels")
    non_numeric <- c(1:ncol(data))[!sapply(data, is.numeric)]
    if(length(non_numeric) > 0) {
      for(v in non_numeric){
        if(any(term_labels == colnames(data)[v])){
          if(auto_encode_factors){
            stop("Auto-encoding failed for column '",
                 colnames(data)[v], "'. Check that it is a factor ",
                 "or character vector.")
          } else {
            stop("Non-numeric columns detected in predictors. ",
                 "Set auto_encode_factors = TRUE to automatically ",
                 "one-hot encode factor/character columns, or ",
                 "convert these to numeric before proceeding. ",
                 "You can use create_onehot() on categorical ",
                 "columns of your dataset to obtain binary ",
                 "indicator variables if you need to. Remember ",
                 "to remove the original categorical variable ",
                 "and append the indicators to your data ",
                 "before adjusting your formula and call ",
                 "to lgspline(...).\n\nSee ?create_onehot for an example.")
          }
        }
      }
    }

    ## Parse rest of formula
    has_3way <- any(attr(terms, "order") == 3)
    if(include_3way_interactions) include_3way_interactions <- has_3way

    ## Check for no intercept specification in formula
    if(inherits(predictors, "formula")) {
      formula_text <- gsub(" ", "", Reduce(paste, deparse(predictors)))
      if(grepl("\\+0|0\\+", formula_text)) {
        no_intercept <- TRUE
      }
    }

    ## Get the factors matrix which shows interaction structure
    factors <- attr(terms, "factors")

    ## Initialize term containers
    spline_terms <- character()
    linear_no_int <- character()
    linear_with_int <- character()

    ## First pass to identify spline terms and extract their variables
    for(term in term_labels) {
      if(grepl("^spl\\(.*\\)$", term)) {
        vars <- gsub("^spl\\((.*)\\)$", "\\1", term)
        spline_terms <- c(spline_terms, trimws(strsplit(vars, ",")[[1]]))
      }
    }

    ## Extract explicit interactions from formula
    formula_interactions <- term_labels[attr(terms, "order") > 1]
    allowed_interactions <- sapply(formula_interactions, function(term) {
      if(grepl(":", term)) {
        vars <- strsplit(term, ":")[[1]]
        # Only keep interactions where NO term is a spline term
        if(!any(vars %in% spline_terms)) {
          return(term)
        }
      }
      return(NULL)
    })

    ## Identify interaction structure using factors matrix
    var_names <- rownames(factors)[-1] # Remove responses

    ## Extract all types of variables in formula
    all_formula_vars <- unique(c(spline_terms, var_names))
    formula_cols <- which(colnames(data) %in% all_formula_vars)

    ## Match custom exclusions nominal to numeric formula columns, if available
    if(inherits(exclude_interactions_for, 'character')){
      exclude_interactions_for <-
        unlist(lapply(exclude_interactions_for,
                      function(var){
                        grep(var,
                             colnames(data)[formula_cols])
                      }))

    }
    if(inherits(just_linear_with_interactions,  'character')){
      just_linear_with_interactions <-
        unlist(lapply(just_linear_with_interactions,
                      function(var){
                        grep(var,
                             colnames(data)[formula_cols])
                      }))

    }
    if(inherits(just_linear_without_interactions, 'character')){
      just_linear_without_interactions <-
        unlist(lapply(just_linear_without_interactions,
                      function(var){
                        grep(var,
                             colnames(data)[formula_cols])
                      }))

    }
    if(length(exclude_these_expansions) > 0){
      for(ii in 1:length(exclude_these_expansions)){
        exps <- exclude_these_expansions[ii]
        if(substr(exps, 1, 1) != "_" |
           substr(exps, length(exps), length(exps)) != "_"){
          for(jj in 1:length(formula_cols)){
            if(grepl(colnames(data)[formula_cols[jj]],
                     exclude_these_expansions[ii])){
              exclude_these_expansions[ii] <- gsub(
                colnames(data)[formula_cols[jj]],
                paste0('_', jj, '_'),
                exclude_these_expansions[ii]
              )
            }
          }
        }
      }
    }
    if(inherits(offset,  'character')){
      offset <-
        unlist(lapply(offset,
                      function(var){
                        grep(var,
                             colnames(data)[formula_cols])
                      }))
    }

    ## Resolve character do_not_cluster_on_these to numeric
    #  column indices within the formula path. Without this,
    #  make_partitions' assign_partition closure fails when indexing
    #  matrices after transf() strips colnames.
    if(length(do_not_cluster_on_these) > 0 &&
       is.character(do_not_cluster_on_these)){
      do_not_cluster_on_these <-
        unlist(lapply(do_not_cluster_on_these,
                      function(var){
                        idx <- which(colnames(data)[formula_cols] == var)
                        if(length(idx) == 0){
                          idx <- grep(var, colnames(data)[formula_cols])
                        }
                        idx
                      }))
      do_not_cluster_on_these <- unique(do_not_cluster_on_these)
    }

    ## Convert factor_groups from column names to predictor-matrix indices.
    #  lgspline.fit uses these positions to build the sum-to-zero constraints.
    #  Groups that cannot be resolved cleanly are dropped silently.
    if(length(factor_groups) > 0){
      factor_groups <- lapply(factor_groups, function(col_names){
        unlist(lapply(col_names, function(nm){
          idx <- which(colnames(data)[formula_cols] == nm)
          if(length(idx) == 0) idx <- grep(paste0("^", nm, "$"),
                                           colnames(data)[formula_cols])
          idx
        }))
      })
      ## Remove any groups that could not be resolved or have only 1 level
      factor_groups <- factor_groups[sapply(factor_groups, length) > 1]
    }

    ## Non-spline variables that appear in interactions (from factors matrix)
    interaction_terms <- term_labels[attr(terms, "order") > 1]
    nonspline_interact_vars <- unique(sapply(interaction_terms,
                                             function(term) {
                                               if(grepl(":", term)) {
                                                 vars <-
                                                   strsplit(term,":")[[1]]
                                                 vars[!vars %in%
                                                        spline_terms]
                                               }
                                             }))

    ## Safe initialization of allowed_interaction_pairs
    if(!exists("allowed_interaction_pairs")){
      allowed_interaction_pairs <- list()
    }

    ## Add explicit interactions for spline terms
    if(length(spline_terms) > 1) {
      ## Generate all possible interactions between spline terms
      spline_interactions <- utils::combn(spline_terms, 2, simplify=FALSE)
      spline_triplets <- if(length(spline_terms) >= 3 &
                            include_3way_interactions) {
        utils::combn(spline_terms, 3, simplify=FALSE)
      } else {
        list()
      }

      ## Add 2-way interactions
      for(pair in spline_interactions) {
        interaction_terms <- c(interaction_terms, paste(pair, collapse=":"))
      }

      ## Add 3-way interactions
      for(triplet in spline_triplets) {
        interaction_terms <- c(interaction_terms, paste(triplet, collapse=":"))
      }

      ## Also register spline triplets as allowed interaction pairs
      ## so they are not excluded by the grouping logic below
      for(triplet in spline_triplets) {
        allowed_interaction_pairs[[length(allowed_interaction_pairs) + 1]] <-
          list(pair = triplet,
               transforms = c(
                 paste(triplet, collapse = "x"),
                 paste(rev(triplet), collapse = "x")
               )
          )
      }
    }

    ## Get indices of variables in raw expansions (after response removed)
    resp_ind <- which(colnames(data) == paste0(terms[[2]]))
    var_positions <- match(var_names, colnames(data[,formula_cols]))

    ## Get allowed interaction pairs
    allowed_pairs <- lapply(interaction_terms[grepl(":", interaction_terms)],
                            function(term) {
                              vars <- strsplit(term, ":")[[1]]
                              match(vars, colnames(data[,formula_cols]))
                            })

    ## Generate exclusion patterns for non-spline vars
    exclude_patterns <- c()

    ## Get all possible 2-way interactions between ANY variables
    vars <- colnames(data[,formula_cols])
    if(include_2way_interactions){
      for(ii in seq_along(vars)) {
        for(jj in seq_along(vars)) {
          if(ii != jj) {
            pattern <- get_interaction_patterns(c(vars[ii], vars[jj]))

            ## Skip if not in formula
            if(any(!(c(vars[ii], vars[jj]) %in% all_formula_vars))){
              next
            }

            ## Skip if any variable in exclude_interactions_for
            if(!is.null(exclude_interactions_for)){
              if(any(c(ii, jj) %in% exclude_interactions_for)){
                next
              }
            }

            ## Only keep spline-spline interactions and explicit interactions
            if(!all(c(vars[ii], vars[jj]) %in% spline_terms) &
               length(spline_terms) > 0 &
               all(c(ii, jj) %in% formula_cols)){
              exclude_patterns <- c(exclude_patterns, pattern)
            }
          }
        }
      }
    }

    ## Add all possible 3-way interactions to exclusions
    if(include_3way_interactions) {
      for(ii in seq_along(vars)) {
        for(jj in seq_along(vars)) {
          for(kk in seq_along(vars)) {
            if(ii != jj && jj != kk && ii != kk) {
              triplet_vars <- c(vars[ii], vars[jj], vars[kk])
              pattern <- get_interaction_patterns(triplet_vars)

              ## Skip if not in formula
              if(any(!(triplet_vars %in% all_formula_vars))){
                next
              }

              ## Skip if any variable in exclude_interactions_for
              if(!is.null(exclude_interactions_for)){
                if(any(c(ii, jj, kk) %in% exclude_interactions_for)){
                  next
                }
              }

              ## Skip if only spline terms
              if(all(triplet_vars %in% spline_terms)){
                next
              }

              ## If ANY var is a spline term but not ALL are spline terms,
              # we should exclude this interaction
              if(any(triplet_vars %in% spline_terms) &&
                 !all(triplet_vars %in% spline_terms)) {
                exclude_patterns <- c(exclude_patterns, pattern)
                next
              }

              ## For non-spline terms, allow explicitly specified interactions
              if(!any(triplet_vars %in% spline_terms)) {
                ## Get interaction terms that could involve these variables
                relevant_terms <- interaction_terms[grepl(paste(triplet_vars,
                                                                collapse="|"),
                                                          interaction_terms)]

                ## Check if this exact triplet exists in any order
                matches_interaction <- any(sapply(strsplit(relevant_terms, ":"),
                                                  function(term) {
                                                    length(term) == 3 &&
                                                      all(sort(triplet_vars) ==
                                                            sort(term))
                                                  }))

                if(!matches_interaction) {
                  exclude_patterns <- c(exclude_patterns, pattern)
                }
              }
            }
          }
        }
      }
    }

    ## Helper-vector for the next step, in isolating which interactions are
    #  not part of the same grouping/block.
    #  The above only works if we have only 1 additive interaction effect.
    #  the code below corrects for when we have multiple
    #  (for example, a*b + c*d)
    #  and want to prevent unspecified interactions
    #  (for example, remove a*d and b*c)
    different_grouping_exclusions <- c()

    ## Explicitly track allowed interactions from * and : terms
    if(!exists("allowed_interaction_pairs")){
      allowed_interaction_pairs <- list()
    }

    ## Track which interactions are explicitly specified
    for(term in term_labels) {
      # Handle * interactions
      if(grepl("\\*", term)) {
        vars <- trimws(strsplit(term, "\\*")[[1]])

        # Ensure exactly 2 variables in the interaction
        if(length(vars) == 2) {
          allowed_interaction_pairs[[length(allowed_interaction_pairs) + 1]] <-
            list(pair = vars,
                 transforms = c(
                   paste(vars, collapse = "x"),
                   paste(rev(vars), collapse = "x")
                 )
            )
        }
      }

      # Handle : interactions
      if(grepl(":", term)) {
        vars <- trimws(strsplit(term, ":")[[1]])

        # Ensure exactly 2 variables in the interaction
        if(length(vars) == 2) {
          allowed_interaction_pairs[[length(allowed_interaction_pairs) + 1]] <-
            list(pair = vars,
                 transforms = c(
                   paste(vars, collapse = "x"),
                   paste(rev(vars), collapse = "x")
                 )
            )
        }
      }
    }

    ## Identify spline block variables
    spline_block_vars <- c()
    for(term in term_labels) {
      if(grepl("^spl\\(.*\\)$", term)) {
        vars <- gsub("^spl\\((.*)\\)$", "\\1", term)
        term_vars <- trimws(strsplit(vars, ",")[[1]])
        spline_block_vars <- c(spline_block_vars, term_vars)
      }
    }

    ## Generate exclusions for any interaction not in the allowed pairs
    vars <- colnames(data[,formula_cols])
    if(include_2way_interactions | include_3way_interactions){
      for(ii in seq_along(vars)) {
        for(jj in seq_along(vars)) {
          if(ii != jj) {

            pair <- c(vars[ii], vars[jj])
            current_transform <- paste(pair, collapse = "x")

            ## Skip if not in formula
            if(any(!(pair %in% all_formula_vars))){
              next
            }

            ## Skip if any variable in exclude_interactions_for
            if(!is.null(exclude_interactions_for)){
              if(any(c(ii, jj) %in% exclude_interactions_for)){
                next
              }
            }

            ## Check if this is an allowed interaction
            is_allowed_interaction <- FALSE
            for(allowed in allowed_interaction_pairs) {
              if(setequal(pair, allowed$pair) ||
                 any(current_transform == allowed$transforms)) {
                is_allowed_interaction <- TRUE
                break
              }
            }

            ## Check if this is a spline block interaction (2-way or 3-way)
            is_spline_block_interaction <- all(pair %in% spline_block_vars)

            ## If not an allowed or spline block interaction,
            # generate exclusion patterns
            if(!is_allowed_interaction && !is_spline_block_interaction) {
              patterns <- get_interaction_patterns(pair)
              different_grouping_exclusions <- c(different_grouping_exclusions,
                                                 patterns)
            }
          }
        }
      }
    }

    ## Remove duplicates
    different_grouping_exclusions <- unique(different_grouping_exclusions)

    ## Remove explicitly allowed interactions from exclusions
    for(term in interaction_terms) {
      vars <- strsplit(term, ":")[[1]]
      allowed <- get_interaction_patterns(vars)
      exclude_patterns <- setdiff(exclude_patterns, allowed)
    }
    exclude_patterns <- unique(c(exclude_patterns,
                                 different_grouping_exclusions))

    ## Convert to positional notation
    vars <- colnames(data)[formula_cols]
    for(ii in seq_along(formula_cols)) {
      exclude_patterns <- gsub(vars[ii], paste0("_", ii, "_"), exclude_patterns)
    }

    ## Append to custom exclusions
    if(!is.null(exclude_these_expansions)){
      exclude_these_expansions <- c(exclude_these_expansions,
                                    exclude_patterns)
    } else if (length(exclude_patterns) > 0){
      exclude_these_expansions <- exclude_patterns
    }

    ## For each variable, determine if linear with or without interactions
    for(var in var_names) {
      if(length(spline_terms) > 0){
        if(var %in% spline_terms) next # Skip spline terms
      }

      ## Check if this variable appears in any interactions
      var_terms <- which(factors[var,] > 0)
      if(length(var_terms) > 0) {
        ## If any term containing this variable has order > 1,
        # its in an interaction
        if(any(attr(terms, "order")[var_terms] > 1)) {
          linear_with_int <- c(linear_with_int, var)
        } else {
          linear_no_int <- c(linear_no_int, var)
        }
      } else {
        linear_no_int <- c(linear_no_int, var)
      }
    }

    ## Remove duplicates and ensure proper separation
    linear_with_int <- unique(linear_with_int)
    linear_no_int <- setdiff(unique(linear_no_int),
                             c(spline_terms, linear_with_int))
    linear_no_int <- linear_no_int[!(substr(linear_no_int,
                                            1,
                                            4) == 'spl(')]

    ## Create predictors matrix and response for compatibility
    predictors <- data[, formula_cols, drop = FALSE]
    y <- data[, resp_ind]

    ## Convert variable names to column indices
    new_just_linear_without_interactions <- match(linear_no_int,
                                                  colnames(predictors))
    new_just_linear_with_interactions <- match(linear_with_int,
                                               colnames(predictors))
    new_just_linear_without_interactions <- new_just_linear_without_interactions[
      !(new_just_linear_without_interactions %in% spline_terms)
    ]
    new_just_linear_with_interactions <- new_just_linear_with_interactions[
      !(new_just_linear_with_interactions %in% spline_terms)
    ]
    if(is.null(just_linear_with_interactions)){
      just_linear_with_interactions <- new_just_linear_with_interactions
    } else{
      just_linear_with_interactions <- unique(c(
        just_linear_with_interactions,
        new_just_linear_with_interactions
      ))
    }
    if(is.null(just_linear_without_interactions)){
      just_linear_without_interactions <- new_just_linear_without_interactions
    } else {
      just_linear_without_interactions <- unique(c(
        just_linear_without_interactions,
        new_just_linear_without_interactions
      ))
    }
  }

  ## Not a formula - try to coerce to matrix
  tryCatch({
    predictors <- as.matrix(predictors)
  }, error = function(e) {
    stop("\n \t Could not coerce predictors to matrix. ",
         "predictors must be either a formula or an object coercible to matrix.",
         "Examples:\n",
         "  Formula: lgspline(y ~ spl(x1, x2) + x3, data = my_data)\n",
         "  Matrix:  lgspline(predictors = Tmat, y = y) \n")
  })

  ## Check numeric type
  if(any(!is.numeric(predictors))){
    stop("\n \t predictors matrix must be numeric. ",
         "Please convert categorical variables to numeric indicators. \n")
  }

  ## Check response for missings
  if(any(is.na(y) | is.nan(y) | !is.finite(y))){
    stop("\n \t NA, NaN, or infinite value detected in response. \n")
  }

  ## Original predictor names
  og_cols <- colnames(predictors)
  if(!any(is.null(og_cols))){
    replace_colnames <- TRUE
  } else {
    replace_colnames <- FALSE
  }

  ## Check nrow of input predictors and matrix coersion
  t <- try({if(nrow(methods::as(predictors,'matrix')) < 3){
    stop('\n \t Need at least 3 observations to fit model \n')
  }}, silent = TRUE)
  if(inherits(t, 'try-error')){
    stop('\n \t Cannot coerce predictors to a matrix \n')
  }

  ## Check if no spline terms - if so, set K = 0
  if(length(unique(c(just_linear_with_interactions,
                     just_linear_without_interactions))) == ncol(predictors)){
    K <- 0
  }

  ## Check if custom knots is not missing, that it can be
  # coerced to a matrix
  if(any(!(is.null(custom_knots)))){
    custom_knots <- try(cbind(custom_knots),
                        silent = TRUE)
    if(any(inherits(custom_knots, 'try-error')) & include_warnings){
      warning('\n \t custom_knots must be a matrix, or should be coercible ',
              'to it. Custom_knots will be ignored. \n')
      custom_knots <- NULL
    }
  }

  ## If ncol(predictors) > 1 and include_quartic_terms is NULL,
  # set include_quartic_terms = TRUE
  # otherwise, if ncol(predictors) == 1 then FALSE
  if(is.null(include_quartic_terms)){
    if(ncol(cbind(predictors)) == 1){
      include_quartic_terms <- FALSE
    } else {
      include_quartic_terms <- TRUE
    }
  }

  ## factor_groups carries one-hot indicator positions back to lgspline.fit.
  #  An empty list means no encoded factor groups were detected.
  list(
    predictors                       = predictors,
    y                                = y,
    og_cols                          = og_cols,
    replace_colnames                 = replace_colnames,
    just_linear_with_interactions    = just_linear_with_interactions,
    just_linear_without_interactions = just_linear_without_interactions,
    exclude_interactions_for         = exclude_interactions_for,
    exclude_these_expansions         = exclude_these_expansions,
    offset                           = offset,
    no_intercept                     = no_intercept,
    do_not_cluster_on_these          = do_not_cluster_on_these,
    observation_weights              = observation_weights,
    K                                = K,
    include_3way_interactions        = include_3way_interactions,
    include_quartic_terms            = include_quartic_terms,
    data                             = data,
    include_constrain_second_deriv   = include_constrain_second_deriv,
    custom_knots                     = custom_knots,
    dummy_fit                        = dummy_fit,
    factor_groups                    = factor_groups
  )
}


