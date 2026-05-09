# S3 class constructor for SCA objects
sca_object <- function(tree, map, predictors, predictants, type, total_nodes, leaf_nodes, cutting_actions, merging_actions, call) {
  structure(
    list(
      Tree = tree,
      Map = map,
      XName = predictors,
      YName = predictants,
      type = type,
      totalNodes = total_nodes,
      leafNodes = leaf_nodes,
      cuttingActions = cutting_actions,
      mergingActions = merging_actions,
      call = call
    ),
    class = "sca"
  )
}

# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
sca <- function(training_data, x, y, nmin, alpha = 0.05, resolution = 1000, verbose = FALSE)
{
  # Store the function call
  call <- match.call()
  
  #: store the start time
  time_stat <- proc.time()

  #: Input validation
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  if (!is.data.frame(training_data) && !is.matrix(training_data)) {
    stop("training_data must be a data frame or matrix")
  }
  
  if (nrow(training_data) == 0) {
    stop("training_data is empty")
  }
  
  if (!all(x %in% colnames(training_data))) {
    missing_vars <- setdiff(x, colnames(training_data))
    stop(sprintf("The following predictors are not found in training_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  if (!all(y %in% colnames(training_data))) {
    missing_vars <- setdiff(y, colnames(training_data))
    stop(sprintf("The following predictants are not found in training_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  if (!is.numeric(nmin) || nmin <= 0) {
    stop("nmin must be a positive number")
  }
  
  if (!is.numeric(resolution) || resolution <= 0) {
    stop("resolution must be a positive number")
  }
  
  # Check for missing values
  if (any(is.na(training_data[, x])) || any(is.na(training_data[, y]))) {
    stop("training_data contains missing values")
  }
  
  # Check data types
  if (!all(sapply(training_data[, x], is.numeric))) {
    non_numeric <- names(which(!sapply(training_data[, x], is.numeric)))
    stop(sprintf("The following predictors are not numeric: %s", 
                paste(non_numeric, collapse = ", ")))
  }
  
  if (!all(sapply(training_data[, y], is.numeric))) {
    non_numeric <- names(which(!sapply(training_data[, y], is.numeric)))
    stop(sprintf("The following predictants are not numeric: %s", 
                paste(non_numeric, collapse = ", ")))
  }

  #: create data structure
  data <- list()
  
  #: store input data
  data$o_sample_data_x <- as.matrix(training_data[, x])
  data$o_sample_data_y <- as.matrix(training_data[, y])
  
  #: store dimensions
  data$n_sample_size <- nrow(data$o_sample_data_x)
  data$n_sample_x_cols <- ncol(data$o_sample_data_x)
  data$n_sample_y_cols <- ncol(data$o_sample_data_y)
  
  # Check minimum sample size
  if (data$n_sample_size < nmin) {
    stop(sprintf("Sample size (%d) is less than nmin (%d)", 
                data$n_sample_size, nmin))
  }
  
  #: store parameters
  data$n_alpha <- alpha
  data$n_mapvalue <- "mean"  
  data$resolution <- resolution
  data$nmin <- nmin
  
  #: do clustering
  result <- do_cluster(data = data, nmin = nmin, resolution = resolution, verbose = verbose)

  #: return the S3 class object
  return(sca_object(
    tree = result$Tree,
    map = result$Map,
    predictors = x,
    predictants = y,
    type = data$n_mapvalue,
    total_nodes = result$totalNodes,
    leaf_nodes = result$leafNodes,
    cutting_actions = result$cuttingActions,
    merging_actions = result$mergingActions,
    call = call
  ))
}

sca_tree_predict <- function(model, testing_data) {
  # Input validation
  if (is.null(model)) {
    stop("model must be an sca object or list")
  }
  
  if (is.null(testing_data)) {
    stop("testing_data must be a data frame or matrix")
  }
  
  if (!is.data.frame(testing_data) && !is.matrix(testing_data)) {
    stop("testing_data must be a data frame or matrix")
  }
  
  if (nrow(testing_data) == 0) {
    stop("testing_data is empty")
  }
  
  # Check if all required predictors are present in test data
  if (!all(model$XName %in% colnames(testing_data))) {
    missing_vars <- setdiff(model$XName, colnames(testing_data))
    stop(sprintf("The following predictors are not found in testing_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  # Initialize Test_X
  Test_X <- testing_data[, model$XName]
  
  # Initialize data structure
  data <- list()
  
  # Store input data
  data$o_sample_data_x <- as.matrix(na.omit(Test_X))
  data$n_sample_size <- nrow(data$o_sample_data_x)
  data$n_sample_x_cols <- ncol(data$o_sample_data_x)
  
  # Store model data
  data$o_result_tree <- model$Tree
  data$n_result_tree_rows <- nrow(data$o_result_tree)
  data$o_mean_y <- model$Map
  data$n_y_cols <- ncol(data$o_mean_y)
  data$n_model_type <- model$type
  
  # Initialize prediction matrix
  data$o_predictants <- matrix(0, data$n_sample_size, data$n_y_cols)
  
  # Perform prediction
  data = f_main_p(data)
  Testing_sim <- data$o_predictants
  
  Testing_sim <- data.frame(do.call(cbind, data.frame(Testing_sim)))
  colnames(Testing_sim) <- model$YName
  
  return(Testing_sim)
}

# S3 Methods for SCA class

# Print method for SCA objects
print.sca <- function(x, ...) {
  cat("Stepwise Cluster Analysis (SCA) Model\n")
  cat("=====================================\n\n")
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Model Structure:\n")
  cat("  Total nodes:", x$totalNodes, "\n")
  cat("  Leaf nodes:", x$leafNodes, "\n")
  cat("  Cutting actions:", x$cuttingActions, "\n")
  cat("  Merging actions:", x$mergingActions, "\n")
  cat("  Mapping type:", x$type, "\n\n")
  
  cat("Variables:\n")
  cat("  Predictors:", paste(x$XName, collapse = ", "), "\n")
  cat("  Predictants:", paste(x$YName, collapse = ", "), "\n")
  
  invisible(x)
}

# Summary method for SCA objects
summary.sca <- function(object, ...) {
  cat("Stepwise Cluster Analysis (SCA) Model Summary\n")
  cat("============================================\n\n")
  
  cat("Model Information:\n")
  cat("  Total nodes:", object$totalNodes, "\n")
  cat("  Leaf nodes:", object$leafNodes, "\n")
  cat("  Internal nodes:", object$totalNodes - object$leafNodes, "\n")
  cat("  Cutting actions:", object$cuttingActions, "\n")
  cat("  Merging actions:", object$mergingActions, "\n")
  cat("  Mapping type:", object$type, "\n\n")
  
  cat("Variables:\n")
  cat("  Number of predictors:", length(object$XName), "\n")
  cat("  Number of predictants:", length(object$YName), "\n")
  cat("  Predictors:", paste(object$XName, collapse = ", "), "\n")
  cat("  Predictants:", paste(object$YName, collapse = ", "), "\n\n")
  
  # Tree depth calculation
  if (object$totalNodes > 0) {
    estimated_depth <- ceiling(log2(object$totalNodes))
    cat("Tree Characteristics:\n")
    cat("  Estimated maximum depth:", estimated_depth, "\n")
    cat("  Average branching factor:", round(object$totalNodes / max(1, object$totalNodes - object$leafNodes), 2), "\n")
  }
  
  invisible(object)
}

# Predict method for SCA objects
predict.sca <- function(object, newdata, ...) {
  # This is a wrapper for sca_tree_predict
  if (missing(newdata)) {
    stop("newdata is required for prediction")
  }
  
  return(sca_tree_predict(model = object, testing_data = newdata))
}

# Importance method for SCA objects
importance.sca <- function(object, digits = 2, ...) {
  # This is a wrapper for sca_importance
  return(sca_importance(model = object, digits = digits))
}

# Evaluate method for SCA objects
evaluate.sca <- function(object, testing_data, training_data, digits = 3, ...) {
  # This is a wrapper for sca_model_evaluation
  if (missing(testing_data)) {
    stop("testing_data is required for evaluation")
  }
  if (missing(training_data)) {
    stop("training_data is required for evaluation")
  }
  
  # Get predictants from the object
  predictant <- object$YName
  
  # Check for extra parameters that are not needed for SCA
  args <- list(...)
  if (length(args) > 0) {
    # Check for any extra parameters
    warning("Extra parameters were provided but are not used for SCA evaluation: ", 
            paste(names(args), collapse = ", "))
  }
  
  # Get predictions using sca_tree_predict
  predictions_testing <- sca_tree_predict(model = object, testing_data = testing_data)
  predictions_training <- sca_tree_predict(model = object, testing_data = training_data)
  
  Testing_GOF <- sca_model_evaluation(
    testing_data = testing_data,
    simulations = predictions_testing,
    predictant = predictant,
    digits = digits
  )
  
  Training_GOF <- sca_model_evaluation(
    testing_data = training_data,
    simulations = predictions_training,
    predictant = predictant,
    digits = digits
  )

  # Create result dataframe with proper rownames
  result_df <- data.frame(
    Training = Training_GOF$Testing,
    Testing = Testing_GOF$Testing
  )
  
  # Preserve the rownames from the evaluation results
  rownames(result_df) <- rownames(Training_GOF)
  
  return(result_df)
}
