# S3 class constructor for SCA objects
SCA_object <- function(tree, map, predictors, predictants, type, total_nodes, leaf_nodes, cutting_actions, merging_actions, call) {
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
    class = "SCA"
  )
}

# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
SCA <- function(Training_data, X, Y, Nmin, alpha = 0.05, resolution = 1000, verbose = FALSE)
{
  # Store the function call
  call <- match.call()
  
  #: store the start time
  time_stat <- proc.time()

  #: Input validation
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  if (!is.data.frame(Training_data) && !is.matrix(Training_data)) {
    stop("Training_data must be a data frame or matrix")
  }
  
  if (nrow(Training_data) == 0) {
    stop("Training_data is empty")
  }
  
  if (!all(X %in% colnames(Training_data))) {
    missing_vars <- setdiff(X, colnames(Training_data))
    stop(sprintf("The following predictors are not found in Training_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  if (!all(Y %in% colnames(Training_data))) {
    missing_vars <- setdiff(Y, colnames(Training_data))
    stop(sprintf("The following predictants are not found in Training_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  if (!is.numeric(Nmin) || Nmin <= 0) {
    stop("Nmin must be a positive number")
  }
  
  if (!is.numeric(resolution) || resolution <= 0) {
    stop("resolution must be a positive number")
  }
  
  # Check for missing values
  if (any(is.na(Training_data[, X])) || any(is.na(Training_data[, Y]))) {
    stop("Training_data contains missing values")
  }
  
  # Check data types
  if (!all(sapply(Training_data[, X], is.numeric))) {
    non_numeric <- names(which(!sapply(Training_data[, X], is.numeric)))
    stop(sprintf("The following predictors are not numeric: %s", 
                paste(non_numeric, collapse = ", ")))
  }
  
  if (!all(sapply(Training_data[, Y], is.numeric))) {
    non_numeric <- names(which(!sapply(Training_data[, Y], is.numeric)))
    stop(sprintf("The following predictants are not numeric: %s", 
                paste(non_numeric, collapse = ", ")))
  }

  #: create data structure
  data <- list()
  
  #: store input data
  data$o_sample_data_x <- as.matrix(Training_data[, X])
  data$o_sample_data_y <- as.matrix(Training_data[, Y])
  
  #: store dimensions
  data$n_sample_size <- nrow(data$o_sample_data_x)
  data$n_sample_x_cols <- ncol(data$o_sample_data_x)
  data$n_sample_y_cols <- ncol(data$o_sample_data_y)
  
  # Check minimum sample size
  if (data$n_sample_size < Nmin) {
    stop(sprintf("Sample size (%d) is less than Nmin (%d)", 
                data$n_sample_size, Nmin))
  }
  
  #: store parameters
  data$n_alpha <- alpha
  data$n_mapvalue <- "mean"  
  data$resolution <- resolution
  data$Nmin <- Nmin
  
  #: do clustering
  result <- do_cluster(data = data, Nmin = Nmin, resolution = resolution, verbose = verbose)

  #: return the S3 class object
  return(SCA_object(
    tree = result$Tree,
    map = result$Map,
    predictors = X,
    predictants = Y,
    type = data$n_mapvalue,
    total_nodes = result$totalNodes,
    leaf_nodes = result$leafNodes,
    cutting_actions = result$cuttingActions,
    merging_actions = result$mergingActions,
    call = call
  ))
}

SCA_tree_predict <- function(model, Testing_data) {
  # Input validation
  if (is.null(model)) {
    stop("model must be an SCA object or list")
  }
  
  if (is.null(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (!is.data.frame(Testing_data) && !is.matrix(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (nrow(Testing_data) == 0) {
    stop("Testing_data is empty")
  }
  
  # Check if all required predictors are present in test data
  if (!all(model$XName %in% colnames(Testing_data))) {
    missing_vars <- setdiff(model$XName, colnames(Testing_data))
    stop(sprintf("The following predictors are not found in Testing_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  # Initialize Test_X
  Test_X <- Testing_data[,model$XName]
  
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
print.SCA <- function(x, ...) {
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
summary.SCA <- function(object, ...) {
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
predict.SCA <- function(object, newdata, ...) {
  # This is a wrapper for SCA_tree_predict
  if (missing(newdata)) {
    stop("newdata is required for prediction")
  }
  
  return(SCA_tree_predict(model = object, Testing_data = newdata))
}

# Importance method for SCA objects
importance.SCA <- function(object, digits = 2, ...) {
  # This is a wrapper for SCA_importance
  return(SCA_importance(model = object, digits = digits))
}

# Evaluate method for SCA objects
evaluate.SCA <- function(object, Testing_data, Training_data, digits = 3, ...) {
  # This is a wrapper for SCA_Model_evaluation
  if (missing(Testing_data)) {
    stop("Testing_data is required for evaluation")
  }
  if (missing(Training_data)) {
    stop("Training_data is required for evaluation")
  }
  
  # Get predictants from the object
  Predictant <- object$YName
  
  # Check for extra parameters that are not needed for SCA
  args <- list(...)
  if (length(args) > 0) {
    # Check for any extra parameters
    warning("Extra parameters were provided but are not used for SCA evaluation: ", 
            paste(names(args), collapse = ", "))
  }
  
  # Get predictions using SCA_tree_predict
  predictions_testing <- SCA_tree_predict(model = object, Testing_data = Testing_data)
  predictions_training <- SCA_tree_predict(model = object, Testing_data = Training_data)
  
  Testing_GOF <- SCA_Model_evaluation(
    Testing_data = Testing_data,
    Simulations = predictions_testing,
    Predictant = Predictant,
    digits = digits
  )
  
  Training_GOF <- SCA_Model_evaluation(
    Testing_data = Training_data,
    Simulations = predictions_training,
    Predictant = Predictant,
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