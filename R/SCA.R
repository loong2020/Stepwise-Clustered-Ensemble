# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
SCA <- function(Training_data, X, Y, Nmin, alpha = 0.05, resolution = 1000)
{
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
  result <- do_cluster(data = data, Nmin = Nmin, resolution = resolution)

  #: return the model
  model <- list(
    Tree = result$Tree,
    Map = result$Map,
    XName = X,
    YName = Y,
    type = data$n_mapvalue,
    totalNodes = result$totalNodes,
    leafNodes = result$leafNodes,
    cuttingActions = result$cuttingActions,
    mergingActions = result$mergingActions
  )
  
  return(model)
}

SCA_tree_predict <- function(Testing_data, model) {
  # Input validation
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
  
  model$Testing_sim <- Testing_sim
  
  return(model)
}