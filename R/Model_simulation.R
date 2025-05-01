#################################################################
# Filename: 	Model_simulation.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
Model_simulation <- function(Testing_data, Training_data, model)
{
  # Input validation
  if (is.null(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (is.null(Training_data)) {
    stop("Training_data must be a data frame or matrix")
  }
  
  if (is.null(model)) {
    stop("model must be a list")
  }
  
  if (!is.data.frame(Testing_data) && !is.matrix(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (!is.data.frame(Training_data) && !is.matrix(Training_data)) {
    stop("Training_data must be a data frame or matrix")
  }
  
  if (!is.list(model)) {
    stop("model must be a list")
  }
  
  if (nrow(Testing_data) == 0) {
    stop("Testing_data is empty")
  }
  
  if (nrow(Training_data) == 0) {
    stop("Training_data is empty")
  }
  
  # Get required predictors from the first tree in the model
  required_predictors <- model[[1]]$Tree_Info$Features
  required_predictants <- model[[1]]$YName
  
  # Check if all required predictors are present in testing data
  if (!all(required_predictors %in% colnames(Testing_data))) {
    missing_vars <- setdiff(required_predictors, colnames(Testing_data))
    stop(sprintf("The following predictors are not found in Testing_data: %s", 
                paste(missing_vars, collapse = ", ")))
  }
  
  # Check for missing values in required predictors
  if (any(is.na(Testing_data[, required_predictors]))) {
    stop("Testing_data contains missing values")
  }
  
  # Check for missing values in required predictants
  if (any(is.na(Testing_data[, required_predictants]))) {
    stop("Testing_data contains missing values")
  } 

  # Check data types of required predictors
  if (!all(sapply(Testing_data[, required_predictors], is.numeric))) {
    stop("The following predictors are not numeric")
  }
  
  # Check data types of required predictants
  if (!all(sapply(Testing_data[, required_predictants], is.numeric))) {
    stop("The following predictants are not numeric")
  }
  
  # Get training simulations using SCE_Prediction
  training_sim <- SCE_Prediction(
    X_sample = Training_data,
    model = model
  )
  training_sim <- as.data.frame(training_sim)
  
  # Get validation (OOB) simulations
  validation_sim <- OOB_validation(model)
  
  # Get testing simulations
  testing_sim <- SCE_Prediction(
    X_sample = Testing_data,
    model = model
  )
  testing_sim <- as.data.frame(testing_sim)
  
  # Return results
  output <- list(
    Training = training_sim,
    Validation = validation_sim,
    Testing = testing_sim
  )
  
  return(output)
}

SCE_Prediction <- function(X_sample, model)
{
  # Input validation
  if (is.null(X_sample)) {
    stop("X_sample must be a data frame or matrix")
  }
  
  if (!is.data.frame(X_sample) && !is.matrix(X_sample)) {
    stop("X_sample must be a data frame or matrix")
  }
  
  if (nrow(X_sample) == 0) {
    stop("X_sample is empty")
  }
  
  # Get required predictors from the first tree in the model
  required_predictors <- model[[1]]$Tree_Info$Features
  
  # Check for missing values in required predictors
  if (any(is.na(X_sample[, required_predictors]))) {
    stop("X_sample contains missing values")
  }
  
  # Check data types of required predictors
  if (!all(sapply(X_sample[, required_predictors], is.numeric))) {
    non_numeric <- names(which(!sapply(X_sample[, required_predictors], is.numeric)))
    stop("The following predictors are not numeric")
  }
  
  # Get model predictions for each tree
  predictions <- lapply(model, function(m) {
    SCA_tree_predict(
      Test_data = X_sample,
      model = m
    )
  })
  
  # Combine predictions with their weights
  weighted_predictions <- mapply(
    function(pred, m) {
      pred$Testing_sim * m$weight
    },
    pred = predictions,
    m = model,
    SIMPLIFY = FALSE
  )
  
  # Get number of predictants and their names
  num_predictants <- length(model[[1]]$YName)
  predictant_names <- model[[1]]$YName
  
  # Calculate ensemble predictions
  ensemble_predictions <- matrix(0, nrow = nrow(X_sample), ncol = num_predictants)
  
  for(i in 1:num_predictants) {
    # Extract predictions for current predictant from all trees
    predictant_predictions <- sapply(weighted_predictions, function(x) x[, i])
    
    # Sum weighted predictions
    ensemble_predictions[, i] <- rowSums(predictant_predictions)
  }
  
  colnames(ensemble_predictions) <- predictant_names
  return(ensemble_predictions)
}

OOB_validation <- function(model)
{
  # Input validation
  if (is.null(model)) {
    stop("model must be a list")
  }
  
  if (!is.list(model)) {
    stop("model must be a list")
  }
  
  # Get number of predictants and their names
  num_predictants <- length(model[[1]]$YName)
  predictant_names <- model[[1]]$YName
  
  # Get OOB predictions for each tree
  oob_predictions <- lapply(model, function(x) {
    # Get OOB information from Tree_Info
    oob_indices <- x$Tree_Info$OOB_Indices
    
    # Create data frame with OOB predictions and weights
    data.frame(
      ID = oob_indices,
      x$OOB_sim,  # This will create columns for each predictant
      weight = x$weight
    )
  })
  
  # Combine all OOB predictions
  oob_data <- do.call(rbind, oob_predictions)
  
  # Rename the columns
  colnames(oob_data) <- c("ID", predictant_names, "weight")
  
  # Get unique IDs
  id_list <- sort(unique(oob_data$ID))
  
  # Calculate weighted means for each predictant
  result <- lapply(predictant_names, function(pred_name) {
    # Calculate weighted means for current predictant
    weighted_means <- sapply(id_list, function(id) {
      # Get subset for this ID
      subset <- oob_data[oob_data$ID == id, ]
      # Calculate weighted mean
      sum(subset[[pred_name]] * subset$weight) / sum(subset$weight)
    })
    
    # Create ordered results for this predictant
    data.frame(
      ID = id_list,
      Value = weighted_means
    )[order(id_list), "Value", drop = FALSE]
  })
  
  # Combine results for all predictants
  result <- do.call(cbind, result)
  colnames(result) <- predictant_names
  
  return(result)
}