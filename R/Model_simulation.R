#################################################################
# Filename: 	Model_simulation.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
Model_simulation <- function(Testing_data, model)
{
  # Input validation
  if (is.null(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (is.null(model)) {
    stop("model must be a list")
  }
  
  if (!is.data.frame(Testing_data) && !is.matrix(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (!is.list(model)) {
    stop("model must be a list")
  }
  
  if (nrow(Testing_data) == 0) {
    stop("Testing_data is empty")
  }
  
  # Get training simulations using SCE_Prediction
  training_sim <- Training_Prediction(model)
  
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
  
  # Get model predictions for each tree
  predictions <- lapply(model, function(m) {
    SCA_tree_predict(
      Testing_data = X_sample,
      model = m
    )
  })
  
  # Combine predictions with their weights
  weighted_predictions <- mapply(
    function(pred, m) {
      pred * m$weight
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

Training_Prediction <- function(model)
{

  # Get number of predictants and their names
  num_predictants <- length(model[[1]]$YName)
  predictant_names <- model[[1]]$YName
  
  # Get model predictions for each tree
  predictions <- lapply(model, function(m) {
    SCA_tree_predict(
      Testing_data = m$Training_data,
      model = m
    )
  })

  #assign the ID for each data point
  predictions <- mapply(
    function(pred, m) {
      data.frame(
        ID = m$Sample,
        pred = pred,
        weight = m$weight
      )
    },
    pred = predictions,
    m = model,
    SIMPLIFY = FALSE
  )

    # Combine all OOB predictions
  combined_predictions <- do.call(rbind, predictions)
  
  # Rename the columns
  colnames(combined_predictions) <- c("ID", predictant_names, "weight")
  
  # Get unique IDs
  id_list <- sort(unique(combined_predictions$ID))
  
  # Calculate weighted means for each predictant
  result <- lapply(predictant_names, function(pred_name) {
    # Calculate weighted means for current predictant
    weighted_means <- sapply(id_list, function(id) {
      # Get subset for this ID
      subset <- combined_predictions[combined_predictions$ID == id, ]
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