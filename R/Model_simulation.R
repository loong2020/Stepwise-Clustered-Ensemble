#################################################################
# Filename: 	Model_simulation.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
model_simulation <- function(model, Testing_data)
{
  # Input validation
  if (is.null(model)) {
    stop("model must be an SCE object or list")
  }
  
  if (is.null(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (!inherits(model, "sce") && !is.list(model)) {
    stop("model must be an SCE object or list")
  }
  
  if (!is.data.frame(Testing_data) && !is.matrix(Testing_data)) {
    stop("Testing_data must be a data frame or matrix")
  }
  
  if (nrow(Testing_data) == 0) {
    stop("Testing_data is empty")
  }
  
  # Handle S3 class objects
  if (inherits(model, "sce")) {
    # Convert SCE object to list format for compatibility
    model_list <- list()
    for (i in seq_along(model$trees)) {
      model_list[[i]] <- model$trees[[i]]
    }
  } else {
    # Legacy list format
    model_list <- model
  }
  
  # Get training simulations using sce_prediction
  training_sim <- training_prediction(model_list)
  
  # Get validation (OOB) simulations
  validation_sim <- oob_validation(model_list)
  
  # Get testing simulations
  testing_sim <- sce_prediction(
    model = model_list,
    X_sample = Testing_data
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

sce_prediction <- function(model, X_sample)
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
    sca_tree_predict(
      model = m,
      Testing_data = X_sample
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

training_prediction <- function(model)
{

  # Get number of predictants and their names
  num_predictants <- length(model[[1]]$YName)
  predictant_names <- model[[1]]$YName
  
  # Get model predictions for each tree
  predictions <- lapply(model, function(m) {
    sca_tree_predict(
      model = m,
      Testing_data = m$Training_data
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

oob_validation <- function(model)
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