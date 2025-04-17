# Recursive Feature Elimination for SCE
# This script implements RFE to identify the most important predictors for SCE models

RFE_SCE <- function(
  Training_data,
  Testing_data,
  Predictors,
  Predictant,
  alpha = 0.05,
  Nmin = 5,
  Ntree = 40,
  resolution = 50,
  metric = "rmse",  # Can be "rmse", "mae", "nse", "log_nse", "R2", "kge"
  step = 1  # Number of predictors to remove at each iteration
) {
  # Input validation
  if (!is.data.frame(Training_data) || !is.data.frame(Testing_data)) {
    stop("All data inputs must be data frames")
  }
  
  if (!is.character(Predictors) || length(Predictors) < 2) {
    stop("Predictors must be a character vector with at least 2 elements")
  }
  
  if (!is.character(Predictant) || length(Predictant) == 0) {
    stop("Predictant must be a non-empty character vector")
  }
  
  if (!is.numeric(step) || step < 1 || step > length(Predictors) - length(Predictant)) {
    stop("step must be a positive integer less than or equal to the number of predictors minus the number of predictants")
  }
  
  # Define which metrics should be maximized vs minimized
  maximize_metrics <- c("nse", "log_nse", "R2", "kge")
  minimize_metrics <- c("rmse", "mae")
  
  if (!metric %in% c(maximize_metrics, minimize_metrics)) {
    stop("Invalid metric. Must be one of: rmse, mae, nse, log_nse, R2, kge")
  }
  
  # if the number of Predictant is more than 1, the metric must be one of "nse", "log_nse", "R2", "kge"
  if (length(Predictant) > 1 && !metric %in% maximize_metrics) {
    stop("For multiple predictants, the metric must be one of: nse, log_nse, R2, kge")
  }
  
  # Initialize variables
  current_predictors <- Predictors
  history <- list(
    summary = data.frame(
      n_predictors = integer(),
      predictors = character(),
      stringsAsFactors = FALSE
    ),
    performances = list(),
    importance_scores = list()
  )
  
  # Main RFE loop
  while (length(current_predictors) > (length(Predictant)+2) ) {
    cat("\nEvaluating model with", length(current_predictors), "predictors:", 
        paste(current_predictors, collapse = ", "), "\n")
    
    # Train SCE model
    model <- SCE(
      Training_data = Training_data,
      X = current_predictors,
      Y = Predictant,
      mfeature = round(0.5*length(current_predictors)),
      Ntree = Ntree,
      alpha = alpha,
      Nmin = Nmin,
      resolution = resolution
    )
    
    # Get predictions
    predictions <- Model_simulation(
      Testing_data = Testing_data,
      Training_data = Training_data,
      model = model
    )
    
    # Evaluate model
    evaluation <- SCE_Model_evaluation(
      Testing_data = Testing_data,
      Training_data = Training_data,
      Simulations = predictions,
      Predictant = Predictant,
      digits = 3
    )
    
    # Store summary and performance
    history$summary <- rbind(history$summary, data.frame(
      n_predictors = length(current_predictors),
      predictors = paste(current_predictors, collapse = ","),
      current_metric = current_metric,
      stringsAsFactors = FALSE
    ))
    
    # Store performance data frames
    if (length(Predictant) == 1) {
      history$performances[[length(history$performances) + 1]] <- evaluation
    } else {
      # For multiple predictants, store each predictant's performance
      for (pred in names(evaluation)) {
        history$performances[[paste0("n_", length(current_predictors), "_", pred)]] <- evaluation[[pred]]
      }
    }
    
    # Calculate importance scores
    importance_scores <- Wilks_importance(
      model = model
    )
    history$importance_scores[[length(history$importance_scores) + 1]] <- importance_scores

    # Remove step number of least important predictors
    least_important <- importance_scores$col_index[order(importance_scores$Importance)[1:min(step, length(current_predictors) - length(Predictant))]]
    current_predictors <- setdiff(current_predictors, least_important)
  }
  
  # Return results
  return(list(
    history = history,
    final_model = model
  ))
}

# Example usage:
# result <- RFE_SCE(
#   Training_data = train_data,
#   Testing_data = test_data,
#   Predictors = c("x1", "x2", "x3", "x4"),
#   Predictant = "y",
#   alpha = 0.05,
#   Nmin = 5,
#   resolution = 50,
#   metric = "rmse",
#   step = 2  # Remove 2 predictors at a time
# ) 
