# Recursive Feature Elimination for SCE
# This script implements RFE to identify the most important predictors for SCE models

RFE_SCE <- function(
  Training_data,
  Testing_data,
  Predictors,
  Predictant,
  Nmin,
  Ntree,
  alpha = 0.05,
  resolution = 1000,
  step = 1,  # Number of predictors to remove at each iteration
  verbose = TRUE,  # Control output verbosity
  parallel = TRUE  # Control parallel processing
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
    if (verbose) {
      message("Evaluating model with ", length(current_predictors), " predictors: ", 
              paste(current_predictors, collapse = ", "))
    }
    
    # Train SCE model
    model <- SCE(
      Training_data = Training_data,
      X = current_predictors,
      Y = Predictant,
      mfeature = round(length(current_predictors)/2),
      Ntree = Ntree,
      alpha = alpha,
      Nmin = Nmin,
      resolution = resolution,
      verbose = verbose,
      parallel = parallel
    )
    
    # Get predictions
    predictions <- Model_simulation(
      Testing_data = Testing_data,
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
    least_important <- importance_scores$Predictor[order(importance_scores$Relative_Importance)[1:min(step, length(current_predictors) - length(Predictant))]]
    current_predictors <- setdiff(current_predictors, least_important)
  }
  
  # Return results
  return(history)
}
