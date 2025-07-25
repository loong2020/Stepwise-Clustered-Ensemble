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
      validation_r2 = numeric(),
      testing_r2 = numeric(),
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
      model = model,
      Testing_data = Testing_data
    )
    
    # Evaluate model
    evaluation <- SCE_Model_evaluation(
      Testing_data = Testing_data,
      Training_data = Training_data,
      Simulations = predictions,
      Predictant = Predictant,
      digits = 3
    )
    
    # Extract R2 values from evaluation
    validation_r2 <- if (is.data.frame(evaluation) && "Validation" %in% colnames(evaluation) && "R2" %in% rownames(evaluation)) {
      as.numeric(evaluation["R2", "Validation"])
    } else {
      NA
    }
    
    testing_r2 <- if (is.data.frame(evaluation) && "Testing" %in% colnames(evaluation) && "R2" %in% rownames(evaluation)) {
      as.numeric(evaluation["R2", "Testing"])
    } else {
      NA
    }
    
    # Store summary and performance
    history$summary <- rbind(history$summary, data.frame(
      n_predictors = length(current_predictors),
      predictors = paste(current_predictors, collapse = ","),
      validation_r2 = validation_r2,
      testing_r2 = testing_r2,
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

# Plot RFE Results
# 
# Creates a plot showing validation and testing R2 values as a function of the number of predictors
# during recursive feature elimination. Uses the R2 values stored directly in the summary dataframe.
Plot_RFE <- function(
  rfe_result,
  main = "OOB Validation and Testing R2 vs Number of Predictors",
  col_validation = "blue",
  col_testing = "red",
  pch = 16,
  lwd = 2,
  cex = 1.2,
  legend_pos = "bottomleft",
  ...
) {
  # Input validation
  if (!is.list(rfe_result) || !all(c("summary") %in% names(rfe_result))) {
    stop("rfe_result must be a list with 'summary' component from RFE_SCE function")
  }
  
  # Check if summary has the required columns
  if (!all(c("n_predictors", "validation_r2", "testing_r2") %in% colnames(rfe_result$summary))) {
    stop("rfe_result$summary must contain 'n_predictors', 'validation_r2', and 'testing_r2' columns")
  }
  
  # Extract data directly from summary
  n_predictors <- rfe_result$summary$n_predictors
  validation_r2 <- as.numeric(rfe_result$summary$validation_r2)
  testing_r2 <- as.numeric(rfe_result$summary$testing_r2)
  
  # Check for valid data
  if (all(is.na(validation_r2)) && all(is.na(testing_r2))) {
    stop("No valid R2 values found in the RFE results")
  }
  
  # Calculate y-axis limits
  y_values <- c(validation_r2, testing_r2)
  y_values <- y_values[!is.na(y_values)]
  ylim <- c(min(y_values), max(y_values))
  
  # Create the plot
  plot(n_predictors, validation_r2, 
       type = "b",  # both points and lines
       col = col_validation,
       pch = pch,
       lwd = lwd,
       cex = cex,
       xlim = rev(range(n_predictors)),  # reverse x-axis
       ylim = ylim,
       xlab = "Number of Predictors",
       ylab = "R2",
       main = main,
       ...)
  
  # Add testing data
  lines(n_predictors, testing_r2, type = "b", col = col_testing, pch = pch, lwd = lwd, cex = cex)
  
  # Add legend
  legend(legend_pos,
         legend = c("OOB Validation", "Testing"),
         col = c(col_validation, col_testing),
         pch = pch,
         lty = 1,
         lwd = lwd)
  
  # Return plot data invisibly
  invisible(list(
    n_predictors = n_predictors,
    validation_r2 = validation_r2,
    testing_r2 = testing_r2
  ))
}
