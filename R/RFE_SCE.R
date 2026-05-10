# Recursive Feature Elimination for SCE
# This script implements RFE to identify the most important predictors for SCE models

rfe_sce <- function(
  training_data,
  testing_data,
  predictors,
  predictant,
  nmin,
  ntree,
  alpha = 0.05,
  resolution = 1000,
  step = 1,  # Number of predictors to remove at each iteration
  verbose = TRUE,  # Control output verbosity
  parallel = TRUE  # Control parallel processing
) {
  # Input validation
  if (!is.data.frame(training_data) || !is.data.frame(testing_data)) {
    stop("All data inputs must be data frames")
  }
  
  if (!is.character(predictors) || length(predictors) < 2) {
    stop("predictors must be a character vector with at least 2 elements")
  }
  
  if (!is.character(predictant) || length(predictant) == 0) {
    stop("predictant must be a non-empty character vector")
  }
  
  if (!is.numeric(step) || step < 1 || step > length(predictors) - length(predictant)) {
    stop("step must be a positive integer less than or equal to the number of predictors minus the number of predictants")
  }
  
  # Initialize variables
  current_predictors <- predictors
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
  while (length(current_predictors) > (length(predictant)+2) ) {
    if (verbose) {
      message("Evaluating model with ", length(current_predictors), " predictors: ", 
              paste(current_predictors, collapse = ", "))
    }
    
    # Train SCE model
    model <- sce(
      training_data = training_data,
      x = current_predictors,
      y = predictant,
      mfeature = round(length(current_predictors)/2),
      ntree = ntree,
      alpha = alpha,
      nmin = nmin,
      resolution = resolution,
      verbose = verbose,
      parallel = parallel
    )
    
    # Get predictions
    predictions <- model_simulation(
      model = model,
      testing_data = testing_data
    )
    
    # Evaluate model
    evaluation <- sce_model_evaluation(
      testing_data = testing_data,
      training_data = training_data,
      simulations = predictions,
      predictant = predictant,
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
    if (length(predictant) == 1) {
      history$performances[[length(history$performances) + 1]] <- evaluation
    } else {
      # For multiple predictants, store each predictant's performance
      for (pred in names(evaluation)) {
        history$performances[[paste0("n_", length(current_predictors), "_", pred)]] <- evaluation[[pred]]
      }
    }
    
    # Calculate importance scores
    importance_scores <- wilks_importance(
      model = model
    )
    history$importance_scores[[length(history$importance_scores) + 1]] <- importance_scores

    # Remove step number of least important predictors
    least_important <- importance_scores$Predictor[order(importance_scores$Relative_Importance)[1:min(step, length(current_predictors) - length(predictant))]]
    current_predictors <- setdiff(current_predictors, least_important)
  }
  
  # Return results
  return(history)
}

# Plot RFE Results
# 
# Creates a plot showing validation and testing R2 values as a function of the number of predictors
# during recursive feature elimination. Uses the R2 values stored directly in the summary dataframe.
plot_rfe <- function(
  rfe_result,
  main = NULL,
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
    stop("rfe_result must be a list with 'summary' component from rfe_sce function")
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
       ylab = expression(R^2),
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
