#################################################################
# Filename: 	SCE.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
# History: 	2019/05/17		created by Kailong Li
#			      2019/09/18		added Wilks feature importance (WFI), by Kailong Li
#	     	    2023/03/10   	enabled weighted and non-weighted options for calculating WFI, by Kailong Li
#           2024/12/19    added S3 class support, by Kailong Li
##################################################################

# reference:
# Li, Kailong, Guohe Huang, and Brian Baetz
# Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling
# Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.
# https://doi.org/10.5194/hess-25-4947-2021

# ---------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------

# S3 class constructor for SCE objects
SCE_object <- function(trees, predictors, predictants, parameters, call) {
  structure(
    list(
      trees = trees,
      predictors = predictors,
      predictants = predictants,
      parameters = parameters,
      call = call
    ),
    class = "SCE"
  )
}

SCE <- function(Training_data, X, Y, mfeature, Nmin, Ntree, alpha = 0.05, resolution = 1000, verbose = FALSE, parallel = TRUE) {
  # Store the function call
  call <- match.call()
  
  # Input validation
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  if (!is.numeric(Nmin) || Nmin <= 0) {
    stop("Nmin must be a positive number")
  }
  
  if (!is.numeric(Ntree) || Ntree <= 0) {
    stop("Ntree must be a positive number")
  }
  
  if (!is.numeric(mfeature) || mfeature <= 0 || mfeature > length(X)) {
    stop(sprintf("mfeature must be between 1 and number of predictors (%d)", length(X)))
  }
  
  if (!is.numeric(resolution) || resolution <= 0) {
    stop("resolution must be a positive number")
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

  # Prepare data (following original structure)
  o_xdata <- as.data.frame(Training_data[, X, drop = FALSE])
  o_ydata <- as.data.frame(Training_data[, Y, drop = FALSE])
  
  # Check if we have enough data
  if (nrow(o_xdata) < Nmin) {
    stop(sprintf("Sample size (%d) is less than Nmin (%d)", 
                nrow(o_xdata), Nmin))
  }

  colnames(o_xdata) <- X
  colnames(o_ydata) <- Y

  # Setup bootstrap and random features (following original approach)
  n_predictors <- ncol(o_xdata)
  n_samples <- nrow(o_xdata)
  
  # Generate random features for all trees at once
  Random_col_matrix <- replicate(Ntree, sort(sample(seq_len(n_predictors), mfeature)))
  Random_col <- split(Random_col_matrix, col(Random_col_matrix))
  
  # Generate bootstrap samples for all trees at once
  tree_list <- replicate(Ntree, sample(seq_len(n_samples), replace = TRUE), simplify = FALSE)
  
  # Create bootstrap list
  Bootst_rep <- Map(function(x, y) list(
    Tree = paste("SCA", x, sep = "_"),
    mfeature = Random_col[[x]],
    sample = tree_list[[x]]
  ), x = seq_len(Ntree))
  
  # Setup parallel processing
  if (parallel) {
    available_cores <- parallel::detectCores()
    if (available_cores == 1) {
      max_cores <- 1
    } else {
      max_cores <- min(available_cores, Ntree)  
    }
    
    if (max_cores > 1) {
      Clus <- parallel::makeCluster(max_cores)
      
      # Export required functions and objects to workers
      parallel::clusterExport(Clus, 
        c("o_xdata", "o_ydata", "Nmin", "alpha", "resolution", "verbose", "find_best_split_iterative", "find_best_split",
          "SCA", "SCA_object", "f_processnode", "f_min_wilks", "f_wilks_statistic",
          "f_cal_chk_f", "f_checkif_leaf", "f_init", "f_main", "do_cluster",
          "SCA_tree_predict", "f_main_p", "f_predict_one", "f_predict", "Inference"),
        envir = environment()
      )
      
      # Load required packages in workers
      parallel::clusterEvalQ(Clus, {
        library(parallel)
      })

      # Parallel processing
      SCE_res <- parallel::parLapply(Clus, Bootst_rep, function(rep) {
        # Get feature names for this tree
        feature_names <- colnames(o_xdata)[rep$mfeature]
        
        # Prepare data for this tree
        tree_data <- cbind(
          o_xdata[rep$sample, rep$mfeature, drop = FALSE],
          o_ydata[rep$sample, , drop = FALSE]
        )
        colnames(tree_data) <- c(feature_names, Y)
        
        # Store the tree data
        tree_info <- list(
          Tree = rep$Tree,
          Features = feature_names,
          Sample_Indices = rep$sample
        )
        
        # Run SCA
        tree_model <- SCA(alpha = alpha, Nmin = Nmin, resolution = resolution, 
                         Training_data = tree_data,
                         X = feature_names,
                         Y = Y,
                         verbose = verbose)
        
        # Calculate OOB error
        all_samples <- seq_len(n_samples)
        sample_counts <- table(factor(rep$sample, levels = all_samples))
        oob_indices <- which(sample_counts == 0)
        
        # Prepare OOB data
        oob_xdata <- o_xdata[oob_indices, rep$mfeature, drop = FALSE]
        colnames(oob_xdata) <- feature_names
        oob_ydata <- o_ydata[oob_indices, , drop = FALSE]
        
        # Store OOB data
        tree_info$OOB_Indices <- oob_indices
        tree_info$OOB_XData <- oob_xdata
        tree_info$OOB_YData <- oob_ydata
        
        # Make predictions on OOB data
        oob_predictions <- SCA_tree_predict(
          model = tree_model,
          Testing_data = oob_xdata
        )
        
        # Calculate R-squared for OOB predictions
        if (ncol(oob_ydata) == 1) {
          oob_ydata_numeric <- as.numeric(oob_ydata[[1]])
          oob_predictions_numeric <- as.numeric(oob_predictions[[1]])
          oob_r2 <- 1 - sum((oob_ydata_numeric - oob_predictions_numeric)^2) / 
                    sum((oob_ydata_numeric - mean(oob_ydata_numeric))^2)
        } else {
          oob_r2 <- mean(sapply(1:ncol(oob_ydata), function(i) {
            oob_ydata_numeric <- as.numeric(oob_ydata[,i])
            oob_predictions_numeric <- as.numeric(oob_predictions[,i])
            1 - sum((oob_ydata_numeric - oob_predictions_numeric)^2) / 
                sum((oob_ydata_numeric - mean(oob_ydata_numeric))^2)
          }))
        }
        
        # Add OOB error to model
        tree_model$OOB_error <- oob_r2
        tree_model$OOB_sim <- oob_predictions
        tree_model$Sample <- rep$sample
        tree_model$Tree_Info <- tree_info
        tree_model$Training_data <- tree_data  # Add training data to output
        return(tree_model)
      })
      
      # Clean up
      parallel::stopCluster(Clus)
    }
  } else {
    # Sequential processing when parallel = FALSE
    SCE_res <- lapply(Bootst_rep, function(rep) {
      # Get feature names for this tree
      feature_names <- colnames(o_xdata)[rep$mfeature]
      
      # Prepare data for this tree
      tree_data <- cbind(
        o_xdata[rep$sample, rep$mfeature, drop = FALSE],
        o_ydata[rep$sample, , drop = FALSE]
      )
      colnames(tree_data) <- c(feature_names, Y)
      
      # Store the tree data
      tree_info <- list(
        Tree = rep$Tree,
        Features = feature_names,
        Sample_Indices = rep$sample
      )
      
      # Run SCA
      tree_model <- SCA(alpha = alpha, Nmin = Nmin, resolution = resolution, 
                       Training_data = tree_data,
                       X = feature_names,
                       Y = Y,
                       verbose = verbose)
      
      # Calculate OOB error
      all_samples <- seq_len(n_samples)
      sample_counts <- table(factor(rep$sample, levels = all_samples))
      oob_indices <- which(sample_counts == 0)
      
      # Prepare OOB data
      oob_xdata <- o_xdata[oob_indices, rep$mfeature, drop = FALSE]
      colnames(oob_xdata) <- feature_names
      oob_ydata <- o_ydata[oob_indices, , drop = FALSE]
      
      # Store OOB data
      tree_info$OOB_Indices <- oob_indices
      tree_info$OOB_XData <- oob_xdata
      tree_info$OOB_YData <- oob_ydata
      
      # Make predictions on OOB data
      oob_predictions <- SCA_tree_predict(
        model = tree_model,
        Testing_data = oob_xdata
      )
      
      # Calculate R-squared for OOB predictions
      if (ncol(oob_ydata) == 1) {
        oob_ydata_numeric <- as.numeric(oob_ydata[[1]])
        oob_predictions_numeric <- as.numeric(oob_predictions[[1]])
        oob_r2 <- 1 - sum((oob_ydata_numeric - oob_predictions_numeric)^2) / 
                  sum((oob_ydata_numeric - mean(oob_ydata_numeric))^2)
      } else {
        oob_r2 <- mean(sapply(1:ncol(oob_ydata), function(i) {
          oob_ydata_numeric <- as.numeric(oob_ydata[,i])
          oob_predictions_numeric <- as.numeric(oob_predictions[,i])
          1 - sum((oob_ydata_numeric - oob_predictions_numeric)^2) / 
              sum((oob_ydata_numeric - mean(oob_ydata_numeric))^2)
        }))
      }
      
      # Add OOB error to model
      tree_model$OOB_error <- oob_r2
      tree_model$OOB_sim <- oob_predictions
      tree_model$Sample <- rep$sample
      tree_model$Tree_Info <- tree_info
      tree_model$Training_data <- tree_data  # Add training data to output
      return(tree_model)
    })
  }
  
  # Calculate weights based on OOB_RSQ
  OOB_RSQ <- sapply(SCE_res, function(x) {
    if (is.null(x$OOB_error)) return(0)
    if (is.na(x$OOB_error)) return(0)
    x$OOB_error
  })
  
  # Handle negative R-squared values by setting them to zero
  OOB_RSQ[OOB_RSQ < 0] <- 0
  
  # Calculate weights
  # First, handle zero values separately
  zero_indices <- which(OOB_RSQ == 0)
  non_zero_indices <- which(OOB_RSQ > 0)
  
  # Calculate weights for non-zero R-squared values
  if (length(non_zero_indices) > 0) {
    weight_OOB <- rep(0, length(OOB_RSQ))
    
    # Prevent OOB_RSQ from being exactly 1 to avoid division by zero
    OOB_RSQ[OOB_RSQ == 1] <- 0.999999
    
    # Calculate log odds ratio safely
    log_odds <- log10(OOB_RSQ[non_zero_indices] / (1 - OOB_RSQ[non_zero_indices]))
    
    # Handle different cases for non-zero weights
    if (length(non_zero_indices) == 1) {
      # If only one non-zero value, give it full weight
      weight_OOB[non_zero_indices] <- 1
    } else if (max(log_odds) == min(log_odds)) {
      # If all non-zero weights are equal, distribute weight equally
      weight_OOB[non_zero_indices] <- 1/length(non_zero_indices)
    } else {
      # Normalize weights safely
      weight_range <- max(log_odds) - min(log_odds)
      if (weight_range == 0) {
        weight_OOB[non_zero_indices] <- 1/length(non_zero_indices)
      } else {
        weight_OOB[non_zero_indices] <- (log_odds - min(log_odds)) / weight_range
        weight_sum <- sum(weight_OOB[non_zero_indices])
        if (weight_sum == 0) {
          weight_OOB[non_zero_indices] <- 1/length(non_zero_indices)
        } else {
          weight_OOB[non_zero_indices] <- weight_OOB[non_zero_indices] / weight_sum
        }
      }
    }
  } else {
    # If all R-squared values are zero, use uniform weights
    weight_OOB <- rep(1/length(OOB_RSQ), length(OOB_RSQ))
  }
  
  SCE_res <- Map(function(x, w) c(x, list(weight = w)), x = SCE_res, w = weight_OOB)
  
  # Create parameters list
  parameters <- list(
    n_trees = Ntree,
    mfeature = mfeature,
    Nmin = Nmin,
    alpha = alpha,
    resolution = resolution,
    n_samples = n_samples,
    n_predictors = length(X),
    n_predictants = length(Y)
  )
  
  # Return S3 class object
  return(SCE_object(
    trees = SCE_res,
    predictors = X,
    predictants = Y,
    parameters = parameters,
    call = call
  ))
}

# S3 Methods for SCE class

#' Generic importance function for S3 method dispatch
#' @param object The object to calculate importance for
#' @param ... Additional arguments passed to methods
#' @export
importance <- function(object, ...) {
  UseMethod("importance")
}

#' Generic evaluate function for S3 method dispatch
#' @param object The object to evaluate
#' @param ... Additional arguments passed to methods
#' @export
evaluate <- function(object, ...) {
  UseMethod("evaluate")
}



#' Print method for SCE objects
#' @param x An SCE object
#' @param ... Additional arguments (not used)
#' @export
print.SCE <- function(x, ...) {
  cat("Stepwise Clustered Ensemble (SCE) Model\n")
  cat("=======================================\n\n")
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Model Parameters:\n")
  cat("  Number of trees:", x$parameters$n_trees, "\n")
  cat("  Features per tree:", x$parameters$mfeature, "\n")
  cat("  Minimum samples per node:", x$parameters$Nmin, "\n")
  cat("  Alpha (significance level):", x$parameters$alpha, "\n")
  cat("  Resolution:", x$parameters$resolution, "\n")
  cat("  Training samples:", x$parameters$n_samples, "\n")
  cat("  Predictors:", x$parameters$n_predictors, "\n")
  cat("  Predictants:", x$parameters$n_predictants, "\n\n")
  
  cat("Predictors:", paste(x$predictors, collapse = ", "), "\n")
  cat("Predictants:", paste(x$predictants, collapse = ", "), "\n\n")
  
  # Calculate average OOB R-squared
  oob_r2 <- sapply(x$trees, function(tree) {
    if (is.null(tree$OOB_error)) return(0)
    if (is.na(tree$OOB_error)) return(0)
    tree$OOB_error
  })
  
  cat("Model Performance:\n")
  cat("  Average OOB R-squared:", round(mean(oob_r2), 4), "\n")
  cat("  Min OOB R-squared:", round(min(oob_r2), 4), "\n")
  cat("  Max OOB R-squared:", round(max(oob_r2), 4), "\n")
  
  invisible(x)
}

#' Summary method for SCE objects
#' @param object An SCE object
#' @param ... Additional arguments (not used)
#' @export
summary.SCE <- function(object, ...) {
  cat("Stepwise Clustered Ensemble (SCE) Model Summary\n")
  cat("==============================================\n\n")
  
  # Basic model info
  cat("Model Information:\n")
  cat("  Number of trees:", object$parameters$n_trees, "\n")
  cat("  Features per tree:", object$parameters$mfeature, "\n")
  cat("  Training samples:", object$parameters$n_samples, "\n")
  cat("  Predictors:", object$parameters$n_predictors, "\n")
  cat("  Predictants:", object$parameters$n_predictants, "\n\n")
  
  # OOB performance statistics
  oob_r2 <- sapply(object$trees, function(tree) {
    if (is.null(tree$OOB_error)) return(0)
    if (is.na(tree$OOB_error)) return(0)
    tree$OOB_error
  })
  
  cat("Out-of-Bag Performance:\n")
  cat("  Mean R-squared:", round(mean(oob_r2), 4), "\n")
  cat("  Median R-squared:", round(median(oob_r2), 4), "\n")
  cat("  Standard deviation:", round(sd(oob_r2), 4), "\n")
  cat("  Min R-squared:", round(min(oob_r2), 4), "\n")
  cat("  Max R-squared:", round(max(oob_r2), 4), "\n")
  cat("  Trees with R-squared > 0.5:", sum(oob_r2 > 0.5), "\n")
  cat("  Trees with R-squared > 0.7:", sum(oob_r2 > 0.7), "\n\n")
  
  # Tree structure statistics
  total_nodes <- sapply(object$trees, function(tree) tree$totalNodes)
  leaf_nodes <- sapply(object$trees, function(tree) tree$leafNodes)
  
  cat("Tree Structure:\n")
  cat("  Average total nodes per tree:", round(mean(total_nodes), 1), "\n")
  cat("  Average leaf nodes per tree:", round(mean(leaf_nodes), 1), "\n")
  cat("  Average tree depth:", round(mean(log2(total_nodes)), 1), "\n\n")
  
  # Weight distribution
  weights <- sapply(object$trees, function(tree) tree$weight)
  cat("Tree Weights:\n")
  cat("  Mean weight:", round(mean(weights), 4), "\n")
  cat("  Weight range:", round(range(weights), 4), "\n")
  cat("  Weight standard deviation:", round(sd(weights), 4), "\n")
  
  invisible(object)
}

#' Predict method for SCE objects
#' @param object An SCE object
#' @param newdata New data for prediction
#' @param ... Additional arguments (not used)
#' @export
predict.SCE <- function(object, newdata, ...) {
  # This is a wrapper for Model_simulation
  if (missing(newdata)) {
    stop("newdata is required for prediction")
  }
  
  # Call Model_simulation which returns Training, Validation, and Testing predictions
  return(Model_simulation(model = object, Testing_data = newdata))
}

#' Importance method for SCE objects
#' @param object An SCE object
#' @param OOB_weight Whether to use OOB weights
#' @param ... Additional arguments (not used)
#' @export
importance.SCE <- function(object, OOB_weight = TRUE, ...) {
  # This is a wrapper for Wilks_importance
  return(Wilks_importance(model = object, OOB_weight = OOB_weight))
}

#' Evaluate method for SCE objects
#' @param object An SCE object
#' @param Testing_data Testing dataset
#' @param Training_data Training dataset
#' @param Predictant Name of predictant variable
#' @param digits Number of digits for output
#' @param ... Additional arguments (not used)
#' @export
evaluate.SCE <- function(object, Testing_data, Training_data, Predictant, digits = 3, ...) {
  # This is a wrapper for SCE_Model_evaluation
  if (missing(Testing_data)) {
    stop("Testing_data is required for evaluation")
  }
  if (missing(Training_data)) {
    stop("Training_data is required for evaluation")
  }
  if (missing(Predictant)) {
    stop("Predictant is required for evaluation")
  }
  
  # Get simulations using Model_simulation
  Simulations <- Model_simulation(model = object, Testing_data = Testing_data)
  
  # Extract validation simulations (which should match Training_data rows)
  Validation_simulations <- Simulations$Validation
  
  # Call SCE_Model_evaluation with validation simulations
  return(SCE_Model_evaluation(
    Testing_data = Testing_data,
    Training_data = Training_data,
    Simulations = Validation_simulations,
    Predictant = Predictant,
    digits = digits
  ))
}

