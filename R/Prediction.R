#################################################################
# Filename: 	Prediction.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Author: 		Kailong Li & Xiuquan Wang
# Email:		  lkl98509509@gmail.com; xiuquan.wang@gmail.com
# ===============================================================
# History: 	2019/06/06		modified the the rsca_predicion.r file from the rsca package developed by Xiuquan Wang
#			      2019/06/08		revised model inference function, by Kailong Li
##################################################################

# ---------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------

#: predict one input
#: input: o_precitor = c(x1, x2, x3, ... )
#: output: o_predictant = c(y1, y2, y3, ... )
f_predict_one <- function(data, o_precitor) {
  # Initialize with NA instead of 0 to detect if prediction fails
  o_predictant <- rep(NA, data$n_y_cols)

  n_tree_index <- 1
  while(n_tree_index <= data$n_result_tree_rows) {
    # Check if current node is a leaf
    if (data$o_result_tree[n_tree_index, 4] == -1 && data$o_result_tree[n_tree_index, 5] == -1) {
      # Get the mean y value from the leaf node
      node_id <- data$o_result_tree[n_tree_index, 1]
      if (node_id > 0 && node_id <= nrow(data$o_mean_y)) {
        o_predictant <- as.numeric(data$o_mean_y[node_id, ])
      } else {
        warning("Invalid node ID encountered: ", node_id)
      }
      break
    }

    # Get the column index for comparison
    n_j <- data$o_result_tree[n_tree_index, 2]
    
    # Handle merged nodes
    if (n_j <= 0) {
      n_tree_index <- data$o_result_tree[n_tree_index, 4]
      next
    }
    
    # Compare predictor value with node's split value
    if (o_precitor[n_j] <= data$o_result_tree[n_tree_index, 3]) {
      n_tree_index <- data$o_result_tree[n_tree_index, 4]
    } else {
      n_tree_index <- data$o_result_tree[n_tree_index, 5]
    }
  }

  # Check if prediction was successful
  if (any(is.na(o_predictant))) {
    warning("Failed to make prediction for input: ", paste(o_precitor, collapse=", "))
  }

  return(o_predictant)
}

#: do prediction to all input
f_predict <- function(data) {
  for (iPredictor in 1:data$n_sample_size) {
    o_vector_predictor <- as.numeric(data$o_sample_data_x[iPredictor, ])
    o_vector_predictant <- f_predict_one(data, o_vector_predictor)
    data$o_predictants[iPredictor, ] <- o_vector_predictant
  }
  return(data)
}

# ---------------------------------------------------------------
# Main function
# ---------------------------------------------------------------
f_main_p <- function(data) {
  # Verify input data
  if (is.null(data$o_result_tree) || nrow(data$o_result_tree) == 0) {
    stop("Model tree is empty or not properly initialized")
  }
  
  if (is.null(data$o_mean_y) || nrow(data$o_mean_y) == 0) {
    stop("Model map is empty or not properly initialized")
  }
  
  # Perform prediction
  data <- f_predict(data)
  
  # Format output based on model type
  if (data$n_model_type == "mean" || data$n_model_type == "max" || 
      data$n_model_type == "min" || data$n_model_type == "median") {
    colnames(data$o_predictants) <- colnames(data$o_mean_y)
  } else if (data$n_model_type == "interval") {
    # Format for interval type
    n_midcol <- data$n_y_cols / 2
    colnames(data$o_predictants) <- c(
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_min"),
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_max")
    )
  } else if (data$n_model_type == "radius") {
    # Format for radius type
    n_midcol <- data$n_y_cols / 2
    colnames(data$o_predictants) <- c(
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_mean"),
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_radius")
    )
  } else if (data$n_model_type == "variation") {
    # Format for variation type
    n_midcol <- data$n_y_cols / 2
    colnames(data$o_predictants) <- c(
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_mean"),
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_sd")
    )
  } else if (data$n_model_type == "random") {
    # Format for random type
    n_midcol <- data$n_y_cols / 2
    colnames(data$o_predictants) <- c(
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_min"),
      paste0(colnames(data$o_mean_y)[1:n_midcol], "_max")
    )
  }
  
  return(data)
}

# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
Inference <- function(x, Weak_L) {
  # Initialize data structure
  data <- list()
  
  # Store input data
  data$o_sample_data_x <- as.matrix(na.omit(x))
  data$n_sample_size <- nrow(data$o_sample_data_x)
  data$n_sample_x_cols <- ncol(data$o_sample_data_x)
  
  # Store model data
  data$o_result_tree <- Weak_L$Tree
  data$n_result_tree_rows <- nrow(data$o_result_tree)
  data$o_mean_y <- Weak_L$Map
  data$n_y_cols <- ncol(data$o_mean_y)
  data$n_model_type <- Weak_L$type
  
  # Initialize prediction matrix
  data$o_predictants <- matrix(NA, data$n_sample_size, data$n_y_cols)
  
  # Perform prediction
  data <- f_main_p(data)
  
  return(data$o_predictants)
}
