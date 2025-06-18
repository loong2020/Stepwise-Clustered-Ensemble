#################################################################
# Filename: 	Training.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Author: 		Kailong Li & Xiuquan Wang
# Email:		  lkl98509509@gmail.com; xiuquan.wang@gmail.com
# ===============================================================
# History: 	2019/06/15		modified the the rsca_training.r file from the rsca package developed by Xiuquan Wang
#			      2019/06/15		fixed the issue of infinite cut-merge loop, by Kailong Li
#			      2019/06/18		changed model input and output from .csv/.txt to data.frame, by Kailong Li
#			      2019/09/15		added minimal Wilks value to the tree output, by Kailong Li
#			      2019/06/08		revised model inference function, by Kailong Li
#			      2025/04/12		added resolution search for the best split point, by Kailong Li
##################################################################

#: calculate wilks statistic value
#: input: top matrix & bot matrix
#: output: wilks value
f_wilks_statistic <- function(o_top_matrix, o_bot_matrix, Nmin)
{
  n_top = nrow(o_top_matrix)
  n_bot = nrow(o_bot_matrix)
  if ((n_top + n_bot) <= (ncol(o_top_matrix) + Nmin))
  {
    return(0)
  }

  o_top_mean = matrix(colMeans(o_top_matrix),nrow=1)
  o_bot_mean = matrix(colMeans(o_bot_matrix),nrow=1)
  o_between_matrix = (n_top*n_bot)/(n_top+n_bot) * crossprod(o_top_mean - o_bot_mean, o_top_mean - o_bot_mean)

  o_top_matrix_mean = matrix(colMeans(o_top_matrix), nrow(o_top_matrix), ncol(o_top_matrix), byrow = TRUE)
  o_bot_matrix_mean = matrix(colMeans(o_bot_matrix), nrow(o_bot_matrix), ncol(o_bot_matrix), byrow = TRUE)
  o_within_top = crossprod(o_top_matrix - o_top_matrix_mean, o_top_matrix - o_top_matrix_mean)
  o_within_bot = crossprod(o_bot_matrix - o_bot_matrix_mean, o_bot_matrix - o_bot_matrix_mean)
  o_within_matrix = o_within_top + o_within_bot

  o_wilks_vaule = 0
  n_det_within = det(o_within_matrix)
  n_det_total = det(o_within_matrix + o_between_matrix)

  # if the determinant of the total matrix is 0, then the wilks value is 0 (no difference between the two groups)
  if (n_det_total == 0)
  {
    # if the determinant of the within matrix is negative, which suggests an invalid or degenerate case in the data.
    if (n_det_within < 0)
      o_wilks_vaule = -1
    # This indicates that the determinant of the within-group matrix is positive, which suggests a valid but homogeneous case where the groups are identical.
    else if (n_det_within > 0)
      o_wilks_vaule = 1
    # This indicates that both the total matrix determinant (n_det_total) and the within-group matrix determinant (n_det_within) are 0, meaning there is no difference between the groups.
    else
      o_wilks_vaule = 0
  }
  else
    o_wilks_vaule = n_det_within / n_det_total

  return(o_wilks_vaule)
}

#: Helper function for iterative refinement
find_best_split_iterative <- function(sorted_y, sorted_x, Nmin, resolution) {
  n <- nrow(sorted_y)
  p <- ncol(sorted_y)
  
  # cat(sprintf("Starting find_best_split_iterative with n=%d, p=%d, resolution=%d\n", n, p, resolution))
  
  # Initialize with all possible split points
  current_indices <- 1:(n-1)
  best_split <- n/2
  best_wilks <- 1
  iteration <- 0
  
  # Early stopping parameters
  patience <- 20
  improvement_threshold <- 0.01  # 1% improvement
  no_improvement_count <- 0
  previous_mean_wilks <- 1
  
  while (length(current_indices) > resolution & length(current_indices) > Nmin) {
    iteration <- iteration + 1
    
    # cat(sprintf("Iteration %d: Current number of indices: %d\n", iteration, length(current_indices)))
    
    # Take resolution evenly spaced points from current indices
    selected_indices <- round(seq(1, length(current_indices), length.out = resolution))
    split_points <- current_indices[selected_indices]
    
    # cat(sprintf("  Number of split points to evaluate: %d\n", length(split_points)))
    
    # Calculate Wilks' lambda for each split point using f_wilks_statistic
    wilks_values <- sapply(split_points, function(i) {
      f_wilks_statistic(
        o_top_matrix = sorted_y[1:i, , drop = FALSE],
        o_bot_matrix = sorted_y[(i+1):n, , drop = FALSE],
        Nmin = Nmin
      )
    })
    
    # cat("  Wilks values range: ", range(wilks_values), "\n")
    
    # Find promising points using only valid Wilks' values
    valid_indices <- which(wilks_values > 0 & wilks_values < 1)   
    valid_wilks <- wilks_values[valid_indices]
    
    # cat(sprintf("  Number of valid splits: %d\n", length(valid_indices)))
    
    # If no valid splits found, check if all values are 0 or 1
    if (length(valid_indices) == 0) {
      # No valid splits found (all values are either 0, 1, or -1)
      return(list(wilks = 0, split = current_indices[length(current_indices) %/% 2]))
    }
    
    mean_wilks <- mean(valid_wilks)
    sd_wilks <- sd(valid_wilks)
    promising_points <- valid_indices[which(valid_wilks < (mean_wilks - 0.5 * sd_wilks))]
    
    # cat(sprintf("  Number of promising points: %d\n", length(promising_points)))
    
    if (length(promising_points) == 0) {
      # cat("  No promising points found, breaking loop\n")
      break
    }
    
    # Check for early stopping
    if (iteration > 1) {
      improvement <- (previous_mean_wilks - mean_wilks) / previous_mean_wilks
      if (improvement < improvement_threshold) {
        no_improvement_count <- no_improvement_count + 1
        # cat(sprintf("  No significant improvement (%.2f%% < %.2f%%), patience count: %d/%d\n", 
        #            improvement * 100, improvement_threshold * 100, no_improvement_count, patience))
      } else {
        no_improvement_count <- 0
        # cat(sprintf("  Improvement: %.2f%%\n", improvement * 100))
      }
      
      if (no_improvement_count >= patience) {
        # cat("  Early stopping: No significant improvement for", patience, "iterations\n")
        break
      }
    }
    previous_mean_wilks <- mean_wilks
    
    # Create new indices from promising regions
    new_indices <- integer(0)
    
    # Sort promising points to handle regions in order
    promising_points <- sort(promising_points)
    
    for (i in seq_along(promising_points)) {
      point_idx <- promising_points[i]
      
      # Determine region boundaries
      left_boundary <- ifelse(i == 1, 
                             max(1, promising_points[i] - 1),
                             promising_points[i-1])
      right_boundary <- ifelse(i == length(promising_points), 
                              min(length(split_points), promising_points[i] + 1),
                              promising_points[i+1])
      
      # Get the actual indices in the original space
      region_start <- split_points[left_boundary] + 1
      region_end <- split_points[right_boundary] - 1
      
      # Ensure we don't go beyond data boundaries
      region_start <- max(1, region_start)
      region_end <- min(n-1, region_end)
      
      # Add all indices from this region
      region_indices <- region_start:region_end
      new_indices <- c(new_indices, region_indices)
      
      # cat(sprintf("  Region %d: point %d (Wilks=%.4f), boundaries [%d, %d], points [%d, %d]\n",
      #            i, split_points[point_idx], wilks_values[point_idx],
      #            region_start, region_end, min(region_indices), max(region_indices)))
    }
    
    # Remove duplicates and sort
    new_indices <- sort(unique(new_indices))
    
    # cat(sprintf("  New number of indices: %d\n", length(new_indices)))
    
    # Update best split if found
    current_best <- which.min(valid_wilks)
    if (valid_wilks[current_best] < best_wilks) {
      best_wilks <- valid_wilks[current_best]
      best_split <- split_points[current_best]
      # cat(sprintf("  New best split found: point %d with Wilks=%.4f\n", best_split, best_wilks))
    }
    
    # Check if we're making progress
    if (length(new_indices) >= length(current_indices)) {
      # cat("  No progress in reducing indices, breaking loop\n")
      break
    }
    
    current_indices <- new_indices
  }
  
  # cat(sprintf("Final number of indices: %d\n", length(current_indices)))
  
  # Final detailed search on remaining points
  final_wilks <- sapply(current_indices, function(i) {
    f_wilks_statistic(
      o_top_matrix = sorted_y[1:i, , drop = FALSE],
      o_bot_matrix = sorted_y[(i+1):n, , drop = FALSE],
      Nmin = Nmin
    )
  })
  
  # Filter out invalid Wilks' values in final search
  valid_final_indices <- which(final_wilks > 0 & final_wilks < 1)
  if (length(valid_final_indices) > 0) {
    final_best <- valid_final_indices[which.min(final_wilks[valid_final_indices])]
    if (final_wilks[final_best] < best_wilks) {
      best_wilks <- final_wilks[final_best]
      best_split <- current_indices[final_best]
      # cat(sprintf("  Final search found better split: point %d with Wilks=%.4f\n", best_split, best_wilks))
    }
  } else {
    # cat("  No valid splits found in final search, returning (wilks = 0, split = middle point)\n")
    return(list(wilks = 0, split = current_indices[length(current_indices) %/% 2]))
  }
  
  # cat(sprintf("Best split found at point %d with Wilks' lambda = %.4f\n", best_split, best_wilks))
  
  return(list(wilks = best_wilks, split = best_split))
}

find_best_split <- function(sorted_y, sorted_x, start, end, Nmin) {
  n <- nrow(sorted_y)
  
  # Calculate all Wilks' values first
  split_points <- start:end
  wilks_values <- sapply(split_points, function(i) {
    f_wilks_statistic(
      o_top_matrix = sorted_y[1:i, , drop = FALSE],
      o_bot_matrix = sorted_y[(i+1):n, , drop = FALSE],
      Nmin = Nmin
    )
  })
  
  # Find the best valid split (0 < wilks < 1)
  valid_indices <- which(wilks_values > 0 & wilks_values < 1)
  if (length(valid_indices) > 0) {
    best_idx <- valid_indices[which.min(wilks_values[valid_indices])]
    return(list(wilks = wilks_values[best_idx], split = split_points[best_idx]))
  }
  
  # If no valid splits, check for splits with wilks = 0, indicating the data is too small or homogeneous to be split
  zero_indices <- which(wilks_values == 0)
  if (length(zero_indices) > 0) {
    return(list(wilks = 0, split = split_points[length(split_points) %/% 2]))
  }
  
  # If all splits return 1, use middle point
  return(list(wilks = 1, split = split_points[length(split_points) %/% 2]))
}

#: calculate minimum wilks value for a matrix indentified by a group of row id [original]
#: input: matrix of rowid [nrows, col=1]
#: output: min_wilks_list(min_wilks_value=1, col_id=1, x_value=1, left_rowids=matrix(), right_rowids=matrix())
f_min_wilks <- function(data, o_matrix_rowid, resolution)
{
  if (nrow(o_matrix_rowid) <= (data$n_sample_y_cols + data$Nmin))
  {
    return(list(min_wilks_value=0, col_id=0, x_value=0, left_rowids=matrix(), right_rowids=matrix()))
  }

  # Pre-calculate y matrix once and ensure it's a matrix
  y_matrix <- as.matrix(data$o_sample_data_y[o_matrix_rowid, , drop = FALSE])
  n <- nrow(y_matrix)
  
  # Initialize results
  results <- matrix(Inf, data$n_sample_x_cols, 4)
  
  for (icol in 1:data$n_sample_x_cols) {
    # Sort data for current column
    sorted_indices <- order(data$o_sample_data_x[o_matrix_rowid, icol])
    sorted_y <- y_matrix[sorted_indices, , drop = FALSE]
    sorted_x <- data$o_sample_data_x[o_matrix_rowid[sorted_indices], icol]
    
    if (n > resolution) {
      # Use iterative refinement
      best_result <- find_best_split_iterative(sorted_y, sorted_x, data$Nmin, resolution)
      best_wilks <- best_result$wilks
      best_split <- best_result$split
    } else {
      # Full search for small datasets
      result <- find_best_split(sorted_y, sorted_x, 1, n-1, data$Nmin)
      best_wilks <- result$wilks
      best_split <- result$split
    }
    
    # Store results - handle case when no valid split is found
    results[icol, ] <- c(best_wilks, icol, sorted_x[best_split], best_split)
  }
  
  # Find best split
  best_col <- which.min(results[, 1])
  best_split <- results[best_col, ]
  
  # Construct return value
  sorted_indices <- order(data$o_sample_data_x[o_matrix_rowid, best_col])
  return(list(
    min_wilks_value = best_split[1],
    col_id = best_split[2],
    x_value = best_split[3],
    left_rowids = matrix(o_matrix_rowid[sorted_indices[1:best_split[4]]], ncol=1),
    right_rowids = matrix(o_matrix_rowid[sorted_indices[(best_split[4]+1):nrow(o_matrix_rowid)]], ncol=1)
  ))
}

#: calculate f statistic value and check if it can be divided or cut
#: input: min_wilks_list
#: output: 0 -> no, 1 -> yes
f_cal_chk_f <- function(data, min_wilks_list)
{
  n_row_left = nrow(min_wilks_list$left_rowids)
  n_row_right = nrow(min_wilks_list$right_rowids)
  if ((n_row_left + n_row_right) <= (data$n_sample_y_cols + data$Nmin))
  {
    return(0)
  }
  if (is.na(as.numeric(min_wilks_list$min_wilks_value)))
  {
    return(1)
  }
  if (as.numeric(min_wilks_list$min_wilks_value) == 0)
  {
    return(0)  # Changed from 1 to 0 - no cut when groups are identical
  }
  if (as.numeric(min_wilks_list$min_wilks_value) < 0 || as.numeric(min_wilks_list$min_wilks_value) > 1)
  {
    return(1)
  }
  o_df_numerator = data$n_sample_y_cols
  o_df_dominator = n_row_left + n_row_right - data$n_sample_y_cols - 1

  o_f_value = ((1 - min_wilks_list$min_wilks_value) / min_wilks_list$min_wilks_value) * (o_df_dominator / o_df_numerator)
  o_f_criterion = qf((1-data$n_alpha), o_df_numerator, o_df_dominator)

  o_check_flag = 0
  if (o_f_value >= o_f_criterion) o_check_flag = 1
  return(o_check_flag)
}

#: check if node is a leaf
#: input: node id
#: output: 1:yes, 0:no
f_checkif_leaf <- function(data, n_nodeid)
{
  # Check if node exists
  if (is.null(data$o_output_tree[[n_nodeid]])) {
    data$o_output_tree[[n_nodeid]]$left <- -1
    data$o_output_tree[[n_nodeid]]$right <- -1
    data$o_output_tree[[n_nodeid]]$left_mat <- -1
    data$o_output_tree[[n_nodeid]]$right_mat <- -1
    data$o_output_tree[[n_nodeid]]$wilk_min <- 0
    return(list(data = data, is_leaf = 1))
  }
  
  # Check if rowids_matrix exists and is a matrix
  if (is.null(data$o_output_tree[[n_nodeid]]$rowids_matrix)) {
    data$o_output_tree[[n_nodeid]]$left <- -1
    data$o_output_tree[[n_nodeid]]$right <- -1
    data$o_output_tree[[n_nodeid]]$left_mat <- -1
    data$o_output_tree[[n_nodeid]]$right_mat <- -1
    data$o_output_tree[[n_nodeid]]$wilk_min <- 0
    return(list(data = data, is_leaf = 1))
  }
  
  if (!is.matrix(data$o_output_tree[[n_nodeid]]$rowids_matrix)) {
    data$o_output_tree[[n_nodeid]]$left <- -1
    data$o_output_tree[[n_nodeid]]$right <- -1
    data$o_output_tree[[n_nodeid]]$left_mat <- -1
    data$o_output_tree[[n_nodeid]]$right_mat <- -1
    data$o_output_tree[[n_nodeid]]$wilk_min <- 0
    return(list(data = data, is_leaf = 1))
  }

  if (nrow(data$o_output_tree[[n_nodeid]]$rowids_matrix) <= (data$n_sample_y_cols + data$Nmin))
  {
    data$o_output_tree[[n_nodeid]]$left <- -1
    data$o_output_tree[[n_nodeid]]$right <- -1
    data$o_output_tree[[n_nodeid]]$left_mat <- -1
    data$o_output_tree[[n_nodeid]]$right_mat <- -1
    data$o_output_tree[[n_nodeid]]$wilk_min <- 0
    return(list(data = data, is_leaf = 1))
  }
  
  if (data$o_output_tree[[n_nodeid]]$left == -1 && data$o_output_tree[[n_nodeid]]$right == -1)
  {
    return(list(data = data, is_leaf = 1))
  }
  
  return(list(data = data, is_leaf = 0))
}

#: process node
f_processnode <- function(data, node_id, Max_merge_iter, resolution, verbose = FALSE)
{
  #: if a leaf, just exit this function
  leaf_check <- f_checkif_leaf(data, node_id)
  data <- leaf_check$data
  if (leaf_check$is_leaf == 1)
  {
    return(data)
  }
  
  #: add this node into stack
  data$o_nodeid_stack_cut[data$n_nodeid_statck_cut_cursor] <- node_id
  data$n_nodeid_statck_cut_cursor <- data$n_nodeid_statck_cut_cursor + 1

  #: loop counter
  n_loop_counter = 0

  #: initiate the loop
  data$n_flag_cut <- 1
  while(data$n_flag_cut == 1)
  {
    n_loop_counter = n_loop_counter + 1

    #: do while only if cutting occured
    #: clear cut flag
    data$n_flag_cut <- 0

    #: do cut till the cut stack is empty
    while(data$n_nodeid_statck_cut_cursor > 1)
    {
      #: read a node from stack
      n_nodeid_cut_temp = data$o_nodeid_stack_cut[data$n_nodeid_statck_cut_cursor - 1]

      #: if this node has been cutted and set as a leaf, continue the next node
      leaf_check <- f_checkif_leaf(data, n_nodeid_cut_temp)
      data <- leaf_check$data
      if (leaf_check$is_leaf >= 1)
      {
        data$o_nodeid_stack_merge[data$n_nodeid_statck_merge_cursor] <- n_nodeid_cut_temp
        data$n_nodeid_statck_merge_cursor <- data$n_nodeid_statck_merge_cursor + 1

        data$n_nodeid_statck_cut_cursor <- data$n_nodeid_statck_cut_cursor - 1
        next
      }

      #: calculate minimum wilks value
      min_wilks_list = f_min_wilks(data, data$o_output_tree[[n_nodeid_cut_temp]]$rowids_matrix, resolution)

      #: check if this node can be cut
      n_cut_flag = f_cal_chk_f(data, min_wilks_list)
      
      # Print cutting information only if verbose is enabled
      if (verbose) {
        message(sprintf("Cutting Node %d: Wilks' lambda = %.4f, Cut Flag = %d", 
                  n_nodeid_cut_temp, min_wilks_list$min_wilks_value, n_cut_flag))
      }
      
      if ((n_cut_flag == 1) && (data$Merge_iter <= Max_merge_iter))
      {
        #: if can be cut, then cut it and process it's sub nodes, respectively
        n_cursor_tree = length(data$o_output_tree) + 1
        
        #: update the current node
        data$o_output_tree[[n_nodeid_cut_temp]]$col_index <- min_wilks_list$col_id
        data$o_output_tree[[n_nodeid_cut_temp]]$value <- min_wilks_list$x_value
        data$o_output_tree[[n_nodeid_cut_temp]]$left <- n_cursor_tree
        data$o_output_tree[[n_nodeid_cut_temp]]$right <- n_cursor_tree + 1
        data$o_output_tree[[n_nodeid_cut_temp]]$wilk_min <- min_wilks_list$min_wilks_value
        data$o_output_tree[[n_nodeid_cut_temp]]$lfMat <- nrow(min_wilks_list$left_rowids)
        data$o_output_tree[[n_nodeid_cut_temp]]$rtMat <- nrow(min_wilks_list$right_rowids)
        
        #: create left node
        o_left_tree_list <- list(
          id = n_cursor_tree,
          col_index = 0,
          value = 0,
          left = 0,
          right = 0,
          rowids_matrix = min_wilks_list$left_rowids,
          left_mat = 0,
          right_mat = 0,
          wilk_min = 0,
          compared = numeric(0)
        )
        data$o_output_tree[[n_cursor_tree]] <- o_left_tree_list
        
        #: create right node
        o_right_tree_list <- list(
          id = n_cursor_tree + 1,
          col_index = 0,
          value = 0,
          left = 0,
          right = 0,
          rowids_matrix = min_wilks_list$right_rowids,
          left_mat = 0,
          right_mat = 0,
          wilk_min = 0,
          compared = numeric(0)
        )
        data$o_output_tree[[n_cursor_tree + 1]] <- o_right_tree_list
        
        data$n_cut_times <- data$n_cut_times + 1

        #: delete this node (replace) and store the 2 sub nodes into cut stack
        data$o_nodeid_stack_cut[data$n_nodeid_statck_cut_cursor - 1] <- n_cursor_tree
        data$o_nodeid_stack_cut[data$n_nodeid_statck_cut_cursor] <- n_cursor_tree + 1
        data$n_nodeid_statck_cut_cursor <- data$n_nodeid_statck_cut_cursor + 1

        #: update cut flag
        data$n_flag_cut <- 1
      }
      else if ((n_cut_flag == 0) && (data$Merge_iter <= Max_merge_iter))
      {
        #: is a leaf, set it as a leaf and add into merge stack
        data$o_output_tree[[n_nodeid_cut_temp]]$left <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$right <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$left_mat <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$right_mat <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$wilk_min <- 0
        data$o_nodeid_stack_merge[data$n_nodeid_statck_merge_cursor] <- n_nodeid_cut_temp
        data$n_nodeid_statck_merge_cursor <- data$n_nodeid_statck_merge_cursor + 1
        #: delete this node from cut stack
        data$n_nodeid_statck_cut_cursor <- data$n_nodeid_statck_cut_cursor - 1
      }
      else
      {
        data$o_output_tree[[n_nodeid_cut_temp]]$left <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$right <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$left_mat <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$right_mat <- -1
        data$o_output_tree[[n_nodeid_cut_temp]]$wilk_min <- 0
        #: update cut flag
        data$n_flag_cut <- 0
      }
    }

    data$o_nodeid_stack_cut <- c()
    data$n_nodeid_statck_cut_cursor <- 1

    #: clear merge flag
    data$n_flag_merge <- 0

    #: do merge if the merge stack is not empty or until merge reach its maximum times
    if(data$n_nodeid_statck_merge_cursor > 1)
    {
      #: initiate the loop
      data$n_flag_merge <- 1
      while(data$n_flag_merge == 1)
      {
        data$n_flag_merge <- 0

        #: bottom index of the merge stack -> used for the no merging case in one for loop
        n_bot_index_merge_stack = 1
        while(data$n_nodeid_statck_merge_cursor > n_bot_index_merge_stack)
        {
          imerge_a = data$n_nodeid_statck_merge_cursor - 1

          #: read a node from stack
          n_nodeid_merge_a = data$o_nodeid_stack_merge[imerge_a]

          #: try to merge with other nodes (including the newly merged nodes)
          imerge_b = imerge_a - 1
          while(imerge_b >= n_bot_index_merge_stack)
          {
            n_nodeid_merge_b = data$o_nodeid_stack_merge[imerge_b]

            #: if the total number of each sample (a, b) is lower than (rSCA.env$n_sample_y_cols + Nmin)
            #: then the wilks value is 0, but now they can not be merged!!!
            if (nrow(data$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix) <= (data$n_sample_y_cols + data$Nmin) || 
                nrow(data$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix) <= (data$n_sample_y_cols + data$Nmin))
            {
              imerge_b = imerge_b - 1
              next
            }

            #: construct the top and bot matrix and calculate wilks value
            o_top_matrix_temp = matrix(0, nrow(data$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix), data$n_sample_y_cols)
            o_top_matrix_temp[1:nrow(o_top_matrix_temp), ] = data.matrix(data$o_sample_data_y[c(data$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix), ])

            o_bot_matrix_temp = matrix(0, nrow(data$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix), data$n_sample_y_cols)
            o_bot_matrix_temp[1:nrow(o_bot_matrix_temp), ] = data.matrix(data$o_sample_data_y[c(data$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix), ])

            n_wilks_value = f_wilks_statistic(o_top_matrix_temp, o_bot_matrix_temp, data$Nmin)

            #: 2> contruct min wilks list
            o_temp_min_wilks_list = list(min_wilks_value=n_wilks_value, col_id=0, x_value=0, 
                                       left_rowids=data$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix, 
                                       right_rowids=data$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix)

            #: 3> check if they can be merged
            n_merge_flag = f_cal_chk_f(data, o_temp_min_wilks_list)

            # Comment out merging information print
            # cat(sprintf("Merging Nodes %d and %d: Wilks' lambda = %.4f, Merge Flag = %d\n", 
            #            n_nodeid_merge_a, n_nodeid_merge_b, n_wilks_value, n_merge_flag))

            if (n_merge_flag == 0)
            {
              #: can be merged
              #: 1> do merge
              n_cursor_tree = length(data$o_output_tree) + 1

              # Print merging information only if verbose is enabled
              if (verbose) {
                message(sprintf("Merging nodes %d and %d into new node %d", 
                         n_nodeid_merge_a, n_nodeid_merge_b, n_cursor_tree))
              }
              
              #: set the cursor for output tree
              o_temp_list = list(id=n_cursor_tree, col_index=0, value=0, left=0, right=0, 
                               rowids_matrix=rbind(o_temp_min_wilks_list$left_rowids, o_temp_min_wilks_list$right_rowids))
              data$o_output_tree[[n_cursor_tree]] <- o_temp_list

              data$o_output_tree[[n_nodeid_merge_a]]$left <- n_cursor_tree
              data$o_output_tree[[n_nodeid_merge_a]]$right <- n_cursor_tree
              data$o_output_tree[[n_nodeid_merge_a]]$left_mat <- -1
              data$o_output_tree[[n_nodeid_merge_a]]$right_mat <- -1
              data$o_output_tree[[n_nodeid_merge_a]]$wilk_min <- 0

              data$o_output_tree[[n_nodeid_merge_b]]$left <- n_cursor_tree
              data$o_output_tree[[n_nodeid_merge_b]]$right <- n_cursor_tree

              #: 2> store the new node into merge stack at [imerge_b]
              data$o_nodeid_stack_merge[imerge_b] <- n_cursor_tree

              #: delete [imerge_a] from merge stack
              data$n_nodeid_statck_merge_cursor <- imerge_a

              #: 3> update merge times variable
              data$n_flag_merge <- 1
              data$n_merge_times <- data$n_merge_times + 1
              #: 4> break the FOR loop
              break
            }
            else
            {
              imerge_b = imerge_b - 1
            }
          }
          if (data$n_flag_merge == 1)
          {
            #: there is merging action occurring, then restart the whole loop
            break
          }
          else
          {
            #: can not be merged with other nodes
            #: exchange imerge_a with the node pointed by n_bot_index_merge_stack
            o_temp_merge = data$o_nodeid_stack_merge[imerge_a]
            data$o_nodeid_stack_merge[imerge_a] <- data$o_nodeid_stack_merge[n_bot_index_merge_stack]
            data$o_nodeid_stack_merge[n_bot_index_merge_stack] <- o_temp_merge
            n_bot_index_merge_stack = n_bot_index_merge_stack + 1
          }
        }
      }

      #: copy all node in the merge stack into the cut stack, continue to cut
      data$o_nodeid_stack_cut <- data$o_nodeid_stack_merge[1:(data$n_nodeid_statck_merge_cursor-1)]
      data$n_nodeid_statck_cut_cursor <- data$n_nodeid_statck_merge_cursor

      data$o_nodeid_stack_merge <- c()
      data$n_nodeid_statck_merge_cursor <- 1
      data$Merge_iter = data$Merge_iter + 1
    }
  }
  
  return(data)
}


# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
do_cluster <- function(data, Nmin, resolution, verbose = FALSE)
{
  #: store the start time
  time_stat <- proc.time()

  #: initialize
  result <- f_init(data)

  #: do main function
  result <- f_main(result, Max_merge_iter=10, Nmin=Nmin, resolution = resolution, verbose = verbose)

  # : calculate the total time used
  time_end <- (proc.time() - time_stat)[[3]]
  Hours <- time_end %/% (60*60)
  Minutes <- (time_end %% 3600) %/% 60
  Seconds <- time_end %% 60
  time_used <- paste(Hours, " h ", Minutes, " m ", Seconds, " s.", sep="")
  
  # Print time information only if verbose is enabled
  if (verbose) {
    message("Time Used: ", time_used)
  }
  
  return(result)
}

# ---------------------------------------------------------------
# Initialization functions
# ---------------------------------------------------------------
f_init <- function(data)
{
  #: matrix to store the sorted results
  data$o_sorted_matrix <- matrix(0, data$n_sample_size, data$n_sample_x_cols)
  data$o_sorted_temp_matrix <- matrix(0, data$n_sample_size, data$n_sample_x_cols)

  #: do sorting
  for (col in 1:data$n_sample_x_cols)
  {
    data$o_sorted_temp_matrix[,col] <- order(data$o_sample_data_x[,col])
  }
  for (col in 1:data$n_sample_x_cols)
  {
    for (row in 1:data$n_sample_size)
    {
      data$o_sorted_matrix[data$o_sorted_temp_matrix[row, col], col] <- row
    }
  }

  #: list for output tree
  data$o_output_tree <- list()

  #: initiate output tree with all required fields
  o_init_tree_list <- list(
    id = 1,
    col_index = 0,
    value = 0,
    left = 0,
    right = 0,
    rowids_matrix = matrix(1:data$n_sample_size, ncol = 1),
    left_mat = 0,
    right_mat = 0,
    wilk_min = 0,
    compared = numeric(0)  # Add compared field for merge operations
  )
  data$o_output_tree[[1]] <- o_init_tree_list

  #: stack to store the unprocessed node ids in the output tree
  data$o_nodeid_stack_cut <- c()
  data$n_nodeid_statck_cut_cursor <- 1
  data$o_nodeid_stack_merge <- c()
  data$n_nodeid_statck_merge_cursor <- 1

  #: cutting and merging flags for loop, 1:do loop, 0:no need
  data$n_flag_cut <- 0
  data$n_flag_merge <- 0

  #: some statistical infos for the results
  data$n_cut_times <- 0
  data$n_merge_times <- 0
  data$Merge_iter <- 0
  data$n_leafnodes_count <- 0

  return(data)
}

# ---------------------------------------------------------------
# Main functions
# ---------------------------------------------------------------
f_main <- function(data, Max_merge_iter, Nmin, resolution, verbose = FALSE)
{
  #: cut from the root node
  data <- f_processnode(data, 1, Max_merge_iter, resolution, verbose)

  #: results matrix structure -> matrix(id, col_id, x_value, left_id, right_id, left_mat, right_mat, wilk_min)
  o_results_matrix = matrix(0, length(data$o_output_tree), 8)

  #: if mapvalue set as "mean" ==> calculate the mean of all samples
  n_mapfile_cols = data$n_sample_y_cols
  if (data$n_mapvalue == "interval" || data$n_mapvalue == "radius" || data$n_mapvalue == "variation")
    n_mapfile_cols = data$n_sample_y_cols * 2

  o_y_results_matrix = matrix(0, length(data$o_output_tree), n_mapfile_cols)

  for(itree in 1:length(data$o_output_tree))
  {
    o_results_matrix[itree,1] = as.numeric(data$o_output_tree[[itree]]$id)
    o_results_matrix[itree,2] = as.numeric(data$o_output_tree[[itree]]$col_index)
    o_results_matrix[itree,3] = as.numeric(data$o_output_tree[[itree]]$value)
    o_results_matrix[itree,4] = as.numeric(data$o_output_tree[[itree]]$left)
    o_results_matrix[itree,5] = as.numeric(data$o_output_tree[[itree]]$right)
    
    # Handle merged nodes (where left == right)
    if (data$o_output_tree[[itree]]$left == data$o_output_tree[[itree]]$right) {
      if (data$o_output_tree[[itree]]$left == -1) {
        # Leaf node
      o_results_matrix[itree,6] = -1
      o_results_matrix[itree,7] = -1
    } else {
        # Merged node - use the size of the combined node
        merged_size = nrow(data$o_output_tree[[data$o_output_tree[[itree]]$left]]$rowids_matrix)
        o_results_matrix[itree,6] = merged_size
        o_results_matrix[itree,7] = merged_size
    }
    } else {
      # Regular split node
      o_results_matrix[itree,6] = ifelse(data$o_output_tree[[itree]]$left == -1, -1, 
                                        as.numeric(data$o_output_tree[[itree]]$lfMat))
      o_results_matrix[itree,7] = ifelse(data$o_output_tree[[itree]]$right == -1, -1, 
                                        as.numeric(data$o_output_tree[[itree]]$rtMat))
    }
    
    o_results_matrix[itree,8] = ifelse(length(as.numeric(data$o_output_tree[[itree]]$wilk_min))==0, 1, 
                                      as.numeric(data$o_output_tree[[itree]]$wilk_min))

    #: if is leaf, calculate y accordingly
    if (data$o_output_tree[[itree]]$left == -1 && data$o_output_tree[[itree]]$right == -1)
    {
      data$n_leafnodes_count <- data$n_leafnodes_count + 1
      n_y_size = nrow(data$o_output_tree[[itree]]$rowids_matrix)
      o_y_matrix = matrix(0, n_y_size, data$n_sample_y_cols)

      for (iy in 1:n_y_size)
      {
        n_row_id_y = data$o_output_tree[[itree]]$rowids_matrix[iy, 1]
        o_y_matrix[iy, ] = as.numeric(data$o_sample_data_y[n_row_id_y, ])
      }

      #: processing mapvalue
      if (data$n_mapvalue == "mean")
      {
        o_y_results_matrix[as.numeric(data$o_output_tree[[itree]]$id), ] = colMeans(o_y_matrix)
      }
      else if (data$n_mapvalue == "interval")
      {
        o_y_results_matrix[as.numeric(data$o_output_tree[[itree]]$id), 1:data$n_sample_y_cols] = apply(o_y_matrix, 2, min)
        o_y_results_matrix[as.numeric(data$o_output_tree[[itree]]$id), (data$n_sample_y_cols + 1):n_mapfile_cols] = apply(o_y_matrix, 2, max)
      }
    }
  }

  #: add matrix column names
  colnames(o_results_matrix) = c("NID", "xCol", "xVal", "lfNID", "rtNID", "lfMat", "rtMat", "wilk_min")

  if (data$n_mapvalue == "mean")
    colnames(o_y_results_matrix) = colnames(data$o_sample_data_y)
  else if (data$n_mapvalue == "interval")
    colnames(o_y_results_matrix) = c(colnames(data$o_sample_data_y), colnames(data$o_sample_data_y))

  #: store the tree and map files
  data$Tree <- o_results_matrix
  data$Map <- o_y_results_matrix
  data$totalNodes <- length(data$o_output_tree)
  data$leafNodes <- data$n_leafnodes_count
  data$cuttingActions <- data$n_cut_times
  data$mergingActions <- data$n_merge_times
  
  return(data)
}