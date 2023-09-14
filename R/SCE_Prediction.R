#: predict one input
#: input: o_precitor = c(x1, x2, x3, ... )
#: output: o_predictant = c(y1, y2, y3, ... )
f_predict_one = function(o_precitor)
{
  o_predictant = c()

  #: initiate as c(0, 0, 0, ... )
  for (ip in 1:rSCA.env$n_y_cols_p)
  {
    o_predictant[ip] = 0
  }

  n_tree_index = 1
  while(n_tree_index <= rSCA.env$n_result_tree_rows_p)
  {
    if (rSCA.env$o_result_tree_p[n_tree_index, 4] == -1 && rSCA.env$o_result_tree_p[n_tree_index, 5] == -1) #make judgement based on 4 and 5 columns of tree file
    {
      #: if meets leaf, then get the mean y
      o_predictant = c(as.numeric(rSCA.env$o_mean_y_p[rSCA.env$o_result_tree_p[n_tree_index, 1], ])) #prediction based on map file
      O_Predictant_ID <- n_tree_index   #use to identify which node the predictant fall into
      break
    }

    #: if not, compare with the Xij
    #: if <= it, go to left node; otherwise go to right node

    #: get j (column index)
    n_j = rSCA.env$o_result_tree_p[n_tree_index, 2] #which column of x to use
    if (n_j <= 0)
    {
      #: if merged, just go to the left node
      n_tree_index = rSCA.env$o_result_tree_p[n_tree_index, 4]
      next
    }
    if (o_precitor[n_j] <= rSCA.env$o_result_tree_p[n_tree_index, 3]) #looking for which node the x belongs
    {
      n_tree_index = rSCA.env$o_result_tree_p[n_tree_index, 4]
    }
    else
    {
      n_tree_index = rSCA.env$o_result_tree_p[n_tree_index, 5]
    }
  }

  #: return
  return(o_predictant)
}

#: do precition to all input
f_predict = function()
{
  for (iPredictor in 1:rSCA.env$n_predictors_rows_p)
  {
    o_vector_predictor = c(as.numeric(rSCA.env$o_predictors_p[iPredictor, ]))
    o_vector_predictant = f_predict_one(o_vector_predictor)

    rSCA.env$o_predictants_p[iPredictor, ] <- o_vector_predictant
  }
}

# ---------------------------------------------------------------
# Main funtion
# ---------------------------------------------------------------
f_main_p = function()
{
  #cat("Processing...\t\t")
  f_predict()
  #cat("SUCCESS!\r\n")

  #: set column names + format output
  if (rSCA.env$n_model_type_p == "mean" || rSCA.env$n_model_type_p == "max" || rSCA.env$n_model_type_p == "min" || rSCA.env$n_model_type_p == "median")
  {
    colnames(rSCA.env$o_predictants_p) = colnames(rSCA.env$o_mean_y_p)
  }
  if (rSCA.env$n_model_type_p == "interval")
  {
    #: print col names
    s_colnames = colnames(rSCA.env$o_mean_y_p)
    n_midcol = rSCA.env$n_y_cols_p / 2
    for (icn in 1:n_midcol)
    {
      cat(s_colnames[icn], "\t\t", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
    cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)

    #: print result matrix
    for (ime in 1:nrow(rSCA.env$o_predictants_p))
    {
      o_vec_min = rSCA.env$o_predictants_p[ime, 1:n_midcol]
      o_vec_max = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
      for (iv in 1:n_midcol)
      {
        cat("[", o_vec_min[iv], ", ", o_vec_max[iv], "]\t", file = rSCA.env$filepath, sep = "", append = TRUE)
      }
      cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
  }
  if (rSCA.env$n_model_type_p == "radius")
  {
    #: print col names
    s_colnames = colnames(rSCA.env$o_mean_y_p)
    n_midcol = rSCA.env$n_y_cols_p / 2
    for (icn in 1:n_midcol)
    {
      cat(s_colnames[icn], "\t\t", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
    cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)

    #: print result matrix
    for (ime in 1:nrow(rSCA.env$o_predictants_p))
    {
      o_vec_mean = rSCA.env$o_predictants_p[ime, 1:n_midcol]
      o_vec_radius = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
      for (iv in 1:n_midcol)
      {
        cat("[", o_vec_mean[iv], " +/- ", o_vec_radius[iv], "]\t", file = rSCA.env$filepath, sep = "", append = TRUE)
      }
      cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
  }
  if (rSCA.env$n_model_type_p == "variation")
  {
    #: print col names
    s_colnames = colnames(rSCA.env$o_mean_y_p)
    n_midcol = rSCA.env$n_y_cols_p / 2
    for (icn in 1:n_midcol)
    {
      cat(s_colnames[icn], "\t\t", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
    cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)

    #: print result matrix
    for (ime in 1:nrow(rSCA.env$o_predictants_p))
    {
      o_vec_mean = rSCA.env$o_predictants_p[ime, 1:n_midcol]
      o_vec_sd = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
      for (iv in 1:n_midcol)
      {
        cat("[", o_vec_mean[iv], " +/- ", o_vec_sd[iv], "]\t", file = rSCA.env$filepath, sep = "", append = TRUE)
      }
      cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
  }
  if (rSCA.env$n_model_type_p == "random")
  {
    #: print col names
    s_colnames = colnames(rSCA.env$o_mean_y_p)
    n_midcol = rSCA.env$n_y_cols_p / 2
    for (icn in 1:n_midcol)
    {
      cat(s_colnames[icn], "\t\t", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
    cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)

    #: print result matrix
    for (ime in 1:nrow(rSCA.env$o_predictants_p))
    {
      o_vec_min = rSCA.env$o_predictants_p[ime, 1:n_midcol]
      o_vec_max = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
      for (iv in 1:n_midcol)
      {
        n_random_value = runif(1, o_vec_min[iv], o_vec_max[iv])
        cat(n_random_value, "\t", file = rSCA.env$filepath, sep = "", append = TRUE)
      }
      cat("\r\n", file = rSCA.env$filepath, sep = "", append = TRUE)
    }
  }
}

# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------

Inference <- function(x,Weak_L)
{
  rSCA.env$o_predictors_p = x
  rSCA.env$n_predictors_rows_p = nrow(x)

  #: read tree and map file
  rSCA.env$o_result_tree_p = Weak_L$Tree
  rSCA.env$n_result_tree_rows_p = nrow(rSCA.env$o_result_tree_p)
  rSCA.env$o_mean_y_p = Weak_L$Map
  rSCA.env$n_y_cols_p = ncol(rSCA.env$o_mean_y_p)
  rSCA.env$n_x_cols_p = ncol(rSCA.env$o_predictors_p)

  #: define the result file
  rSCA.env$filepath = paste("rsl_", Weak_L$ID, ".txt", sep = "")

  #: set model type
  rSCA.env$n_model_type_p = Weak_L$type

  #: start prediction
  rSCA.env$o_predictants_p <- matrix(0, rSCA.env$n_predictors_rows_p, rSCA.env$n_y_cols_p)

  f_main_p()
  o_result = rSCA.env$o_predictants_p
  return(o_result)
}
