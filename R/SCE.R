#################################################################
# Filename: 	SCE_Main.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
# History: 	2019/05/17		created by Kailong Li
#			      2019/09/18		added Wilks feature importance (WFI), by Kailong Li
#	     	    2023/03/10   	enabled weighted and non-weighted options for calculating WFI, by Kailong Li
##################################################################

# reference:
# Li, Kailong, Guohe Huang, and Brian Baetz
# Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling
# Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.
# https://doi.org/10.5194/hess-25-4947-2021

# ---------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------
SCE <- function(Training_data, X, Y, mfeature, Nmin, Ntree, alpha = 0.05) {
  # New environment
  rSCA.env <- new.env()

  # Define X and Y
  o_xdata <- data.frame(na.omit(Training_data[, X]))
  o_ydata <- data.frame(na.omit(Training_data[, Y]))

  colnames(o_xdata) <- X
  colnames(o_ydata) <- Y

  # Setup bootstrap
  Tree_name <- list()
  for(i in 1:Ntree) {
    Tree_index <- i
    Tree_name[[i]] <- paste("SCA_", Tree_index, sep = "_")
  }

  # Random feature
  sample_matrix <- matrix(rep(1:ncol(o_xdata), Ntree), Ntree, ncol(o_xdata), byrow = TRUE)
  Random_col_matrix <- matrix(0, nrow = Ntree, ncol = mfeature) # which columns will be used in the prediction
  Random_col <- list()
  set.seed(10)
  for (i in 1:Ntree) {
    Random_col_matrix[i,] <- sort(sample(sample_matrix[i,], mfeature))
    Random_col[[i]] <- Random_col_matrix[i,]
  }

  # Bootstrap
  temp <- o_xdata
  set.seed(10)  # Set a seed for reproducibility
  tree_list <- list()
  for (i in 1:Ntree) {
    # Create a bootstrap sample
    tree_list[[i]] <- sample(1:nrow(temp), replace = TRUE)
  }
  Bootst_rep <- mapply(function(x, y, z) list(Tree = x, mfeature = y, sample = z), x = Tree_name, y = Random_col, z = tree_list, SIMPLIFY = FALSE)

  # Parallel computing
  numcores <- parallel::detectCores()
  Clus <- parallel::makeCluster(numcores)

  # Source training and prediction scripts within the cluster
  parallel::clusterEvalQ(Clus, {
    source("Training.R", local = TRUE)
    source("Prediction.R", local = TRUE)
  })

  # Export variables to the workers
  parallel::clusterExport(Clus, c("o_xdata", "o_ydata", "Bootst_rep", "Nmin", "alpha", "rSCA.env"), envir = environment())

  # Evaluate variables on workers
  parallel::clusterEvalQ(Clus, {
    o_xdata
    o_ydata
    Bootst_rep
    Nmin
    alpha
    rSCA.env
  })

  # Parallel processing
  SCE_res <- parallel::parLapply(Clus, Bootst_rep, SCA, alpha = alpha, Nmin = Nmin)
  parallel::stopCluster(Clus)

  # Calculate the weights based on OOB_RSQ
  OOB_RSQ <- lapply(SCE_res, function(x) x$OOB_error)
  weight_OOB <- lapply(OOB_RSQ, function(x) log10(x / (1 - x)))
  weight_OOB <- data.frame(do.call(rbind, weight_OOB))

  # Re-scale the weight into [0,1]
  for(i in 1:ncol(weight_OOB)) {
    weight_OOB[, i] <- (weight_OOB[, i] - min(weight_OOB[, i])) / (max(weight_OOB[, i]) - min(weight_OOB[, i]))
    weight_OOB[, i] <- weight_OOB[, i] / sum(weight_OOB[, i])
  }

  weight_OOB <- apply(weight_OOB, 1, function(x) list(as.numeric(x)))
  weight_OOB <- lapply(weight_OOB, function(x) x[[1]])

  SCE_res <- mapply(function(x, y) c(x, list(weight = y)), x = SCE_res, y = weight_OOB, SIMPLIFY = FALSE)
  return(SCE_res)
}

SCE_Prediction <- function(X_sample,model)
{
  #:do the prediction
  Predicted <- lapply(model,function(m)SCA_Prediction(X_sample,m))
  model <- mapply(function(x,y) c(x,list(predicted=y)),x=model,y=Predicted,SIMPLIFY = FALSE)
  #:Number of predictants
  Num_Pre <- length(model[[1]]$weight)
  #:Name of predictants
  Nam_Pre <- model[[1]]$YName
  #:ensemble the prediction
  Weighted_prediction <- lapply(model,function(x)x$predicted*x$weight)
  Predictants <- matrix(0,nrow = nrow(X_sample),ncol = Num_Pre)
  for(i in 1:Num_Pre)
  {
    Predictants[,i] <- apply(do.call(cbind,lapply(Weighted_prediction,function(x)x[[i]])),1,sum)
  }
  colnames(Predictants) <- Nam_Pre
  return(Predictants)
}

OOB_validation <- function(model)
{
  #:Number of predictants
  Num_Pre <- length(model[[1]]$weight)
  #:Name of predictants
  Nam_Pre <- model[[1]]$YName
  #:OOB prediction
  OOB_predicted <- lapply(model,function(x) data.frame(ID=setdiff(1:length(x[["Sample"]]),unique(x[["Sample"]])),x[["OOB_sim"]],weight=matrix(x[["weight"]],nrow=1)))
  OOB_predicted <- do.call(rbind,OOB_predicted)
  OOB <- list()
  for (i in 1:Num_Pre) {
    DT_OOB <- OOB_predicted[,c(1,1+i,1+i+Num_Pre)]
    DT_OOB <- data.frame(DT_OOB)
    colnames(DT_OOB) <- c("ID", "Predicted", "weight")
    Predicted <- DT_OOB$Predicted
    weight <- DT_OOB$weight
    # Calculate weighted means for each ID manually
    ID_list <- unique(DT_OOB$ID)
    weighted_means <- sapply(ID_list, function(id) {
      subset <- DT_OOB[DT_OOB$ID == id, ]
      sum(subset$Predicted * subset$weight) / sum(subset$weight)
    })
    # Create a data frame for the results
    OOB[[i]] <- data.frame(ID = ID_list, predicted = weighted_means)
    OOB[[i]] <- OOB[[i]][order(OOB[[i]]$ID), ]
  }
  # Combine the results into a single data frame
  OOB <- data.frame(do.call(cbind,lapply(OOB,function(x)x[,-1])))
  colnames(OOB) <- Nam_Pre
  return(OOB)
}

SCA <- function(alpha,Input,Nmin)
{
  #: retrieve the out-of-bag data
  OOB <- setdiff(1:nrow(o_xdata),unique(Input[[3]]))
  OOB_X <- o_xdata[OOB,Input[[2]]]
  OOB_y <- data.frame(o_ydata[OOB,])
  colnames(OOB_y) <- colnames(o_ydata)

  #: Training set
  Trainging_X <- o_xdata[,Input[[2]]]

  #: generate bootstrap replicate
  o_xdata <- data.frame(o_xdata[Input[[3]],Input[[2]]])
  o_ydata <- data.frame(o_ydata[Input[[3]],])
  colnames(o_ydata) <- colnames(OOB_y)

  #: assign name
  if(is.null(colnames(o_ydata))==TRUE)
  {
    o_ydata <- data.frame(response=o_ydata)
  }else{
    o_ydata <- data.frame(o_ydata)
  }

  #: add randomize noise
  ex.runif <- function(n, exclmin, exclmax, min, max) # increase sample variance
  {
    # ex.runif() excludes the specific range 'exclmin' to 'exclmax'
    temp <- runif(n, min = min, max = max)
    for (i in 1:length(temp))
    {
      while (((exclmin <= temp[i]) && (temp[i] <= exclmax))==TRUE)
      {
        temp[i] <- runif(1, min = min, max = max)
      }
    }
    return(temp)
  }

  #:load required functions
  rsq <- function(x, y) summary(lm(y~x))$r.squared #R squared function
  Random_number <- rep(ex.runif(nrow(o_xdata), -0.01, 0.01, -0.05, 0.05), ncol(o_xdata)) #random number generator
  Random_number_Matrix <- matrix(Random_number,nrow=nrow(o_xdata),ncol=ncol(o_xdata))
  x_mean <- as.numeric(apply(o_xdata,2,mean))
  y_mean <- as.numeric(apply(o_ydata,2,mean))
  x_random <- t(apply(Random_number_Matrix,1,function(x) x*x_mean))
  if(ncol(o_ydata)==1)
  {
    y_random <- data.frame(dat=Random_number_Matrix[,1]*y_mean)
  } else {
    y_random <- t(data.frame(dat=apply(data.frame(dat=Random_number_Matrix[,1:ncol(o_ydata)]),1,function(x) x*y_mean))) # build random data for each column
  }
  o_xdata <- o_xdata+x_random
  o_ydata <- o_ydata+y_random
  o_xdata[is.na(o_xdata)] <- 0 #if sample is too small, bootstrap will not cover all values, change NA to zero
  o_ydata[is.na(o_ydata)] <- 0

  #: basic info
  n_xdatarow = nrow(o_xdata)
  n_ydatarow = nrow(o_ydata)

  rSCA.env$n_alpha = alpha
  rSCA.env$n_mapvalue = "mean"

  #: model list
  o_model = list()

  #: start training
  rSCA.env$o_sample_data_x = o_xdata
  rSCA.env$n_sample_size = nrow(rSCA.env$o_sample_data_x)
  rSCA.env$n_sample_x_cols = ncol(rSCA.env$o_sample_data_x)
  rSCA.env$o_sample_data_y = o_ydata
  rSCA.env$n_sample_y_cols = ncol(rSCA.env$o_sample_data_y)
  rSCA.env$s_tree_file = paste("tree_",Input[[1]],".txt",sep = "")
  rSCA.env$s_map_file = paste("map_",Input[[1]],".txt",sep = "")

  #: start SCA modeling
  do_cluster(Nmin)
  o_model$XName = colnames(o_xdata)
  o_model$YName = colnames(o_ydata)
  o_model$Tree = rSCA.env$Tree
  o_model$Map = rSCA.env$Map
  o_model$type = rSCA.env$n_mapvalue
  o_model$totalNodes = length(rSCA.env$o_output_tree)
  o_model$leafNodes = rSCA.env$n_leafnodes_count
  o_model$cuttingActions = rSCA.env$n_cut_times
  o_model$mergingActions = rSCA.env$n_merge_times
  o_model$Feature = Input[[2]]
  o_model$Sample = Input[[3]]
  #: return the model
  OOB_sim <- Inference(x=OOB_X,Weak_L=o_model)
  OOB_sim <- data.frame(OOB_sim)
  RSQ_OOB <- do.call(c,mapply(rsq,x=OOB_sim,y=data.frame(OOB_y),SIMPLIFY = FALSE))
  o_model$OOB_error = as.numeric(RSQ_OOB)
  o_model$OOB_sim = OOB_sim

  #: start inference
  Training_sim <- Inference(x=Trainging_X,Weak_L=o_model)
  Training_sim <- data.frame(do.call(cbind,data.frame(Training_sim)))
  colnames(Training_sim) <- o_model$YName
  o_model$Training_sim = Training_sim
  return(o_model)
}

SCA_Prediction <- function(X_sample,model)
{
  #:assign an environment to the global environment
  rSCA.env = new.env()
  assign("rSCA.env",rSCA.env,envir=globalenv())
  #: Prediction set
  Prediction_X <- X_sample[,model$Feature]
  Predicted <- Inference(x=Prediction_X,Weak_L=model)
  Predicted <- data.frame(Predicted)
  colnames(Predicted) <- model$YName
  #: remove the environment from the global environment
  rm(rSCA.env,envir=globalenv())
  return(Predicted)
}
