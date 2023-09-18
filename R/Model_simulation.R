#################################################################
# Filename: 	Model_simulation.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
Model_simulation <- function(Testing_data, X, model)
{
  #:Number of predictants
  Num_Pre <- length(model[[1]]$weight)
  #:Name of predictants
  Nam_Pre <- model[[1]]$YName
  #:Number of rows for training
  Num_Tr <- nrow(model[[1]]$Training_sim)
  #:ensemble the prediction
  Training_sim <- lapply(model,function(x)x$Training_sim*x$weight)
  Training_sim_mat <- matrix(0,nrow = Num_Tr,ncol = Num_Pre)
  for(i in 1:Num_Pre)
  {
    Training_sim_mat[,i] <- apply(do.call(cbind,lapply(Training_sim,function(x)x[[i]])),1,sum)
  }
  colnames(Training_sim_mat) <- Nam_Pre
  Training_sim <- as.data.frame(Training_sim_mat)
  Validation_sim <- OOB_validation(model)
  Testing_sim <- SCE_Prediction(X_sample=Testing_data[,X],model)
  Testing_sim <- as.data.frame(Testing_sim)
  Output <- list(Training_sim,Validation_sim,Testing_sim)
  names(Output) <- c("Training","Validation","Testing")
  return(Output)
}