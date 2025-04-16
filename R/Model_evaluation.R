#################################################################
# Filename: 	Model_evaluation.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
#: load the function
rsq <- function(x, y) summary(lm(y~x))$r.squared #R squared function
Adj_rsq <- function(x,y,z) 1-(1-rsq(x,y))*(length(y)-1)/(length(y)-z-1) #adjusted R squared function
NSE_equation <- function(x, y) {1-(sum((x - y)^2)/ sum((x - mean(x))^2))} #obs = x, sim = y
KGE_equation <- function(x,y) {
  # x = observed, y = simulated
  r <- cor(x, y)  # correlation coefficient
  alpha <- sd(y)/sd(x)  # ratio of standard deviations
  beta <- mean(y)/mean(x)  # ratio of means
  KGE <- 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2)
  return(KGE)
} #obs = x, sim = y

GOF <- function(obs,sim,Num_predictor,digits=4)
{
  sim[sim<=0] <- 0.0001
  obs[obs<=0] <- 0.0001
  mae <- mean(abs(obs-sim))
  rmse <- sqrt(mean((obs-sim)^2))
  nse <- NSE_equation(x=obs,y=sim)
  log_nse <- NSE_equation(x=log(obs),y=log(sim))
  Adj_R2 <- Adj_rsq(obs,sim,Num_predictor)
  kge <- KGE_equation(x=obs,y=sim)
  names(Adj_R2) <- "Adj_R2"
  names(log_nse) <- "Log.NSE"
  GOF <- rbind(mae,rmse,nse,log_nse,Adj_R2,kge)
  GOF <- format(GOF, scientific = FALSE, digits = digits)
  GOF <- as.matrix(GOF)
  colnames(GOF) <- "GOF"
  return(GOF)
}

SCA_Model_evaluation <- function(Testing_data, Simulations, Predictant, Num_predictor, digits=3)
{
  # Input validation
  if (!is.data.frame(Testing_data)) {
    stop("Testing_data must be a data frame")
  }
  
  if (!is.list(Simulations) || !"Testing_sim" %in% names(Simulations)) {
    stop("Simulations must be a list with 'Testing_sim' component")
  }
  
  if (!is.character(Predictant) || length(Predictant) == 0) {
    stop("Predictant must be a non-empty character vector")
  }
  
  if (!is.numeric(Num_predictor) || Num_predictor < 0) {
    stop("Num_predictor must be a positive number")
  }
  
  # Check if predictants exist in data
  missing_predictants <- Predictant[!Predictant %in% colnames(Testing_data)]
  if (length(missing_predictants) > 0) {
    stop(sprintf("The following predictants are not found in Testing_data: %s", 
                paste(missing_predictants, collapse = ", ")))
  }
  
  # Check if predictants exist in simulations
  if (!all(Predictant %in% colnames(Simulations$Testing_sim))) {
    missing_predictants <- Predictant[!Predictant %in% colnames(Simulations$Testing_sim)]
    stop(sprintf("The following predictants are not found in Testing_sim: %s", 
                paste(missing_predictants, collapse = ", ")))
  }
  
  # Check row count match
  if (nrow(Testing_data) != nrow(Simulations$Testing_sim)) {
    stop("The number of rows in Testing_data must be equal to the number of rows in Testing_sim")
  }
  
  # For SCA, we only evaluate testing performance
  all_results <- list()
  for(pred in Predictant) {
    Testing_GOF <- GOF(obs=Testing_data[,pred], sim=Simulations$Testing_sim[,pred], Num_predictor=Num_predictor, digits=digits)
    GOF_res <- data.frame(Testing=Testing_GOF)
    colnames(GOF_res) <- "Testing"
    all_results[[pred]] <- GOF_res
  }
  
  if(length(Predictant) == 1) {
    return(all_results[[1]])
  }
  return(all_results)
}

SCE_Model_evaluation <- function(Testing_data, Training_data, Simulations, Predictant, Num_predictor, digits=3)
{
  # Input validation
  if (!is.data.frame(Testing_data) || !is.data.frame(Training_data)) {
    stop("Testing_data and Training_data must be data frames")
  }
  
  if (!is.list(Simulations) || !all(c("Training", "Validation", "Testing") %in% names(Simulations))) {
    stop("Simulations must be a list with 'Training', 'Validation', and 'Testing' components")
  }
  
  if (!is.character(Predictant) || length(Predictant) == 0) {
    stop("Predictant must be a non-empty character vector")
  }
  
  if (!is.numeric(Num_predictor) || Num_predictor < 0) {
    stop("Num_predictor must be a positive number")
  }
  
  # Check if predictants exist in data
  missing_predictants <- Predictant[!Predictant %in% colnames(Testing_data)]
  if (length(missing_predictants) > 0) {
    stop(sprintf("The following predictants are not found in the data: %s", 
                paste(missing_predictants, collapse = ", ")))
  }
  
  # Check if predictants exist in simulations
  for (pred in Predictant) {
    for (set in c("Training", "Validation", "Testing")) {
      if (!pred %in% colnames(Simulations[[set]])) {
        stop(sprintf("Predictant '%s' not found in %s simulations", pred, set))
      }
    }
  }

  # Total training data rows much be equal to the number of rows in the training simulations
  if (nrow(Training_data) != nrow(Simulations[["Training"]])) {
    stop("The number of rows in Training_data must be equal to the number of rows in the Training simulations")
  }

  # total training data rows must be equal to the number of rows in the validation simulations
  if (nrow(Training_data) != nrow(Simulations[["Validation"]])) {
    stop("The number of rows in Training_data must be equal to the number of rows in the Validation simulations")
  }

  # total testing data rows must be equal to the number of rows in the testing simulations
  if (nrow(Testing_data) != nrow(Simulations[["Testing"]])) {
    stop("The number of rows in Testing_data must be equal to the number of rows in the Testing simulations")
  }
  
  # Initialize a list to store results for each predictant
  all_results <- list()
  
  # Loop through each predictant
  for(pred in Predictant) {
    # Calculate GOF for each predictant
    Training_GOF <- GOF(obs=Training_data[,pred],sim=Simulations[["Training"]][,pred],Num_predictor=Num_predictor,digits=digits)
    Validation_GOF <- GOF(obs=Training_data[,pred],sim=Simulations[["Validation"]][,pred],Num_predictor=Num_predictor,digits=digits)
    Testing_GOF <- GOF(obs=Testing_data[,pred],sim=Simulations[["Testing"]][,pred],Num_predictor=Num_predictor,digits=digits)
    
    # Combine results for this predictant
    GOF_res <- data.frame(Training_GOF,Validation_GOF,Testing_GOF)
    colnames(GOF_res) <- c("Training","Validation","Testing")
    
    # Add to results list with predictant name
    all_results[[pred]] <- GOF_res
  }
  
  # If there's only one predictant, return the single result
  if(length(Predictant) == 1) {
    return(all_results[[1]])
  }
  
  # For multiple predictants, return the list of results
  return(all_results)
}