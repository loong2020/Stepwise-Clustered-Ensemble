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

GOF <- function(obs,sim,Num_predictor,digits=4)
{
  sim[sim<=0] <- 0.1
  obs[obs<=0] <- 0.1
  GOF <- gof(obs=obs,sim=sim,digits=digits)
  log_nse <- NSE_equation(x=log(obs),y=log(sim))
  Adj_R2 <- Adj_rsq(obs,sim,Num_predictor)
  names(Adj_R2) <- "Adj_R2"
  names(log_nse) <- "Log.NSE"
  GOF <- rbind(GOF,Adj_R2,log_nse)
  colnames(GOF) <- "GOF"
  return(GOF)
}

Model_evaluation <- function(Predictant,obs_training,obs_testing,sim,Num_predictor,digits=3)
{
  Training_GOF <- GOF(obs=obs_training[,Predictant],sim=sim[["Training"]][,Predictant],Num_predictor=Num_predictor,digits=digits)
  Validation_GOF <- GOF(obs=obs_training[,Predictant],sim=sim[["Validation"]][,Predictant],Num_predictor=Num_predictor,digits=digits)
  Testing_GOF <- GOF(obs=obs_testing[,Predictant],sim=sim[["Testing"]][,Predictant],Num_predictor=Num_predictor,digits=digits)
  GOF_res <- data.frame(Training_GOF,Validation_GOF,Testing_GOF)
  colnames(GOF_res) <- c("Training","Validation","Testing")
  GOF_res <- format(GOF_res, scientific = FALSE)
  return(GOF_res)
}