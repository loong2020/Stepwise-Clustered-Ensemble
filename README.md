# Stepwise-Clustered-Ensemble

SCE (Stepwise Clustered Ensemble) is a versatile tool for statistical multivariate prediction and inference. This R package facilitates the generation of an ensemble of tree-structured models that analyze the interrelationships between independent (input) variables and the variance in dependent (output) variables. Each ensemble member within SCE is  designed to evaluate these relationships. When the number of ensemble members increases, the variation of the unstable predictions is reduced, thereby increasing the predictive accuracy and inference robustness.

One of the standout features of SCE is its ability to simultaneously handle multiple inputs to predict multiple outputs, offering a comprehensive analysis that goes a step further than conventional tree-structured tools like random forest. Another key feature of SCE is that it offers a sophisticated ranking of input variables based on a statistical F-test, serving as a reliable method for determining the dominant controlling effects within a study system. This approach stands as a beacon of reliability as it doesn't hinge on evaluation metrics for ranking, significantly reducing the potential bias often introduced by these metrics. 

Author:
Kailong li  <lkl98509509@gmail.com>

references:
Li, Kailong, Guohe Huang, and Brian Baetz. Development of a Wilks feature importance method with improved variable rankings for supporting hydrological inference and modelling. Hydrology and Earth System Sciences 25.9 (2021): 4947-4966.

How to use:
## Install SCE package 
library(devtools)

devtools::install_github("loong2020/Stepwise-Clustered-Ensemble")

## Load SCE package and supporting packages 
library(SCE)

#: parallel computing package

library(parallel)

#: goodness-of-fit package

library(hydroGOF) 

## Training data file
data("Training_input")
## Testing data file
data("Testing_input")

## Define independent (x) and dependent (y) variables
Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2")

Predictants <- c("swvl3","swvl4")

## Build the SCE model
Model <- SCE(Training_data=Training_input, X=Predictors, Y=Predictants, mfeature=round(0.5*length(Predictors)),Nmin=5,Ntree=48,alpha=0.05)

## Make predictions
Results <- Model_simulation(Testing_data=Testing_input, X=Predictors, model=Model)

## Evaluate the model over a specified predictant
Evaluation_swvl3 <- Model_evaluation(Predictant="swvl3",obs_training=Training_input,obs_testing=Testing_input,sim=Results,Num_predictor=length(Predictors),digits=2)

Evaluation_swvl4 <- Model_evaluation(Predictant="swvl4",obs_training=Training_input,obs_testing=Testing_input,sim=Results,Num_predictor=length(Predictors),digits=2)

## Rank the importance of predictors
Importance_score <- Wilks_importance(Model,OOB_weight=FALSE)
