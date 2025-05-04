# SCE: Stepwise Clustered Ensemble

## Overview

The SCE (Stepwise Clustered Ensemble) package provides implementation of Stepwise Clustered Ensemble (SCE) and Stepwise Cluster Analysis (SCA) methods for multivariate data analysis. These methods are particularly useful for handling complex, high-dimensional datasets and building robust predictive models.

## Installation

You can install the development version of SCE from GitHub:

```r
# install.packages("devtools")
devtools::install_github("loong2020/Stepwise-Clustered-Ensemble")
```

## Main Functions

- `SCE`: Main function for building a Stepwise Clustered Ensemble model
- `SCA`: Stepwise Cluster Analysis (ensemble member of SCE)
- `Model_simulation`: Perform SCE model prediction
- `SCA_tree_predict`: Perform SCA model prediction
- `SCA_Model_evaluation`: Evaluate model performance for SCA 
- `SCE_Model_evaluation`: Evaluate model performance for SCE
- `RFE_SCE`: Recursive Feature Elimination for SCE
- `Wilks_importance`: Calculate variable importance for SCE using Wilks' lambda
- `SCA_importance`: Calculate variable importance for a single SCA tree

## Usage Examples

First, load the required packages and data:

```r
# Load required packages
library(SCE)
library(parallel)

# Load the example datasets
data(Training_input)
data(Testing_input)
```

### SCA (Single tree) Analysis
```r
# Define predictors and predictants
Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3","swvl4")
Predictants <- c("Flow")

# Perform SCA
set.seed(123)
model <- SCA(alpha = 0.05, 
            Training_data = Training_input, 
            X = Predictors, 
            Y = Predictants, 
            Nmin = 5, 
            resolution = 100)

# Calculate variable importance
importance <- SCA_importance(model)
print(importance)

# Make predictions
prediction <- SCA_tree_predict(Testing_data = Testing_input, model = model)

# Evaluate performance
performance <- SCA_Model_evaluation(Testing_data = Testing_input,
                                  Simulations = prediction,
                                  Predictant = Predictants)
print(performance)
```

### SCE (Tree ensemble) Analysis
```r
# Build SCE model
set.seed(123)
Ensemble <- SCE(Training_data = Training_input,
               X = Predictors,
               Y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               Nmin = 5,
               Ntree = 40,
               alpha = 0.05,
               resolution = 100)

# Make predictions
Simulations <- Model_simulation(Testing_data = Testing_input,
                              Training_data = Training_input,
                              model = Ensemble)

# Evaluate model performance
Evaluation <- SCE_Model_evaluation(Testing_data = Testing_input,
                                 Training_data = Training_input,
                                 Simulations = Simulations,
                                 Predictant = Predictants,
                                 digits = 2)

# Calculate variable importance
Importance_ranking <- Wilks_importance(Ensemble)
print(Evaluation)
```

### Multiple Predictants Case
```r
# Define predictors and multiple predictants
Predictors <- c("Prcp","SRad","Tmax","Tmin","VP","smlt","swvl1","swvl2","swvl3")
Predictants <- c("Flow","swvl4")

# Build and evaluate model
set.seed(123)
Ensemble <- SCE(Training_data = Training_input,
               X = Predictors,
               Y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               Nmin = 5,
               Ntree = 40,
               alpha = 0.05,
               resolution = 100)

Simulations <- Model_simulation(Testing_data = Testing_input,
                              Training_data = Training_input,
                              model = Ensemble)

Evaluation <- SCE_Model_evaluation(Testing_data = Testing_input,
                                 Training_data = Training_input,
                                 Simulations = Simulations,
                                 Predictant = Predictants)
print(Evaluation)
```

### Recursive Feature Elimination
```r
# Perform RFE
result <- RFE_SCE(
  Training_data = Training_input,
  Testing_data = Testing_input,
  Predictors = Predictors,
  Predictant = Predictants,
  Nmin = 5,
  Ntree = 48,
  mfeature = round(0.5 * length(Predictors)),
  alpha = 0.05,
  resolution = 1000,
  step = 1  # Number of predictors to remove at each iteration
)
```

## Documentation

Full documentation is available through the R help system:

```r
# Core functions
?SCE
?SCA
?Model_simulation
?SCA_tree_predict

# Evaluation functions
?SCA_Model_evaluation
?SCE_Model_evaluation

# Feature selection and importance
?RFE_SCE
?Wilks_importance
?SCA_importance
```

## License

This package is licensed under the GPL-3 License.

## Authors

- Kailong Li (lkl98509509@gmail.com) 
- Xiuquan Wang (xxwang@upei.ca)