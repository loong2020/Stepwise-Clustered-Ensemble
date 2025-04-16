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
- `Wilks_importance`: Calculate variable importance using Wilks' lambda

## Usage Examples

First, load the required data:

```r
# Load the example datasets
data(Training_input)
data(Testing_input)
```

### Basic SCA Analysis
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

# Make predictions
prediction <- SCA_tree_predict(Test_data = Testing_input,
                             X = Predictors,
                             model = model)

# Evaluate performance
performance <- SCA_Model_evaluation(Testing_data = Testing_input,
                                  Simulations = prediction,
                                  Predictant = Predictants,
                                  Num_predictor = length(Predictors),
                                  digits = 2)
print(performance)
```

### SCE Ensemble Analysis
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
                                 Num_predictor = length(Predictors),
                                 digits = 2)

# Calculate variable importance
Importance_ranking <- Wilks_importance(Ensemble)
print(Evaluation)
```

### Multiple Predictants
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
                                 Predictant = Predictants,
                                 Num_predictor = length(Predictors),
                                 digits = 2)
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
  alpha = 0.05,
  Nmin = 5,
  Ntree = 48,
  resolution = 100,
  metric = "nse",
  step = 1
)
```

## Documentation

Full documentation is available through the R help system:

```r
?SCE
?SCA
```

## License

This package is licensed under the GPL-3 License.

## Authors

- Kailong Li (lkl98509509@gmail.com) 
- Xiuquan Wang (xxwang@upei.ca)