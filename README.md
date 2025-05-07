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
```

### SCA (Single tree) Analysis
```r
# Load the example datasets
data(Streamflow_training_10var)
data(Streamflow_testing_10var)

# Define predictors and predictants
Predictors <- c("Prcp", "SRad", "Tmax", "Tmin", "VP", "smlt", "swvl1", "swvl2", "swvl3", "swvl4")
Predictants <- c("Flow")

# Perform SCA
set.seed(123)
model <- SCA(alpha = 0.05, 
            Training_data = Streamflow_training_10var, 
            X = Predictors, 
            Y = Predictants, 
            Nmin = 5, 
            resolution = 100)

# Calculate variable importance
importance <- SCA_importance(model)
print(importance)

# Make predictions
prediction <- SCA_tree_predict(Testing_data = Streamflow_testing_10var, model = model)

# Evaluate performance
performance <- SCA_Model_evaluation(Testing_data = Streamflow_testing_10var,
                                  Simulations = prediction,
                                  Predictant = Predictants)
print(performance)

Importance_ranking_sorted <- importance[order(-importance$Relative_Importance), ]
barplot(
  Importance_ranking_sorted$Relative_Importance,
  names.arg = Importance_ranking_sorted$Predictor,
  las = 2, # vertical labels
  col = "skyblue",
  main = "Variable Importance (SCE)",
  ylab = "Importance",
  xlab = "Predictor"
)
```

### SCE (Tree ensemble) Analysis
```r
# Build SCE model
set.seed(123)
Ensemble <- SCE(Training_data = Streamflow_training_10var,
               X = Predictors,
               Y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               Nmin = 5,
               Ntree = 40,
               alpha = 0.05,
               resolution = 100)

# Make predictions
Simulations <- Model_simulation(Testing_data = Streamflow_testing_10var, model = Ensemble)

# Evaluate model performance
Evaluation <- SCE_Model_evaluation(Testing_data = Streamflow_testing_10var,
                                 Training_data = Streamflow_training_10var,
                                 Simulations = Simulations,
                                 Predictant = Predictants,
                                 digits = 2)

# Calculate variable importance
importance <- Wilks_importance(Ensemble)
print(Evaluation)

Importance_ranking_sorted <- importance[order(-importance$Relative_Importance), ]
barplot(
  Importance_ranking_sorted$Relative_Importance,
  names.arg = Importance_ranking_sorted$Predictor,
  las = 2, # vertical labels
  col = "skyblue",
  main = "Variable Importance (SCE)",
  ylab = "Importance",
  xlab = "Predictor"
)
```

### Multiple Predictants Case
```r
# Define predictors and multiple predictants
# Load the example datasets
data(Air_quality_training)
data(Air_quality_testing)

Predictors <- c("SO2", "NO2", "CO", "O3", "TEMP", "PRES", "DEWP", "RAIN", "WSPM")
Predictants <- c("PM2.5", "PM10")

# Build and evaluate model
set.seed(123)
Ensemble <- SCE(Training_data = Air_quality_training,
               X = Predictors,
               Y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               Nmin = 5,
               Ntree = 40,
               alpha = 0.05,
               resolution = 100)

Simulations <- Model_simulation(Testing_data = Air_quality_testing, model = Ensemble)

Evaluation <- SCE_Model_evaluation(Testing_data = Air_quality_testing,
                                 Training_data = Air_quality_training,
                                 Simulations = Simulations,
                                 Predictant = Predictants)
print(Evaluation)

importance <- Wilks_importance(Ensemble)
print(Evaluation)

Importance_ranking_sorted <- importance[order(-importance$Relative_Importance), ]
barplot(
  Importance_ranking_sorted$Relative_Importance,
  names.arg = Importance_ranking_sorted$Predictor,
  las = 2, # vertical labels
  col = "skyblue",
  main = "Variable Importance (SCE)",
  ylab = "Importance",
  xlab = "Predictor"
)
```

### Recursive Feature Elimination
```r
# Load the example datasets
data(Streamflow_training_22var)
data(Streamflow_testing_22var)

# Define predictors and predictants
Predictors <- c(
  "Precipitation", "Radiation", "Tmax", "Tmin", "VP",
  "Precipitation_2Mon", "Radiation_2Mon", "Tmax_2Mon", "Tmin_2Mon", "VP_2Mon",
  "PNA", "Nino3.4", "IPO", "PDO",
  "PNA_lag1", "Nino3.4_lag1", "IPO_lag1", "PDO_lag1",
  "PNA_lag2", "Nino3.4_lag2", "IPO_lag2", "PDO_lag2"
)
Predictants <- c("Flow")

# Perform RFE
set.seed(123)
result <- RFE_SCE(
  Training_data = Streamflow_training_22var,
  Testing_data = Streamflow_testing_22var,
  Predictors = Predictors,
  Predictant = Predictants,
  Nmin = 5,
  Ntree = 48,
  alpha = 0.05,
  resolution = 1000,
  step = 3  # Number of predictors to remove at each iteration
)

# Plot Testing R² results
library(ggplot2)

# Extract Validation and Testing R² values
validation_r2 <- sapply(result[["performances"]], function(x) x["R2", "Validation"])
testing_r2 <- sapply(result[["performances"]], function(x) x["R2", "Testing"])
n_predictors <- result[["summary"]][["n_predictors"]]

# Create base R plot
plot(n_predictors, validation_r2, 
     type = "b",  # both points and lines
     col = "blue",
     pch = 16,    # filled circle point type
     xlim = rev(range(n_predictors)),  # reverse x-axis
     ylim = c(min(c(validation_r2, testing_r2)), max(c(validation_r2, testing_r2))),  # explicit y-axis limits
     xlab = "Number of Predictors",
     ylab = "R²",
     main = "Validation and Testing R² vs Number of Predictors")

# Add testing data
lines(n_predictors, testing_r2, type = "b", col = "red", pch = 16)

# Add legend
legend("bottomleft",
       legend = c("Validation", "Testing"),
       col = c("blue", "red"),
       pch = 16,
       lty = 1)
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
