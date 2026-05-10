# SCE: Stepwise Clustered Ensemble

## Overview

The SCE (Stepwise Clustered Ensemble) package provides implementation of Stepwise Clustered Ensemble (SCE) and Stepwise Cluster Analysis (SCA) methods for multivariate data analysis. These methods are particularly useful for handling complex, high-dimensional datasets and building robust predictive models.

The package supports proper S3 object-oriented programming, providing dedicated output classes with associated methods for `print`, `summary`, `predict`, `importance`, and `evaluate`.

## Installation

Install SCE from CRAN:

```r
install.packages("SCE")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("loong2020/Stepwise-Clustered-Ensemble")
```

## Core Functions

### Main Modeling Functions
- `sce()`: Build a Stepwise Clustered Ensemble model
- `sca()`: Build a Stepwise Cluster Analysis model (single tree)

### Prediction and Evaluation
- `model_simulation()`: Perform SCE model prediction (also invoked via `predict()` on `sce` objects)
- `sca_tree_predict()`: Perform SCA model prediction (also invoked via `predict()` on `sca` objects)
- `sce_model_evaluation()`: Evaluate SCE model performance
- `sca_model_evaluation()`: Evaluate SCA model performance

### Feature Selection and Importance
- `rfe_sce()`: Recursive Feature Elimination for SCE
- `wilks_importance()`: Calculate variable importance for SCE using Wilks' lambda
- `sca_importance()`: Calculate variable importance for a single SCA tree

## S3 Classes and Methods

The package provides S3 classes for both SCE and SCA models with convenient methods:

### SCE Class Methods
- `print()`: Display model information and performance metrics
- `summary()`: Detailed model summary with statistics
- `predict()`: Make predictions on new data (returns Training, Validation, and Testing predictions)
- `importance()`: Calculate variable importance using Wilks' lambda
- `evaluate()`: Evaluate model performance (training, validation, and testing)

### SCA Class Methods
- `print()`: Display tree structure and variable information
- `summary()`: Detailed tree summary with statistics
- `predict()`: Make predictions on new data
- `importance()`: Calculate variable importance
- `evaluate()`: Evaluate model performance (training and testing)

### Quick Start with S3 Methods
```r
# Build models
sce_model <- sce(training_data = data, x = predictors, y = predictants, ...)
sca_model <- sca(training_data = data, x = predictors, y = predictants, ...)

# Use S3 methods
print(sce_model)           # Display model info
summary(sce_model)         # Detailed summary
predictions <- predict(sce_model, newdata)  # Make predictions
imp_ranking <- importance(sce_model)  # Calculate variable importance
evaluation <- evaluate(sce_model, testing_data, training_data)  # Evaluate model

# Check available methods
methods(class = "sce")
methods(class = "sca")
```

## Available Datasets

The package includes several datasets for demonstration and testing:

### Streamflow Datasets
- **Basic (10 variables)**: `streamflow_training_10var`, `streamflow_testing_10var` — hydrological and meteorological predictors.
- **Extended (22 variables)**: `streamflow_training_22var`, `streamflow_testing_22var` — the basic predictors plus lagged climate indices (IPO, Nino3.4, PDO, PNA).

### Air Quality Datasets
- `air_quality_training`, `air_quality_testing` — eleven variables from air-quality monitoring stations.

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
data(streamflow_training_10var)
data(streamflow_testing_10var)

# Define predictors and predictants
Predictors <- c("Prcp", "SRad", "Tmax", "Tmin", "VP", "smlt", "swvl1", "swvl2", "swvl3", "swvl4")
Predictants <- c("Flow")

# Perform SCA
set.seed(123)
model <- sca(alpha = 0.05, 
            training_data = streamflow_training_10var, 
            x = Predictors, 
            y = Predictants, 
            nmin = 5, 
            resolution = 100)

# Use S3 methods
print(model)
summary(model)

# Calculate variable importance
Imp_ranking <- importance(model, digits = 2)
print(Imp_ranking)

# Make predictions
prediction <- predict(model, streamflow_testing_10var)

# Evaluate performance
performance <- evaluate(
  object = model,
  testing_data = streamflow_testing_10var,
  training_data = streamflow_training_10var
)
print(performance)

Importance_ranking_sorted <- Imp_ranking[order(-Imp_ranking$Relative_Importance), ]
barplot(
  Importance_ranking_sorted$Relative_Importance,
  names.arg = Importance_ranking_sorted$Predictor,
  las = 2, # vertical labels
  col = "skyblue",
  main = "Variable Importance (SCA)",
  ylab = "Importance",
  xlab = "Predictor"
)
```

### SCE (Tree ensemble) Analysis
```r
# Build SCE model
set.seed(123)
Ensemble <- sce(training_data = streamflow_training_10var,
               x = Predictors,
               y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               nmin = 5,
               ntree = 40,
               alpha = 0.05,
               resolution = 100)

# Use S3 methods
print(Ensemble)
summary(Ensemble)

# Make predictions
predictions <- predict(Ensemble, streamflow_testing_10var)
cat("Prediction components:", names(predictions), "\n")
cat("Testing predictions dimensions:", dim(predictions$Testing), "\n")

# Calculate variable importance
Imp_ranking <- importance(Ensemble, digits = 2)

# Evaluate model performance
evaluation <- evaluate(
  object = Ensemble,
  testing_data = streamflow_testing_10var,
  training_data = streamflow_training_10var,
  digits = 3
)
print(evaluation)

Importance_ranking_sorted <- Imp_ranking[order(-Imp_ranking$Relative_Importance), ]
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
data(air_quality_training)
data(air_quality_testing)

Predictors <- c("SO2", "NO2", "CO", "O3", "TEMP", "PRES", "DEWP", "RAIN", "WSPM")
Predictants <- c("PM2.5", "PM10")

# Build and evaluate model
set.seed(123)
Ensemble <- sce(training_data = air_quality_training,
               x = Predictors,
               y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               nmin = 5,
               ntree = 40,
               alpha = 0.05,
               resolution = 100)

# Use S3 methods
print(Ensemble)
summary(Ensemble)

# Make predictions
predictions <- predict(Ensemble, air_quality_testing)

# Calculate variable importance
Imp_ranking <- importance(Ensemble, digits = 2)

# Evaluate model performance
evaluation <- evaluate(
  object = Ensemble,
  testing_data = air_quality_testing,
  training_data = air_quality_training
)
print(evaluation)

Importance_ranking_sorted <- Imp_ranking[order(-Imp_ranking$Relative_Importance), ]
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
data(streamflow_training_22var)
data(streamflow_testing_22var)

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
set.seed(1)
result <- rfe_sce(
  training_data = streamflow_training_22var,
  testing_data = streamflow_testing_22var,
  predictors = Predictors,
  predictant = Predictants,
  nmin = 5,
  ntree = 48,
  alpha = 0.05,
  resolution = 1000,
  step = 3  # Number of predictors to remove at each iteration
)

# Plot RFE results
plot_rfe(result)
```

## Documentation

Full documentation is available through the R help system:

```r
# Core functions
?sce
?sca
?rfe_sce
?plot_rfe

# S3 methods (see also ?predict, ?importance, ?evaluate, ?print, ?summary)
?predict.sce
?predict.sca
?importance.sce
?importance.sca
?evaluate.sce
?evaluate.sca
?print.sce
?print.sca
?summary.sce
?summary.sca
```

## License

This package is licensed under the GPL-3 License.

## Authors

- Kailong Li (lkl98509509@gmail.com) 
- Xiuquan Wang (xxwang@upei.ca)
