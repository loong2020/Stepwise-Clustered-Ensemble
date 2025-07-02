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
- `SCE()`: Build a Stepwise Clustered Ensemble model
- `SCA()`: Build a Stepwise Cluster Analysis model (single tree)

### Prediction and Evaluation
- `Model_simulation()`: Perform SCE model prediction
- `SCA_tree_predict()`: Perform SCA model prediction
- `SCE_Model_evaluation()`: Evaluate SCE model performance
- `SCA_Model_evaluation()`: Evaluate SCA model performance

### Feature Selection and Importance
- `RFE_SCE()`: Recursive Feature Elimination for SCE
- `Wilks_importance()`: Calculate variable importance for SCE using Wilks' lambda
- `SCA_importance()`: Calculate variable importance for a single SCA tree

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
- `evaluate()`: Evaluate model performance (testing only)

### Quick Start with S3 Methods
```r
# Build models
sce_model <- SCE(Training_data = data, X = predictors, Y = predictants, ...)
sca_model <- SCA(Training_data = data, X = predictors, Y = predictants, ...)

# Use S3 methods
print(sce_model)           # Display model info
summary(sce_model)         # Detailed summary
predictions <- predict(sce_model, newdata)  # Make predictions
imp_ranking <- importance(sce_model)  # Calculate variable importance
evaluation <- evaluate(sce_model, Testing_data, Training_data, Predictant)  # Evaluate model

# Check available methods
methods(class = "SCE")
methods(class = "SCA")
```

## Available Datasets

The package includes several datasets for demonstration and testing:

### Streamflow Datasets
- **Basic datasets (10 variables)**: `Streamflow_training_10var`, `Streamflow_testing_10var`
  - Contains hydrological and meteorological variables
  - Suitable for introductory examples and basic modeling
- **Extended datasets (22 variables)**: `Streamflow_training_22var`, `Streamflow_testing_22var`
  - Includes climate indices (IPO, Nino3.4, PDO, PNA) with lagged versions
  - Suitable for advanced modeling and research applications

### Air Quality Datasets
- `Air_quality_training`, `Air_quality_testing`
  - Contains air quality monitoring data
  - Useful for environmental modeling examples

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

# Use S3 methods
print(model)
summary(model)

# Calculate variable importance
Imp_ranking <- importance(model)
print(Imp_ranking)

# Make predictions
prediction <- predict(model, Streamflow_testing_10var)

# Evaluate performance
performance <- evaluate(
  object = model,
  Testing_data = Streamflow_testing_10var,
  Predictant = Predictants
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
Ensemble <- SCE(Training_data = Streamflow_training_10var,
               X = Predictors,
               Y = Predictants,
               mfeature = round(0.5 * length(Predictors)),
               Nmin = 5,
               Ntree = 40,
               alpha = 0.05,
               resolution = 100)

# Use S3 methods
print(Ensemble)
summary(Ensemble)

# Make predictions
predictions <- predict(Ensemble, Streamflow_testing_10var)
cat("Prediction components:", names(predictions), "\n")
cat("Testing predictions dimensions:", dim(predictions$Testing), "\n")

# Calculate variable importance
Imp_ranking <- importance(Ensemble)

# Evaluate model performance
evaluation <- evaluate(
  object = Ensemble,
  Testing_data = Streamflow_testing_10var,
  Training_data = Streamflow_training_10var,
  Predictant = Predictants,
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

# Use S3 methods
print(Ensemble)
summary(Ensemble)

# Make predictions
predictions <- predict(Ensemble, Air_quality_testing)

# Calculate variable importance
Imp_ranking <- importance(Ensemble)

# Evaluate model performance
evaluation <- evaluate(
  object = Ensemble,
  Testing_data = Air_quality_testing,
  Training_data = Air_quality_training,
  Predictant = Predictants
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

# Convert to numeric (remove any formatting/spaces)
validation_r2 <- as.numeric(validation_r2)
testing_r2 <- as.numeric(testing_r2)

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

# S3 methods
?predict.SCE
?predict.SCA
?importance.SCE
?importance.SCA
?evaluate.SCE
?evaluate.SCA
?print.SCE
?print.SCA
?summary.SCE
?summary.SCA

# Traditional functions (for advanced users)
?Model_simulation
?SCA_tree_predict
?SCA_Model_evaluation
?SCE_Model_evaluation
?RFE_SCE
?Wilks_importance
?SCA_importance
```

## License

This package is licensed under the GPL-3 License.

## Authors

- Kailong Li (lkl98509509@gmail.com) 
- Xiuquan Wang (xxwang@upei.ca)