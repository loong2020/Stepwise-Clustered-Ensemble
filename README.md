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
data(Streamflow_training_10var)
data(Streamflow_testing_10var)
```

### SCA (Single tree) Analysis
```r
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

Importance_ranking_sorted <- importance[order(-importance$Importance), ]
barplot(
  Importance_ranking_sorted$Importance,
  names.arg = Importance_ranking_sorted$col_index,
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
Simulations <- Model_simulation(Testing_data = Streamflow_testing_10var,
                              Training_data = Streamflow_training_10var,
                              model = Ensemble)

# Evaluate model performance
Evaluation <- SCE_Model_evaluation(Testing_data = Streamflow_testing_10var,
                                 Training_data = Streamflow_training_10var,
                                 Simulations = Simulations,
                                 Predictant = Predictants,
                                 digits = 2)

# Calculate variable importance
importance <- Wilks_importance(Ensemble)
print(Evaluation)

Importance_ranking_sorted <- importance[order(-importance$Importance), ]
barplot(
  Importance_ranking_sorted$Importance,
  names.arg = Importance_ranking_sorted$col_index,
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

Simulations <- Model_simulation(Testing_data = Air_quality_testing,
                              Training_data = Air_quality_training,
                              model = Ensemble)

Evaluation <- SCE_Model_evaluation(Testing_data = Air_quality_testing,
                                 Training_data = Air_quality_training,
                                 Simulations = Simulations,
                                 Predictant = Predictants)
print(Evaluation)

importance <- Wilks_importance(Ensemble)
print(Evaluation)

Importance_ranking_sorted <- importance[order(-importance$Importance), ]
barplot(
  Importance_ranking_sorted$Importance,
  names.arg = Importance_ranking_sorted$col_index,
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

# Create data frame for plotting
plot_data <- data.frame(
  n_predictors = rep(n_predictors, 2),
  r2 = c(validation_r2, testing_r2),
  type = rep(c("Validation", "Testing"), each = length(n_predictors))
)

# Create line plot with both Validation and Testing R²
ggplot(plot_data, aes(x = n_predictors, y = r2, color = type, group = type)) +
  geom_line() +
  geom_point() +
  scale_x_reverse(breaks = n_predictors) +  # Reverse x-axis and set breaks
  scale_color_manual(values = c("Validation" = "blue", "Testing" = "red")) +
  labs(
    x = "Number of Predictors",
    y = "R²",
    title = "Validation and Testing R² vs Number of Predictors",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
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
