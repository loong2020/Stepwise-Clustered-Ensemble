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

- `SCE()`: Main function for Stepwise Clustered Ensemble analysis
- `SCA()`: Stepwise Cluster Analysis function
- `Model_evaluation()`: Evaluate model performance
- `Model_simulation()`: Perform model simulations
- `RFE_SCE()`: Recursive Feature Elimination for SCE
- `Training()`: Training function for the ensemble model
- `Prediction()`: Make predictions using the trained model
- `Wilks_importance()`: Calculate variable importance using Wilks' lambda

## Usage

For detailed usage examples and tutorials, please refer to the package vignettes.

## Documentation

Full documentation is available through the R help system:

```r
?SCE
?SCA
```

## License

This package is licensed under the GPL-3 License.

## Author

- Kailong Li (lkl98509509@gmail.com) 