# SCE: Stepwise Clustered Ensemble

## Overview

The SCE (Stepwise Clustered Ensemble) package provides implementation of Stepwise Clustered Ensemble (SCE) and Stepwise Cluster Analysis (SCA) methods for multivariate data analysis. These methods are particularly useful for handling complex, high-dimensional datasets and building robust predictive models.

## Installation

You can install the development version of SCE from GitHub:

```r
# install.packages("devtools")
devtools::install_github("loong2020/Stepwise-Clustered-Ensemble", 
                         build_vignettes = TRUE, 
                         force = TRUE)
```

## Main Functions

- `SCE()`: Main function for building a Stepwise Clustered Ensemble model
- `SCA()`: Stepwise Cluster Analysis (ensemble member of SCE)
- `Model_simulation()`: Perform SCE model preidction
- `SCA_tree_predict()`: Perform SCA model preidction
- `SCA_Model_evaluation()`: Evaluate model performance for SCA 
- `SCE_Model_evaluation()`: Evaluate model performance for SCA 
- `RFE_SCE()`: Recursive Feature Elimination for SCE
- `Wilks_importance()`: Calculate variable importance using Wilks' lambda

## Usage

For detailed usage examples and tutorials, please refer to the package vignettes. To view the vignettes:

```r
vignette("SCE-introduction", package = "SCE")
```

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