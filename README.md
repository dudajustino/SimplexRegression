
# SimplexRegression: A Package for Simplex Regression with parametric and fixed mean link function

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R package and associated data for the scientific article:

“Simplex regression with a flexible logit link: inference and application to cross-country impunity data” by Justino, M. E. C. and Cribari-Neto, F.

## 📑 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [✨ Key Features](#-key-features)
- [📂 Repository Structure](repository-structure)
- [🛠️ Installation](#-installation)
- [Quick Start](#-quick-start)
- [Usage Examples](#-usage-examples)
  - [Basic Model](#basic-model)
  - [Parametric Mean Link Functions](#parametric-mean-link-functions)
  - [Model Diagnostics](#model-diagnostics)
  - [Model Selection](#model-selection)
- [Vignettes](#-vignettes)
- [Real Data Application](#-real-data-application)
- [Functions Reference](#-functions-reference)
- [🤝 Contributing](#-contributing)
- [References](#-references)
- [📄 License](#-license)
- [Citation](#-citation)
- [📬 Contact](#-contact)

---

## Installation

You can install the development version of SimplexRegression from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("dudajustino/SimplexRegression")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SimplexRegression)
## Fit simplex regression model
model <- simplexreg(y ~ x1 + x2, data = your_data, 
                    link.mu = "plogit1", 
                    link.sigma2 = "log")

# Model summary
summary(model)
```

## Code of Conduct

Please note that the SimplexRegression project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
