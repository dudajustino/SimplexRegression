
# SimplexRegression: A Package for Simplex Regression with parametric and fixed mean link function

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R package and associated data for the scientific article:

“Simplex regression with a flexible logit link: inference and application to cross-country impunity data” by Justino, M. E. C. and Cribari-Neto, F.

## 📑 Table of Contents

- [🎯 Overview](#-Overview)
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

## 🎯 Overview

Simplex regression is a powerful statistical framework for modeling bounded continuous responses in (0,1), such as proportions, rates and indices.

Traditional approaches use **fixed mean link functions** (logit, probit, log-log, complementar log-log, cauchit). This package extends these models by introducing **parametric mean link functions** (plogit1, plogit2), which include an additional parameter λ estimated from the data, providing greater flexibility to the model.

### Why Parametric Mean Link Functions?

- ✅ **Data-driven flexibility**: The link parameter λ is estimated from the data, not imposed
- ✅ **Captures asymmetry**: plogit1 and plogit2 accommodate different directions of asymmetry
- ✅ **Nests standard links**: When λ = 1, plogit1 and plogit2 reduce to the logit link
- ✅ **Testable specification**: Formal score tests evaluate whether standard links are adequate
- ✅ **Better predictive performance**: Often outperforms fixed link specifications in practice

---

## 🌟 Key Features

### Parametric Mean Link Functions
- **plogit1**: `g(μ, λ) = log((1-μ)^(-λ) - 1)`
- **plogit2**: `g(μ, λ) = log(μ^λ / (1-μ^λ))`
- **Data-driven selection**: Choose between plogit1 and plogit2 using model selection criteria

### Fixed Mean Link Functions
- `logit`, `probit`, `cloglog`, `loglog`, `cauchit`

### Dispersion Modeling
- Model heterogeneity with covariates in the dispersion submodel
- Logarithmic, square root, or identity dispersion links

### Comprehensive Diagnostics
- **Residual analysis**: Quantile, standardized weighted, deviance, and bias-corrected residuals
- **Visual diagnostics**: Half-normal plots with simulated envelopes, Q-Q plots, worm plots
- **Influence measures**: Cook's distance, leverage (hat values), and local influence analysis

### Model Selection Tools
- **Scout Score (SS)** criterion with optional penalty for parametric links
- **Penalized information criteria**: AIC^(λ), BIC^(λ), HQIC^(λ)
- **Score tests**: Test H₀: λ = 1 (logit link) vs. H₁: λ ≠ 1

### Global and Local Influence Analysis
- **Case-weight perturbation**: Identify observations with high leverage
- **Response perturbation**: Assess sensitivity to changes in response values
- **Curvature-based measures**: Detect jointly influential observations

---

## 📦 Installation

### Development version from GitHub

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install SimplexRegression
devtools::install_github("dudajustino/SimplexRegression")
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
