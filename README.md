
# SimplexRegression: A Package for Simplex Regression with parametric and fixed mean link function

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R package and associated data for the article:

“Simplex regression with a flexible logit link: inference and application to cross-country impunity data” by Justino, M. E. C. and Cribari-Neto, F.

## 📑 Table of Contents

- [🎯 Overview](#-Overview)
- [🌟Key Features](#-key-features)
- [📂 Repository Structure](repository-structure)
- [🎯 Code of Conduct](#-code-of-conduct)
- [🛠️ Installation](#-installation)
- [🚀 Getting Started & Examples](#-getting-start)
- [📚 Vignettes](#-vignettes)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📖 Citation](#-citation)
- [📧 Contact](#-contact)

---

## 🎯 Overview

Simplex regression is a powerful statistical framework for modeling bounded continuous responses in (0,1), such as proportions, rates and indices.

Traditional approaches use **fixed mean link functions** (logit, probit, log-log, complementary log-log, cauchit). This package extends these models by introducing **parametric mean link functions** (plogit1, plogit2), which include an additional parameter λ estimated from the data, providing greater flexibility to the model.

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
- **Penalized information criteria**: AIC<sup>(λ)</sup>, BIC<sup>(λ)</sup>, HQIC<sup>(λ)</sup>
- **Score tests**: Test H₀: λ = 1 (logit link) vs. H₁: λ ≠ 1

### Global and Local Influence Analysis
- **Case-weight perturbation**: Identify influential observations using case-weight perturbation
- **Response perturbation**: Identify influential observations using response perturbation
- **Curvature-based measures**: Detect jointly influential observations

---

📂 Repository Structure

```
.
├── R/                  # Source code for all R functions.
├── data/               # Processed data included in the package (.rda).
├── man/                # R package documentation files for functions.
├── vignettes/          # Detailed tutorials and case studies (.Rmd).
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

--- 

## 🎯 Code of Conduct

Please note that the SimplexRegression project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

---

## 🛠 Installation

### Development version from GitHub

This research compendium can be installed as an R package directly from GitHub. This is the recommended method as it handles all dependencies automatically.

First, ensure you have the remotes package. If not, install it from CRAN:

``` r
if (!require("remotes")) {
  install.packages("remotes")
}
```

Then, install the package from GitHub:

``` r
remotes::install_github("dudajustino/SimplexRegression", 
                        dependencies = TRUE,
                        build_vignettes = TRUE)
```
Prerequisites

This package requires the following external R packages. You can run the command below to ensure all dependencies are installed on your system before proceeding.

``` r
install.packages(c("expint", "gamlss", "graphics", "stats", "moments", "tseries"))

```

Last Tested Environment The scripts were last successfully tested on:

R version: 4.4.1
Platform: x86_64-w64-mingw32 (64-bit)

---

## 🚀 Getting Started & Examples

```r
library(SimplexRegression)

# Load the dataset
data("impunity_dataset")

# Fit a simplex regression model with parametric mean link
model1 <- simplexreg(
  Impunity ~ EconomicFreedom + I(GDP^1.15) + I(HDI*HealthSpending) + I(Democracy*Press) + dnordic | 
             EconomicFreedom + HealthSpending + Democracy,
  data = impunity_dataset,
  link.mu = "plogit1",
  link.sigma2 = "log"
)

# View the results
summary(model1)

# Model diagnostics
plot(model1, which = 1:7, type = "quantile")

# Half-normal plot with simulated envelopes
hnp.simplexreg(model1, type = "sweighted2")

# Fit a simplex regression with fixed log-log mean link:

# Simplex model with logit link
model2 <- simplexreg(
  Impunity ~ EconomicFreedom + I(GDP^1.15) + I(HDI*HealthSpending) + I(Democracy*Press) + dnordic | 
             EconomicFreedom + HealthSpending + Democracy,
  data = impunity_dataset,
  link.mu = "loglog",
  link.sigma2 = "log"
)

summary(model2)
```

---

## 📚 Vignettes

Detailed tutorials and applications:

```r
# List all vignettes
browseVignettes("SimplexRegression")

# View specific vignette
vignette("impunity-analysis", package = "SimplexRegression")
```

### Available vignettes:

1. **Impunity Analysis**: Complete analysis of impunity data across 119 countries
   - Parametric vs fixed mean link functions
   - Model selection strategies
   - Comprehensive diagnostics
   - Influence analysis

---

## 🤝 Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📖 Citation

If you use this package in your research, please cite:

```bibtex
<<<<<<< HEAD
@article{Justino+Cribari_2025,
  title = {Simplex Regression with a Flexible Logit Link: Inference and Application to Cross-Country Impunity Data},
  author = {Maria Eduarda da Cruz Justino and Francisco Cribari-Neto},
=======
@Article{Justino+Cribari_2025,
  title = {Simplex regression with a flexible logit link: Inference and application to cross-country impunity data},
  author = {Maria E. C. Justino and Francisco Cribari-Neto},
>>>>>>> f41e0bec68c086d8869d4f20af4a8a02e42e1b42
  year = {2025},
  url = {https://github.com/dudajustino/SimplexRegressionAMM},
}
```

---

## 📧 Contact

**Maria Eduarda da Cruz Justino**
- 📫 Email: eueduardacruz@gmail.com
- 🐙 GitHub: [@dudajustino](https://github.com/dudajustino)

**Francisco Cribari-Neto**
- 🏛️ Departamento de Estatística, UFPE

