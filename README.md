
# SimplexRegression: A Package for Simplex Regression with parametric and fixed mean link function

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R package and associated data for the scientific article:

“Simplex regression with a flexible logit link: inference and application to cross-country impunity data” by Justino, M. E. C. and Cribari-Neto, F.

## 📑 Table of Contents

- [🎯 Overview](#-Overview)
- [🌟Key Features](#-key-features)
- [📂 Repository Structure](repository-structure)
- [🛠️ Installation](#-installation)
- [🚀 Getting Started & Examples](#-getting-start)
- [📚 Vignettes](#-vignettes)
- [🎯 Code of Conduct](#-code-of-conduct)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📖 Citation](#-citation)
- [📧 Contact](#-contact)

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

📂 Repository Structure

.
├── R/                  # Source code for all R functions.
├── data/               # Processed data included in the package (.rda).
├── man/                # R package documentation files for functions.
├── vignettes/          # Detailed tutorials and case studies (.Rmd).
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
├── LICENSE             # MIT License file.
└── README.md           # This file.

--- 

## 🛠 Installation

### Development version from GitHub

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install SimplexRegression
devtools::install_github("dudajustino/SimplexRegression")
```

## 🚀 Getting Started & Examples

```r
library(SimplexRegression)

# Load example dataset
data("impunity_dataset")

# Fit a simplex regression model with parametric mean link
model1 <- simplexreg(
  Impunity ~ EconomicFreedom + I(GDP^1.15) + I(HDI*HealthSpending) + I(Democracy*Press) + dnordic | 
             EconomicFreedom + HealthSpending + Democracy,
  data = impunity_dataset,
  link.mu = "plogit1",
  link.sigma2 = "log"
)

# View results
summary(model1)

# Model diagnostics
plot(model1, which = 1:7, type = "quantile")

# Half-normal plot with simulated envelope
hnp.simplexreg(model1, type = "sweighted2")

# Fit a simplex regression with fixed log-log mean link:

# Simple model with logit link
model2 <- simplexreg(
  Impunity ~ EconomicFreedom + I(GDP^1.15) + I(HDI*HealthSpending) + I(Democracy*Press) + dnordic | 
             EconomicFreedom + HealthSpending + Democracy,
  data = impunity_dataset,
  link.mu = "loglog",
  link.sigma2 = "log"
)

summary(model2)
```

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

## 🎯 Code of Conduct

Please note that the SimplexRegression project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

---

## 🤝 Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📖 Citation

If you use this package in your research, please cite:

```r
citation("SimplexRegression")
```

```bibtex
@Article{Justino+Cribari_2025,
  title = {Simplex Regression with a Flexible Logit Link: Inference and Application to Cross-Country Impunity Data},
  author = {Maria Eduarda da Cruz Justino and Francisco Cribari-Neto},
  year = {2025},
  url = {https://github.com/dudajustino/SimplexRegression},
}
```

---

## 📧 Contact

**Maria Eduarda da Cruz Justino**
- 📫 Email: eueduardacruz@gmail.com
- 🐙 GitHub: [@dudajustino](https://github.com/dudajustino)

**Francisco Cribari-Neto**
- 🏛️ Departamento de Estatística, UFPE

