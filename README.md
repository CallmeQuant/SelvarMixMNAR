# SelvarMixMNAR

`SelvarMixMNAR` provides Gaussian model-based clustering with variable selection for complete and incomplete data, including a class-dependent MNARz missingness model. It implements and extends the methods developed in the NeurIPS 2025 paper *A Unified Framework for Variable Selection in Model-Based Clustering with Missing Not at Random* by Binh H. Ho, Long Nguyen Chi, TrungTin Nguyen, Binh T. Nguyen, Van Ha Hoang, and Christopher Drovandi ([arXiv preprint](https://arxiv.org/abs/2505.19093)).

The package combines three components:

- a regularized Gaussian-mixture procedure that ranks variables using sparse component means and precision matrices;
- the SRUW variable-role model, which partitions variables into clustering variables (`S`), redundant variables (`U`) explained by a subset `R` of `S`, and independent variables (`W`);
- a parsimonious MNARz model in which each latent component has one missingness probability shared across the modeled variables.

The fitted missingness model is class-only MNARz. The package does not fit MNARy mechanisms or the variable-specific `MNARz_j` model. In the decoupled workflow, `mnarz_control$mecha = "mixed"` accepts a fixed logical mask `is_mnar`. Missingness indicators for coordinates marked `TRUE` contribute to the MNARz likelihood and share the component-specific probability `rho_k`; coordinates marked `FALSE` are treated as ignorable. Thus, this option provides a class-only sensitivity analysis rather than an `MNARz_j` model.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "CallmeQuant/SelvarMixMNAR",
  dependencies = TRUE
)

library(SelvarMixMNAR)
```

Using `dependencies = TRUE` installs the packages required by the optional missing-data and covariance-distance workflows. The required `gcimputeR` revision is pinned in `DESCRIPTION`. Because the package contains C++ source code, Windows users need the version of Rtools corresponding to their R installation.

The equivalent `devtools` command is

```r
devtools::install_github(
  "CallmeQuant/SelvarMixMNAR",
  dependencies = TRUE
)
```

A local source tree can be installed from the directory containing `DESCRIPTION`:

```r
remotes::install_local(
  "path/to/SelvarMixMNAR",
  dependencies = TRUE,
  build_vignettes = FALSE
)
```

For incomplete data, the ranking stage first constructs a completed working matrix. The default imputation method uses `gcimputeR`; set `use_copula = FALSE` to use `missRanger`. Non-Euclidean covariance distances require the optional `shapes` package. The Euclidean/Frobenius distance is implemented internally.

## Quick start

The following example contains two clustering variables, two redundant variables, and two independent noise variables. The simulated class labels are used only for external evaluation and are not supplied to the estimator.

```r
library(SelvarMixMNAR)

set.seed(2025)
n <- 90
latent_class <- rep(1:3, each = n / 3)

signal_1 <- rnorm(n, c(-2, 0, 2)[latent_class], 0.7)
signal_2 <- rnorm(n, c(1.5, -1.5, 0)[latent_class], 0.7)

x <- cbind(
  signal_1 = signal_1,
  signal_2 = signal_2,
  redundant_1 = 0.8 * signal_1 + rnorm(n, sd = 0.35),
  redundant_2 = -0.7 * signal_2 + rnorm(n, sd = 0.35),
  noise_1 = rnorm(n),
  noise_2 = rnorm(n)
)

set.seed(2026)
fit <- SelvarClustLasso(
  x = x,
  nbcluster = 2:4,
  models = "VVI",
  num_vals_penalty = 3,
  hsize = 2,
  nbcores = 1,
  workflow = "imputed_sruw"
)

fit
summary(fit)
plot(fit)
```

The fitted object contains the disjoint role sets `S`, `U`, and `W`; the regression subset `R` of `S`; the selected number of components; posterior component-membership probabilities; cluster assignments; fitted parameters; the completed data when applicable; and numerical diagnostics.

## Incomplete-data workflows

For incomplete input, `workflow = "auto"` selects `decoupled_mnarz`. Three workflows are available:

| Workflow | Description |
|---|---|
| `imputed_sruw` | Completes the data, then performs variable ranking and SRUW role selection without modeling the missingness process. |
| `decoupled_mnarz` | Selects SRUW roles on the completed data and then fits an MNARz mixture separately to the original incomplete data. This is the default for incomplete input. |
| `joint_mnarz` | Uses a completed matrix for ranking, then jointly scores the SRUW role model and the class-only MNARz model on the original incomplete data. Forward and reverse role updates are compared by observed-data BIC. |

Use `decoupled_mnarz` for a less expensive two-stage analysis. Use `joint_mnarz` when variable-role decisions should be based on the same observed-data criterion as the MNARz fit.

```r
set.seed(2025)
x_miss <- x
rho <- c(0.02, 0.08, 0.15)

mask_probability <- matrix(
  rho[latent_class],
  nrow = nrow(x_miss),
  ncol = ncol(x_miss)
)
mask <- matrix(
  runif(length(x_miss)),
  nrow = nrow(x_miss)
) < mask_probability
x_miss[mask] <- NA_real_

fit_decoupled <- SelvarClustLasso(
  x = x_miss,
  nbcluster = 2:4,
  workflow = "decoupled_mnarz",
  models = "VVI",
  num_vals_penalty = 3,
  hsize = 2,
  nbcores = 1
)

fit_joint <- SelvarClustLasso(
  x = x_miss,
  nbcluster = 2:4,
  workflow = "joint_mnarz",
  models = "VVI",
  num_vals_penalty = 3,
  hsize = 2,
  nbcores = 1
)
```

## Initialization and penalty paths

Available initialization methods include hierarchical clustering from `mclust`, `mclust` EM, k-means, random starts, a user-supplied state, a previous fit, and deterministic multistart initialization. Hierarchical `mclust` initialization is the default. Alternative initializations are useful for sensitivity analysis because Gaussian-mixture likelihoods are non-convex.

Penalty-grid points are fitted independently by default. Set `rank_control$warm_start` to `"inner"`, `"outer"`, or `"both"` to enable continuation:

- `"inner"` warm-starts the graphical-lasso solver between EM iterations;
- `"outer"` carries the complete fitted state along a deterministic penalty path;
- `"both"` applies both forms of continuation.

Points on the same continuation path are evaluated serially.

A grid fit is excluded from the ranking score if it does not converge, has a non-finite objective, produces a non-positive-definite covariance matrix, or decreases the objective. Excluded grid points remain available in the diagnostics and are not treated as zero-activity fits. A candidate number of components is scored only when the fraction of valid grid fits is at least `min_scorable_fraction`.

## Scope and limitations

- Input variables must be numeric, and the clustering model is Gaussian.
- The MNARz model is class-only: all modeled coordinates share one missingness probability within each latent component.
- Under `mnarz_control$mecha = "mixed"`, the `is_mnar` mask is fixed by the analyst; the package does not estimate which coordinates are MNAR.
- The default SRUW backend is `mclust` with covariance model `VVI`.

For the model factorization, interpretation of the workflows, and reproducible control settings, see

```r
vignette("selvarmix-mnar")
```

## Relationship to SelvarMix

The regularized ranking and SRUW role-selection procedures build on `SelvarMix` 1.2.1 by Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau. The SRUW formulation originates in the variable-role model of Maugis, Celeux, and Martin-Magniette, while the regularized ranking follows the later SelvarMix methodology of Celeux, Maugis-Rabusseau, and Sedki.

`SelvarMixMNAR` retains the main clustering interface of `SelvarMix` and adds class-only MNARz estimation, incomplete-data workflows, numerical diagnostics, and initialization and penalty-path metadata.

## Citation and license

Use

```r
citation("SelvarMixMNAR")
```

to obtain the recommended citation. `SelvarMixMNAR` is distributed under GPL (>= 3).