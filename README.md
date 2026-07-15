# SelvarMixMNAR

`SelvarMixMNAR` fits Gaussian model-based clustering with variable selection
when observations may be missing not at random. It is the direct R
implementation accompanying the NeurIPS 2025 paper *A Unified Framework for
Variable Selection in Model-Based Clustering with Missing Not at Random* by
Binh H. Ho, Long Nguyen Chi, TrungTin Nguyen, Binh T. Nguyen, Van Ha Hoang,
and Christopher Drovandi ([see arxiv version here](https://arxiv.org/abs/2505.19093)).

The package combines three statistical components:

- a penalized Gaussian-mixture ranking based on sparse component means and
  precision matrices;
- an SRUW decomposition, in which variables are clustering-relevant (`S`),
  redundant (`U`) through a regression on `R` contained in `S`, or independent
  noise (`W`);
- a parsimonious MNARz model with one missingness probability for each latent
  component, shared across variables.

Under class-only MNARz, each latent component has one missingness probability
shared across variables. The fitting interface is restricted to this
parsimonious mechanism; MNARy and variable-specific MNARzj require different
likelihoods. For sensitivity analyses in the separate decoupled fit,
`mnarz_control$mecha = "mixed"` accepts a user-specified logical `is_mnar`
mask. Coordinates marked `TRUE` share the same class-specific probability
`rho_k`; mask factors for coordinates marked `FALSE` are treated as ignorable
and omitted from the modeled likelihood. The analyst fixes the mask, and all
selected coordinates share the same class effect; this is a class-only
sensitivity model rather than MNARzj.

## Installation

Install the development version directly from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "CallmeQuant/SelvarMixMNAR",
  dependencies = TRUE
)

library(SelvarMixMNAR)
```

Setting `dependencies = TRUE` installs the optional packages used by the
missing-data and covariance-distance workflows. It also resolves
`gcimputeR` at the commit pinned in `DESCRIPTION`. Installation compiles the
package's C++ source, so Windows users need the version of Rtools that matches
their R installation.

Users who already have `devtools` installed may equivalently run
`devtools::install_github("CallmeQuant/SelvarMixMNAR", dependencies = TRUE)`.

### Installation from a local checkout

Package developers can install a local checkout by supplying the directory
that contains `DESCRIPTION`:

```r
remotes::install_local(
  "path/to/SelvarMixMNAR",
  dependencies = TRUE,
  build_vignettes = FALSE
)
```

For release validation, build from the parent directory and install the
resulting source archive:

```sh
R CMD build SelvarMixMNAR
R CMD INSTALL SelvarMixMNAR_0.1.0.9000.tar.gz
```

Incomplete-data fits require an imputation method for the ranking stage. The
default uses `gcimputeR`; `missRanger` is available through
`use_copula = FALSE`. Non-Euclidean covariance distances require the optional
`shapes` package. The default Euclidean/Frobenius distance is implemented
internally.

## A first example

The following example has two clustering variables, two redundant
variables, and two independent noise variables. `latent_class` is retained
only for external evaluation and is not passed to the estimator.

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

The result stores the disjoint role partition `S`/`U`/`W`, the regression
predictor subset `R` contained in `S`, the selected component count, posterior
probabilities, partition, fitted parameters, imputed data when applicable,
and numerical diagnostics.

## Missing-data workflows

For an incomplete input, `workflow = "auto"` resolves to the decoupled
workflow. The available targets are:

| Workflow | Statistical operation |
|---|---|
| `imputed_sruw` | Rank variables and fit the SRUW model on a completed working data set without a missingness likelihood. |
| `decoupled_mnarz` | Select SRUW roles on the completed working data, then fit a separate role-agnostic MNARz mixture. This is the default for incomplete input. |
| `joint_mnarz` | Rank on a completed matrix, then fit the SRUW and class-only MNARz terms jointly to the original incomplete matrix. Forward and reverse role updates are scored by total observed-data BIC. |

The decoupled estimator separates imputed-data role selection from MNARz
estimation and is computationally less expensive. The joint estimator couples
the SRUW and missingness terms within one observed-data objective. The
appropriate choice depends on whether computational economy or joint modeling
of variable roles and missingness is the primary consideration.

```r
set.seed(2025)
x_miss <- x
rho <- c(0.02, 0.08, 0.15)
mask_probability <- matrix(
  rho[latent_class], nrow = nrow(x_miss), ncol = ncol(x_miss)
)
mask <- matrix(runif(length(x_miss)), nrow = nrow(x_miss)) < mask_probability
x_miss[mask] <- NA_real_

fit_decoupled <- SelvarClustLasso(
  x_miss,
  nbcluster = 2:4,
  workflow = "decoupled_mnarz",
  models = "VVI",
  num_vals_penalty = 3,
  hsize = 2,
  nbcores = 1
)

fit_joint <- SelvarClustLasso(
  x_miss,
  nbcluster = 2:4,
  workflow = "joint_mnarz",
  models = "VVI",
  num_vals_penalty = 3,
  hsize = 2,
  nbcores = 1
)
```

## Initialization, penalty paths, and failures

Gaussian initializers are constructed in R and passed to native routines as a
common numeric state. Available methods include hierarchical Mclust,
Mclust EM, k-means, random starts, user-supplied state, a previous fit, and
deterministic multistart. Hierarchical Mclust is the data-driven default;
alternative initializers support sensitivity analysis when the mixture
likelihood has several well-separated local maxima.

Penalty-grid fits are cold by default. `rank_control$warm_start` may be
`"inner"`, `"outer"`, or `"both"`. Inner continuation warm-starts the
graphical-lasso solver between EM iterations. Outer continuation transfers the
complete fitted state along a deterministic serial penalty path. Dependent
points on one warm path are never parallelized.

Nonconverged fits and fits with a non-finite objective, a non-positive-definite
covariance, or an objective decrease are omitted from the ranking
score. Their grid positions remain recorded rather than being replaced by
zero activity. A component-count candidate is scored only when the retained
grid fraction reaches `min_scorable_fraction`.

## Statistical scope

- The input variables must be numeric, and the clustering model is
  Gaussian.
- MNARz assumes that missingness depends on latent component membership
  through one component-specific probability shared across coordinates.
- The optional mixed mask is specified by the analyst for the separate
  decoupled fit, and estimation conditions on that fixed partition.
- `joint_mnarz` is typically more expensive than the decoupled workflow.
- Mclust and Rmixmod are the supported SRUW backends. MixAll tokens stop before
  native execution because an STK failure may terminate the R session. The
  default is Mclust with model `VVI`.

See the package vignette, `vignette("selvarmix-mnar")`, for the model
factorization, workflow interpretation, and reproducible control settings.

## Relationship to SelvarMix

The penalized ranking and SRUW role-selection layer descends from SelvarMix
1.2.1 by Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau. This
package retains the principal clustering function names while adding
class-only MNARz workflows, numerical diagnostics, and reproducible
initialization and path metadata.

## Citation and license

Use `citation("SelvarMixMNAR")` to cite the accompanying paper. The package is
distributed under GPL (>= 3).
