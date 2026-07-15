# SelvarMixMNAR 0.1.0.9000

## Statistical models

- The package now fits the class-only MNARz mechanism, with one missingness
  probability per latent component shared across variables.
- Incomplete inputs use `decoupled_mnarz` by default: SRUW roles are selected
  on an imputed working data set and a separate role-agnostic MNARz mixture is
  fitted afterward.
- `joint_mnarz` fits the SRUW and class-only MNARz terms under one
  observed-data criterion. Forward and reverse role updates are scored by
  total observed-data BIC.
- MNARy and variable-specific MNARzj are outside the fitted-model interface.

## Ranking and role selection

- Penalty grids are constructed separately for each candidate component
  count. Failed grid positions are retained but omitted from activity scores.
- The mean update solves the precision-coupled penalized block by coordinate
  descent and reports objective and KKT diagnostics.
- Graphical-lasso updates check solver termination, positive definiteness,
  covariance-precision consistency, and the relevant KKT conditions.
- Cold, inner, outer, and combined warm-start modes govern whether solver or
  complete mixture state is reused. Complete fitted state is transferred
  along each outer path.
- The SRUW stopping parameter implements a rolling count of consecutive
  nonpositive criterion increments. The former block rule remains available
  only as `stopping = "legacy_block"`.
- Empty-`R` candidates refit `U` together with `W` under the independent model
  instead of relabeling variables without refitting.

## Initialization and numerical diagnostics

- Hierarchical Mclust, Mclust EM, k-means, random, user-supplied,
  previous-fit, and deterministic-multistart initializers return a common
  numeric state with seed and objective metadata.
- Pure, mixed-mask, and joint MNARz routines use updated-mean covariance
  centering and coherent final likelihood, posterior, and imputation states.
- Tiny components, objective decreases, covariance adjustments,
  nonconvergence, and degenerate missingness patterns produce classed
  conditions or numerical diagnostics.
- Repeated missingness patterns are evaluated through cached Cholesky factors
  and streamed sufficient statistics in the joint fit.

## Interface

- Fitted models use a versioned `selvarmix` result structure with standard
  `print()`, `summary()`, and `plot()` methods.
- Results retain role sets, posterior probabilities, partitions, fitted and
  imputed data, preprocessing transformations, initialization provenance,
  penalty-grid status, selection diagnostics, and workflow-specific fit
  information.
- `true_labels` is ignored during estimation and may be used only after fitting
  to compute external evaluation measures.
- Likelihood ranking and supervised-lasso fitting raise classed
  unavailable-feature conditions.
- The unused supervised native estimator and `glasso` dependency have been
  removed.
- MixAll model specifications are rejected before backend execution. Mclust
  `VVI` is the supported default SRUW backend.
- Advanced diagnostics use descriptive names including `fit_metadata`,
  `component_mass_rule`, `mechanism_specification`, and
  `initialization_schema_version`.
