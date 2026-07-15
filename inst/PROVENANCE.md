# Source and data provenance

This file records the origin of material distributed with SelvarMixMNAR. It
does not replace the copyright holders' licensing decisions.

## Software lineage

SelvarMixMNAR is derived from SelvarMix 1.2.1 by Mohammed Sedki, Gilles
Celeux, and Cathy Maugis-Rabusseau. The local reference copy was extracted
from the archived CRAN package, whose `DESCRIPTION` declares `GPL (>= 3)`:

<https://cran.r-project.org/src/contrib/Archive/SelvarMix/SelvarMix_1.2.1.tar.gz>

The current source retains the original exported names while adding MNARz
estimation, explicit numerical diagnostics, deterministic initialization and
penalty paths, and the joint and decoupled workflows described in the package
documentation. The MNARz formulation was checked against the local NeurIPS
2025 manuscript and the accompanying reference R programs supplied in
`Clustering-MNAR-Sportisse/`. Those reference programs are not included in the
package archive, and no license statement was found in their local directory.

The package does not vendor a third-party numerical library. Rcpp,
RcppArmadillo, glassoFast, Rmixmod, mclust, matrixStats, and igraph are linked
or loaded as declared dependencies. The headers under `inst/include/` are part
of the SelvarMix-derived native interface, not copies of those dependencies.
MixAll is not an installation dependency: its historical model tokens are
recognized only to return the documented error before its namespace or native
backend can be entered.

`inst/include/selectRegGen.hpp` is inherited from SelvarMix 1.2.1. Before its
encoding was normalized, it was byte-identical to the archived file (SHA-256
`8A96C5C17581FDD3BBC304B683A50672839048881AC8DA0263E7F1B152C8BBDA`).
The UTF-8-normalized file has SHA-256
`00DB34FF6A43B4C0DCAF6E8AFF153510DE190EB425FB5F2595C04421ABE4B13B`.
Its original GPL and Lucent permission text is retained verbatim. The scope of
the Lucent notice and the apparent absence of its referenced warranty
disclaimer require copyright-holder review before a tagged public release.

The standard Windows/Linux workflow follows the current
[`r-lib/actions` v2 examples](https://github.com/r-lib/actions/tree/v2/examples),
including explicit Pandoc and TinyTeX setup. Those examples are released under
CC0. The sanitizer workflow follows the
[R-hub container guidance](https://r-hub.github.io/containers/gha.html) and
uses the published Clang ASAN/UBSAN image by immutable digest. These repository
files are excluded from the R source archive by `.Rbuildignore`. The sanitizer
job verifies that the pinned R configuration contains
`-fsanitize=address,undefined` before dependencies are installed or checks are
run.

Before a public release, the copyright holders must approve the attribution
and licensing of the post-1.2.1 modifications, and a current maintainer must
approve the contact recorded in `DESCRIPTION`. The present metadata preserves
the authors and maintainer listed by SelvarMix 1.2.1; it does not assign
authorship for later contributions.

## Bundled data

The package archive contains no data objects. The legacy `wine.rda` file was
byte-identical to SelvarMix 1.2.1 (SHA-256
`F019E6E2C287C4C8BBDD0C8815509BE59DD52B8DAD6766CD484A5165AE29E41C`),
but its help page gave neither a primary source nor data-specific reuse terms.
It and its help page were therefore removed from the package source. The stale
`scenarioCor` help topic was also removed because no corresponding data object
was present. Historical copies remain only in the private development archive;
they cannot enter the public repository or source archive without independent
provenance and license evidence.

## Dependency evidence

`provenance/dependencies.csv` records the versions and licenses observed in
the audited Windows library. `gcimputeR` is pinned to the exact installed
[GitHub commit](https://github.com/udellgroup/gcimputeR/tree/0d7c0e6fd3be9ff7d1c3a3f8100d91e721dad468)
in both that manifest and `DESCRIPTION`. CRAN dependencies remain
version-resolved by the installation environment; the observed versions are
evidence for the audited configuration, not asserted minimum versions.

`provenance/ci.csv` records the intended operating systems, R versions,
dependency scopes, check arguments, and sanitizer image digest. Every row is
marked `configured_not_run`; only a completed GitHub Actions log can change
that status.

`provenance/artifacts.csv` gives a machine-readable summary of shipped source
and data categories. License compatibility has not been independently
adjudicated. In particular, `gcimputeR` declares GPL-2 while this package
declares GPL (>= 3); no `gcimputeR` source is bundled here, but the final
distributor remains responsible for confirming the dependency-license
position.

`provenance/repository_files.csv` records the SHA-256 and size of every
publishable repository file at the audited snapshot, together with its
byte-level relationship to the extracted SelvarMix 1.2.1 tree. The manifest
excludes itself to avoid a recursive hash and excludes ignored IDE and compiled
artifacts.
