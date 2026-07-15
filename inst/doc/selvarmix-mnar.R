## ----setup, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)
options(width = 78)

## ----data-------------------------------------------------------------------
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
dim(x)

## ----complete-fit-----------------------------------------------------------
set.seed(2026)
fit <- SelvarClustLasso(
  x = x,
  nbcluster = 3,
  models = "VVI",
  rank = seq_len(ncol(x)),
  hsize = 2,
  nbcores = 1,
  workflow = "imputed_sruw",
  init_control = list(method = "mclust_hc", seed = 2026),
  selection_control = list(stopping = "consecutive", seed = 2026)
)

## ----missing-data, eval=FALSE-----------------------------------------------
# set.seed(2027)
# x_miss <- x
# rho <- c(0.02, 0.08, 0.15)
# mask_probability <- matrix(
#   rho[latent_class], nrow = nrow(x_miss), ncol = ncol(x_miss)
# )
# mask <- matrix(runif(length(x_miss)), nrow = nrow(x_miss)) < mask_probability
# x_miss[mask] <- NA_real_
# 
# common <- list(
#   nbcluster = 2:4,
#   models = "VVI",
#   num_vals_penalty = 3,
#   hsize = 2,
#   nbcores = 1,
#   init_control = list(method = "mclust_hc", seed = 2027),
#   selection_control = list(stopping = "consecutive", seed = 2027)
# )
# 
# fit_decoupled <- do.call(
#   SelvarClustLasso,
#   c(list(x = x_miss, workflow = "decoupled_mnarz"), common)
# )
# 
# fit_joint <- do.call(
#   SelvarClustLasso,
#   c(list(x = x_miss, workflow = "joint_mnarz"), common)
# )

## ----warm-controls, eval=FALSE----------------------------------------------
# fit_warm <- SelvarClustLasso(
#   x,
#   nbcluster = 2:4,
#   workflow = "imputed_sruw",
#   nbcores = 1,
#   rank_control = list(warm_start = "both")
# )

## ----methods, results='hide', fig.keep='none'-------------------------------
fit
summary(fit)
plot(fit)

