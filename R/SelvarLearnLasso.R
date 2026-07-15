SelvarLearnLasso <- function(
    x,
    z,
    lambda = seq(20, 100, by = 10),
    rho = seq(1, 2, length = 2),
    type = "lasso",
    rank,
    hsize = 3,
    models = mixmodGaussianModel(listModels = c(
      "Gaussian_pk_L_C", "Gaussian_pk_Lk_C",
      "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck"
    )),
    rmodel = c("LI", "LB", "LC"),
    imodel = c("LI", "LB"),
    xtest = x,
    ztest = z,
    nbcores = min(2, detectCores(all.tests = FALSE, logical = FALSE))) {
  .supervised_route_abort()
}
