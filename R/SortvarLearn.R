SortvarLearn <- function(
    x,
    z,
    type = "lasso",
    lambda = seq(20, 100, by = 10),
    rho = seq(1, 2, length = 2),
    nbcores = min(2, detectCores(all.tests = FALSE, logical = FALSE))) {
  .supervised_route_abort()
}
