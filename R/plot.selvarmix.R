.selvarmix_missingness_probabilities <- function(object) {
  parameters <- object$parameters
  engine_parameters <- if (is.list(object$parametersMNARz)) {
    object$parametersMNARz$parameters
  } else {
    NULL
  }
  rho <- .selvarmix_value_or(parameters$rho, engine_parameters$rho)
  if (!is.null(rho) && length(rho)) {
    rho <- as.matrix(rho)
    if (nrow(rho) != object$nbcluster) {
      stop(
        "Stored missingness probabilities do not match nbcluster.",
        call. = FALSE
      )
    }
    if (ncol(rho) > 1L &&
        max(abs(sweep(rho, 1L, rho[, 1L], "-"))) > 1e-12) {
      stop(
        "Class-only MNARz requires one missingness probability per component.",
        call. = FALSE
      )
    }
    return(as.numeric(rho[, 1L]))
  }

  alpha <- .selvarmix_value_or(parameters$alpha, engine_parameters$alpha)
  if (is.null(alpha) || !length(alpha)) return(NULL)
  alpha <- as.matrix(alpha)
  if (nrow(alpha) != object$nbcluster) {
    stop("Stored missingness parameters do not match nbcluster.",
         call. = FALSE)
  }
  probabilities <- stats::pnorm(alpha)
  if (ncol(probabilities) == 1L ||
      max(abs(sweep(
        probabilities,
        1L,
        probabilities[, 1L],
        "-"
      ))) <= 1e-12) {
    return(as.numeric(probabilities[, 1L]))
  }
  stop(
    "Class-only MNARz requires one missingness probability per component.",
    call. = FALSE
  )
}

#' Plot a SelvarMixMNAR fit
#'
#' @param x A fitted object inheriting from class `selvarmix`.
#' @param type Display variable-role counts, posterior classification
#'   probabilities, or fitted class-specific missingness probabilities.
#' @param main,xlab,ylab Optional axis and title labels.
#' @param col Optional plotting colors.
#' @param ... Additional graphical parameters passed to the base graphics
#'   function used for the selected display.
#' @return The plotted numerical values, invisibly.
#' @export
plot.selvarmix <- function(
  x,
  type = c("roles", "classification", "missingness"),
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  col = NULL,
  ...
) {
  .validate_selvarmix_result(x, strict = FALSE)
  type <- match.arg(type)

  if (identical(type, "roles")) {
    counts <- c(
      S = length(x$S),
      U = length(x$U),
      W = length(x$W),
      `R within S` = length(x$R)
    )
    if (is.null(main)) main <- "Variable roles (R is a subset of S)"
    if (is.null(ylab)) ylab <- "Number of variables"
    if (is.null(xlab)) xlab <- "Role"
    if (is.null(col)) {
      col <- grDevices::hcl.colors(4L, palette = "Dark 3")
    }
    graphics::barplot(
      counts,
      main = main,
      xlab = xlab,
      ylab = ylab,
      col = col,
      ...
    )
    plotted <- data.frame(
      role = names(counts),
      count = as.integer(counts),
      primaryPartition = names(counts) != "R within S",
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    return(invisible(plotted))
  }

  if (identical(type, "classification")) {
    probabilities <- as.matrix(x$proba)
    if (!is.numeric(probabilities) || !nrow(probabilities) ||
        !ncol(probabilities)) {
      stop("Posterior probabilities are unavailable for this result.",
           call. = FALSE)
    }
    ordering <- order(x$partition, seq_along(x$partition))
    probabilities <- probabilities[ordering, , drop = FALSE]
    if (is.null(main)) {
      main <- "Posterior probabilities ordered by assigned component"
    }
    if (is.null(xlab)) xlab <- "Ordered observation"
    if (is.null(ylab)) ylab <- "Posterior probability"
    if (is.null(col)) {
      col <- grDevices::hcl.colors(ncol(probabilities), palette = "Dark 3")
    }
    graphics::matplot(
      seq_len(nrow(probabilities)),
      probabilities,
      type = "l",
      lty = 1L,
      ylim = c(0, 1),
      main = main,
      xlab = xlab,
      ylab = ylab,
      col = col,
      ...
    )
    if (ncol(probabilities) <= 12L) {
      graphics::legend(
        "topright",
        legend = paste("Component", seq_len(ncol(probabilities))),
        col = col,
        lty = 1L,
        bty = "n"
      )
    }
    attr(probabilities, "observationOrder") <- ordering
    return(invisible(probabilities))
  }

  probabilities <- .selvarmix_missingness_probabilities(x)
  if (is.null(probabilities) || !length(probabilities) ||
      any(!is.finite(probabilities))) {
    stop(
      "Class-specific missingness probabilities are unavailable for this result.",
      call. = FALSE
    )
  }
  if (is.null(main)) main <- "Class-specific missingness probabilities"
  if (is.null(ylab)) ylab <- "Missingness probability"
  if (is.null(xlab)) xlab <- "Component"

  probabilities <- as.numeric(probabilities)
  names(probabilities) <- paste("Component", seq_along(probabilities))
  if (is.null(col)) {
    col <- grDevices::hcl.colors(length(probabilities), palette = "Dark 3")
  }
  graphics::barplot(
    probabilities,
    ylim = c(0, 1),
    main = main,
    xlab = xlab,
    ylab = ylab,
    col = col,
    ...
  )
  invisible(probabilities)
}

#' @export
plot.selvarmix_collection <- function(x, which = 1L, ...) {
  .validate_selvarmix_collection(x)
  if (is.character(which)) {
    if (length(which) != 1L || is.na(which) || !(which %in% names(x))) {
      stop("which must name one result in the collection.", call. = FALSE)
    }
    index <- match(which, names(x))
  } else {
    if (!is.numeric(which) || length(which) != 1L || !is.finite(which) ||
        which != round(which) || which < 1L || which > length(x)) {
      stop("which must select one result in the collection.", call. = FALSE)
    }
    index <- as.integer(which)
  }
  plot(x[[index]], ...)
}
