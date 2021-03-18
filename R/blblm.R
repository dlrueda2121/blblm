#' @aliases blblm-package
#' @import purrr
#' @import stats
#' @import furrr
#' @import future
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))



#' Builds blblm object
#'
#' Creates linear regression objects with little bag of bootstraps
#'
#' @param formula formula
#' @param data dataframe
#' @param m integer
#' @param B integer
#' @param workers integer
#' @param fastlm boolean
#'
#' @return object class blblm
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, fastlm = TRUE)
#' @export
blblm <- function(formula, data, m = 10, B = 5000, workers = 1, fastlm = TRUE) {
  if(workers > 1){
    suppressWarnings(future::plan(multiprocess, workers = workers))
    options(future.rng.onMisuse = "ignore")
    data_list <- split_data(data, m)
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, fastlm = fastlm))
    res <- list(estimates = estimates, formula = formula)

  }else{
    data_list <- split_data(data, m)
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, fastlm = fastlm))
    res <- list(estimates = estimates, formula = formula)

  }
  class(res) <- "blblm"
  invisible(res)

}

#' Splits data
#'
#' Splits data into m parts of approximated equal sizes
#'
#' @param data dataframe
#' @param m integer
#'
#' @return list
#' @examples
#' split_data(data = mtcars, m = 3)
#' @export
split_data <- function(data, m) {
  set.seed(141)
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' Computes estimates
#'
#' Computes the regression model estimates, and replicate lm1 process
#'
#' @param formula regression formula
#' @param data dataframe
#' @param n integer
#' @param B integer
#' @param fastlm boolean
#'
#' @return ‘pdb’ with replicated atomic coordinates
#' @examples
#' lm_each_subsample(mpg ~ wt * hp, data = mtcars, n = 10, B = 100, fastlm = TRUE)
#' @export
lm_each_subsample <- function(formula, data, n, B, fastlm) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n, fastlm), simplify = FALSE)
}

#' Computes estimates blb
#'
#' Computes the regression estimates for a blb dataset
#'
#' @param X matrix
#' @param y vector
#' @param n integer
#' @param fastlm boolean
#'
#' @return list
#' @examples
#' n <- 7 ; p <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#' lm1(X, y, 10, fastlm = TRUE)
#' @export
lm1 <- function(X, y, n, fastlm) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  if(fastlm){
    fit <- RcppEigen::fastLm(X, y, weigths = freqs)
    list(coef = blbcoef(fit), sigma = blbsigma(fit))
  }else{
    fit <- lm.wfit(X, y, w = freqs)
    list(coef = blbcoef(fit), sigma = blbsigma(fit))
  }
}

#' BLB Coefficients
#'
#' Computes the coefficients from model fit
#'
#'
#' @param fit fitted blblm model
#'
#' @return vector
#' @export
#' @examples
#' blbcoef(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
blbcoef <- function(fit) {
  coef(fit)
}

#' Sigma
#'
#' Computes sigma from fit
#'
#' @param fit blblm object
#'
#' @return double
#' @examples
#' blbsigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
#' @export
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
