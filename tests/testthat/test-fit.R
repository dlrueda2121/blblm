test_that("fit multiworkers w/fastlm  works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)

  expect_that(fit, is_a("blblm"))
})

test_that("fit single thread w/lm works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 1, fastlm = FALSE)

  expect_that(fit, is_a("blblm"))
})

test_that("coef(fit) fastlm/multiworker works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)
  coefs <- coef(fit)

  expect_equal(length(coefs), 3)
})

test_that("sigma(fit) fastlm/multiworker works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)
  sigma = sigma(fit)

  expect_that(sigma, is_a("numeric"))
})

test_that("confint(fit, c('wt','hp')) fastlm/multiworker works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)
  cint <- confint(fit)

  expect_equal(dim(cint), c(2,2))
})

test_that("sigma(fit, confidence = TRUE) fastlm/multiworker works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)
  sigma = sigma(fit, confidence = TRUE)

  expect_equal(length(sigma), 3)
})

test_that("predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170))) fastlm/multiworker works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)
  prediction =  predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))

  expect_equal(length(prediction), 2)
})

test_that("predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE) fastlm/multiworker works", {

  fit <- blblm(mpg ~ wt + hp, data = mtcars, m = 3, B = 100, workers = 2, fastlm = TRUE)
  prediction =  predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)

  expect_equal(dim(prediction), c(2,3))
})




