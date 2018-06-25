context("Inference")

#Initialize parameters
n = 100
iters = 1000
p = 2
prior_var = 9

#Generate sample data
b_vec = rnorm(p, mean = 0, sd = 5)
X = matrix(rnorm(n * p), ncol = p)
p_vec = pnorm(X %*% b_vec)
prior_var_vec = rep(prior_var, p)

#Sample data
y = rbinom(n, 1,
           prob =p_vec  )
#Inference
mcmc_iters = probitFit(y, X, prior_var_vec,
                       iters = iters)
mcmc_iters_k = probitFitk(y, X, prior_var_vec,
                          iters = iters)

beta_hat = apply(mcmc_iters, 2, median)
beta_hat_k = apply(mcmc_iters_k, 2, median)

test_that("Approximatly Daniel K's results", {
  expect_equal(beta_hat, beta_hat_k, tolerance = 0.25)
})

test_that('Check outputs', {
  expect_equal(sum(!is.finite(mcmc_iters)),
               0) #All output elements finite
  expect_gt(min(apply(mcmc_iters, 2, var)),
            0) #Parameters should have positive posterior variance
})


##Make sure it works with sparce matrices
p = 10
n = 100
y = rbinom(n, 1, 0.5)
X = Matrix::Matrix( rbinom(n * p, 1, 0.5), ncol = p)
prior_var_vec = rep(10, p)

test_that('Sparce Matrix works', {
  expect_error(probitFit(y, X, prior_var_vec), NA)
})


##Test bayesian regression
test_that('Linear Regression Works', {
  expect_error(bayesianRegression(y, X, prior_var_vec, iters = iters), NA)
})
