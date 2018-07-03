context('DG Functions')

p = 10
n = 1000

#Simulate coefficients
dat_mat = matrix( rnorm(p * n), ncol = p)
n_vec = c(50, 20, 1, 3, 5, 10, 20, 25, 100, 9)

dg_scores = dgConstraint(dat_mat, n_vec)
dg_z_scores = toDG(dat_mat, n_vec, baseline = FALSE)
baseline_z_scores = toDG(dat_mat, n_vec, baseline = TRUE)

test_that('D-G functions work as expected', {
  expect_equal( max(dg_scores %*% n_vec), 0)
  expect_equal( max(baseline_z_scores[ ,1]), 0)
  expect_equal( max(apply(baseline_z_scores, 2, var)), 1 )
})
