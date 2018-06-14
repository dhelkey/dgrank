context('helperfunctions')


##Make sure Cholesky Solver works
y = ChickWeight$weight
X = cbind(1, ChickWeight$Time)
XtX = t(X) %*% X
Xty = t(X) %*% y
beta_ls = solve(XtX, Xty) #Naieve least squares estimate
beta_chol = linCholSolver( chol(XtX), Xty)

test_that('cholesky solver works', {
	expect_equal(beta_ls, beta_chol)
})