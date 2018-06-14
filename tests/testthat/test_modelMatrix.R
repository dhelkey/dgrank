context("modelMatrix")

##Show how modelMatrix handles variables
df_full = data.frame(outcome = c(1,1,0,0,0,1,0),
                    id = c(1, 1, 3, 2, 5, 2, 3), #7 observations
                     color = c('red', 'red', 'blue', 'blue', 'red', 'red', 'red'),
                     tooHot = c(0, 0, 0, 1, 0, 1, 0),
                     favL = c('L', 'A', 'A', 'Z', 'A', 'D', 'D'),
                     w = c(1,1,1,1,1,1,1),
                weird = c(NaN, NA, NA, 1, Inf, NaN, 1 ))
df = df_full[-7] #Version w/o weird values
n = dim(df)[1] #Number of observations

test_that("Handles edge cases and weird values", {
  expect_equal(modelMatrix(NULL), NULL)
  expect_error(modelMatrix(df_full))
})
#

test_that('Works on a variety of variable types',{
  expect_equal( c(modelMatrix(df$outcome, sparse = FALSE)),
                df$outcome)
  expect_equal( c(modelMatrix(df$id, sparse = FALSE)),
               df$id)
  expect_equal( c(modelMatrix(df$color, sparse = FALSE)),
                 c(1,1,0,0,1,1,1))
})

#Categoricals behave as expected
hot_mat = modelMatrix(df['tooHot'])
favL_mat = modelMatrix(df['favL'])
test_that("Categorical tests", {
  expect_equal(dim(hot_mat), c(7,1))
  expect_equal(dim(favL_mat), c(7,3))
  expect_equal(colnames(hot_mat),'tooHot')
  expect_equal(colnames(favL_mat),
                c('favLD','favLL','favLZ'))
  expect_equal(as.vector(hot_mat[ ,1]),
               df$tooHot)
  expect_equal(as.vector(favL_mat[,1]),
              c(0,0,0,0,0,1,1))
})

test_that("class", {
  expect_is(modelMatrix(df,  sparse = FALSE),
            'matrix')
  expect_is(modelMatrix(df, sparse = TRUE),
            'dgCMatrix')
})

test_that("intercept", {
  expect_equal(as.vector(modelMatrix(df, intercept=TRUE)[,1]),
               rep(1,n))
})
