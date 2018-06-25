context("designBased")

outcome_vec = c(1, 0, 0, 1, 0)
pcf_cat_vec = c('A', 'B', 'B', 'A', "C")
inst_id_vec = c(1,1,2,2,2)
gamma = 0.5
designBased(0, outcome_vec, pcf_cat_vec, inst_id_vec)

test_that('Design Based Scoring Works', {
  expect_error( designBased(gamma, c(1,1,1,1,2), pcf_cat_vec, inst_id_vec))
  expect_equal( designBased( 0, outcome_vec, pcf_cat_vec, inst_id_vec)$Z, c(NaN, NaN) )
  expect_equal( designBased(1, outcome_vec, pcf_cat_vec, inst_id_vec)$Z, c(0,0))
})
