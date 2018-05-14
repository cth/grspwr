test_that("randomSwap permutes elements but retains its set of elements", {
  expect_false(isTRUE(all.equal(randomSwap(1:5),1:5)))
  expect_equal(sort(randomSwap(1:5)),1:5)
})

exampleGRS <- data.frame(SNP = c("rs1", "rs2", "rs3"),
                         Weight = c(0.01, 0.01, 0.2),
                         EAF = c(0.42, 0.42, 0.42),
                         INFO= c(0.9, 0.9, 0.2))


calculated_power <- grspwr(exampleGRS,n=100, alpha = 0.05)

test_that("Power is on the interval [0-1]", {
  expect_more_than(calculated_power$power, 0)
  expect_less_than(calculated_power$power, 1)
})
