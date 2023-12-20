#Tests for SplitHalf

#First for an even set of items
#with an obvious split

R1 <- matrix(.36,3,3)
R2 <- matrix(.25,3,3)
g <- matrix(.16,6,6)
R12 <- superMatrix(R1,R2)
R <- g + R12
diag(R ) <- 1
dimnames(R) <-list( paste0("V",1:6),paste0("V",1:6))
sp <- splitHalf(R)

worst.cov <-  sum(R[c(1,2,3),c(4,5,6)])
worst.v1 <-   sum(R[c(1,2,3), c(1,2,3)])
worst.v2 <-   sum(R[c(4,5,6), c(4,5,6)])
worst.r <- worst.cov/sqrt(worst.v1*worst.v2)
split.min  <-  2 * worst.r /(1+worst.r )

best.cov <-    sum(R[c(1,4,5),c(2,3,6)])
best.v1 <-   sum(R[ c(1,4,5), c(1,4,5)])
best.v2 <-   sum(R[c(2,3,6), c(2,3,6)])
best.r <- best.cov/sqrt(best.v1*best.v2)
split.max <- best.r*2/(1+best.r)

test_that("Max split half ",
expect_equal(sp$maxrb,split.max,tolerance=.001))

test_that("Min split half",
expect_equal(sp$minrb,split.min,tolerance=.001))
test_that("average r )",
expect_equal(mean(R[lower.tri(R)]),sp$av.r))

#now find alpha from variances
 k <- 6 
alpha <- (sum(R) -k)*(k/(k-1))/sum(R)
test_that("alpha from variances", 
expect_equal(alpha,sp$alpha))

# see https://www.r-bloggers.com/2019/11/automated-testing-with-testthat-in-practice/ for help