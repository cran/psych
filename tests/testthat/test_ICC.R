#Tests for ICC  -- does it match the numbers from Shrout and Fleiss
sf <- matrix(c(
9,    2,   5,    8,
6,    1,   3,    2,
8,    4,   6,    8,
7,    1,   2,    6,
10,   5,   6,    9,
6,   2,   4,    7),ncol=4,byrow=TRUE)
colnames(sf) <- paste("J",1:4,sep="")
rownames(sf) <- paste("S",1:6,sep="")
sf  #example from Shrout and Fleiss (1979)
temp <- ICC(sf,lmer=FALSE) 

test_that("Shrout and Fleiss ICC(1,1)",
     expect_equal(round(temp$results[1,2],2),.17))
test_that("Shrout and Fleiss ICC(2,1)",
expect_equal(round(temp$results[2,2],2),.29))
test_that("Shrout and Fleiss ICC(3,1)",
expect_equal(round(temp$results[3,2],2),.71))
#NOW TEST THE K CONDITION
test_that("Shrout and Fleiss ICC(1,4)",
     expect_equal(round(temp$results[4,2],2),.44))
test_that("Shrout and Fleiss ICC(2,4)",
expect_equal(round(temp$results[5,2],2),.62))
test_that("Shrout and Fleiss ICC(3,4)",
expect_equal(round(temp$results[6,2],2),.91))

# see https://www.r-bloggers.com/2019/11/automated-testing-with-testthat-in-practice/ for help