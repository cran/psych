#tests for various correlation functions
#corPlot

R <- corPlot(attitude)
r <- cor(attitude)
test_that("correlations from CorPlot match cor", 
    expect_equivalent(R,r))
    
R <- corPlot(ability)
r <- cor(ability,use="pairwise")
test_that("correlations from CorPlot match cor with missing data", 
    expect_equivalent(R,r))  
    
R <- lowerCor(ability,show=FALSE)
r <- cor(ability,use="pairwise")
test_that("correlations from lowerCor match cor with missing data", 
    expect_equivalent(R,r))  