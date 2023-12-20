#testthat for lowerCor

R <- lowerCor(ability)
r <- cor(ability,use="pairwise")
test_that("correlations from lowerCor match cor with missing data", 
    expect_equivalent(R,r))  
    
    
    