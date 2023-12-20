#testthat for tetrachoric
#Alsto tests that comorbidity and phi work
#data from Kirk, 1973

Q1 = .5
Q2 = .5
t0 <- comorbidity(.5,.2275,.11375)  # 0
t4 = comorbidity(.5,.5,.3154950)  # .4
t5= comorbidity(.5,.5,.3333333) # .5
t8 <- comorbidity(.5,.5,.39758450) #.8

test_that("tetrachoric matchs Kirk", 
    expect_equivalent(t0$tetra$rho,0)) 
     
 test_that("tetrachoric matchs Kirk", 
    expect_equivalent(t4$tetra$rho,.4,tolerance=10^-4))   
test_that("tetrachoric matchs Kirk", 
    expect_equivalent(t5$tetra$rho,.5,tolerance=10^-5))
test_that("tetrachoric matchs Kirk", 
    expect_equivalent(t8$tetra$rho,.8,tolerance=10^-5))    
    
