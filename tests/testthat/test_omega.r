#Tests for omega and for factor analysis
#this is a test for sim.hierarchical, fa, and omega all in one

#Uses the data model from Jensen and Weng (1994)
#we simulate their data for a 3 lower order 1 higher factor using sim.hierarchical

sim.jen <- sim.hierarchical()  #the default is to simulate the Jensen-Weng model

om <- omega(sim.jen,plot=FALSE)   #do a factor analysis and then a higher level factoring 
loadings <- om$schmid$sl[,1:5]

#from Jensen and Weng
jensen <- data.frame(g= c(.72, .63, .54, .56, .48,.40,.42,.35,.28),
                    F1 = c(.3487, .3051,.2615, 0,0,0,0,0,0),
                    F2 = c(0,0,0,.42,.36,.30, 0,0,0),
                    F3 = c(0,0,0,0,0,0,.428,.3570,.2856),
                    h2= c(.64,.49,.36,.49,.36,.25,.36,.25, .16)
                    )
diff <- loadings - jensen

sum.diff <- sum(abs(diff))
test_that("omega output matches Jensen-Weng)",
     expect_equivalent(sum.diff,0,tolerance=.001))

