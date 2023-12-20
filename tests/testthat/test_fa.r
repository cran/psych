#Tests for fa
#compare with published solutions
#and also compare to factanal
#Various output from Harman 

#Harman.Burt page 206
#minres solution to a singular matrix

minres <- data.frame(F1= c(.982,.935,.833,.720,.676,.526,.515,.355),
                     F2 = c(.065,-.111,-.550,-.091,.330, .164,.583,-.128),
                     h2 = c(.968, .888,.997,.526,.566,.304,.605,.142))
                     
f2 <- fa(Harman.Burt,2, rotate="none",warnings=FALSE)
#throws warnings about a bad matrix
diff <- f2$loadings - minres[,1:2]
test_that("Minresidual solution to Harman.Burt",
    expect_equal(max(abs(diff)),0,tolerance=.001))

#communalities
h2.diff <- f2$communalities -  rowSums(minres[,1:2]^2)  #find the communalities based upon the tabled factor loadings

test_that("Minres solution to Harman.Burt h2",
  expect_equal(max(abs(h2.diff)),0,tolerance=.002))
  
h2.diff1 <- f2$communalities - minres[,3]          #use the communalities from the book
 test_that("Minres solution to Harman.Burt h2  from formula",
  expect_equal(max(abs(h2.diff1)),0,tolerance=.0025)) 
####

#compare with factanal 
f3 <- fa(Thurstone,3, fm="mle",rotate="varimax")
f <- factanal(covmat=Thurstone,factors=3,rotation="varimax")
test_that("fa loadings match those from factanal", 
    expect_equivalent(f3$loadings,f$loading,tolerance=10^-5))  
#further tests are done by testing omega which also involves factoring and rotations

