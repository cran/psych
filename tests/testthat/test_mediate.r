#test the mediate function against published results from Hayes ()

# A simple mediation example is the Tal_Or data set (pmi for Hayes)
#The pmi data set from Hayes is available as the Tal_Or data set. 
mod4 <- mediate(reaction ~ cond + (pmi), data =Tal_Or,n.iter=50) 
summary(mod4)


#traditional regression  (c')
b.coeff <-c(.5249,.2544, .5064)
t.coeff <-c(.9585,.9943,5.2183)
test_that("traditional regression  (c')",
  expect_equivalent( b.coeff, mod4$cprime.reg$beta[,1],tolerance=.002))
  test_that("traditional regression  (c')  t ",
   expect_equivalent( t.coeff, mod4$cprime.reg$t,tolerance=.001))

#total effect reaction (c ) 3.25   .4957
b.coeff <- c( 3.25,   .4957)
t.coeff <- c( 17.0525, 1.7860)
test_that("c path (direct effect-- intercept and slope)",
        expect_equivalent(mod4$total.reg$beta,b.coeff,tolerance=.001)) 
    test_that("c path (direct effect-- t values)",
        expect_equivalent(mod4$total.reg$t,t.coeff,tolerance=.001)) 
    
#a effect X on M 
b.coeff <- c( 5.3769 , .4765 )
t.coeff <- c( 33.2222, 2.0218)
test_that("a path (direct effect-- intercept and slope)",
        expect_equivalent(mod4$a.reg$beta,b.coeff,tolerance=.001)) 
 test_that("a path (direct effect-- t values)",
        expect_equivalent(mod4$a.reg$t,t.coeff,tolerance=.001))   
        
#b effect (M on Y) 
b.coeff <- c(.5064 )
t.coeff <- c(5.2185)
test_that("b path M on Y) -- slope)",
        expect_equivalent(mod4$b.reg$beta[-1],b.coeff,tolerance=.001)) 
 test_that("a path (direct effect-- t values)",
        expect_equivalent(mod4$b.reg$t[-1],t.coeff,tolerance=.001))   
       
#ab effect (the indirect effect)) 
b.coeff <- c(.2413 )

test_that("ab  on Y) -- slope)",
        expect_equivalent(mod4$ab,b.coeff,tolerance=.001)) 
 test_that("a path (direct effect-- t values)",
        expect_equivalent(mod4$b.reg$t[-1],t.coeff,tolerance=.001))   



################
#Two mediators (from Hayes model 6 (chapter 5))
mod6 <- mediate(reaction ~ cond + (pmi) + (import), data =Tal_Or,n.iter=50) 
summary(mod6)
#note this is the model on page 135 


#traditional regression  (c')
b.coeff <-c(-.1498, .1034,.3965, .3244 )
t.coeff <-c(-.28287, .4324, 4.2645, 4.5857)
test_that("traditional regression  (c')",
  expect_equivalent( b.coeff, mod6$cprime.reg$beta[,1],tolerance=.002))
  test_that("traditional regression  (c')  t ",
   expect_equivalent( t.coeff, mod6$cprime.reg$t,tolerance=.001))

#total effect reaction (c ) 
b.coeff <- c( 3.25,   .4957)
t.coeff <- c( 17.0525, 1.7860)
test_that("c path (direct effect-- intercept and slope)",
        expect_equivalent(mod6$total.reg$beta,b.coeff,tolerance=.001)) 
    test_that("c path (direct effect-- t values)",
        expect_equivalent(mod6$total.reg$t,t.coeff,tolerance=.001)) 
    
#a effect X on M 
b.coeff <- c(5.3769,.47653,3.9077, .6268)
t.coeff <- c(33.222,2.0218, 18.3704, 2.0234)
test_that("a path (direct effect-- intercept and slope)",
        expect_equivalent(mod6$a.reg$beta,b.coeff,tolerance=.001)) 
 test_that("a path (direct effect-- t values)",
        expect_equivalent(mod6$a.reg$t,t.coeff,tolerance=.001))   
        
#b effect (M on Y) 
b.coeff <- c(.5064 )
t.coeff <- c(5.2185)
test_that("b path M on Y) -- slope)",
        expect_equivalent(mod4$b.reg$beta[-1],b.coeff,tolerance=.001)) 
 test_that("a path (direct effect-- t values)",
        expect_equivalent(mod4$b.reg$t[-1],t.coeff,tolerance=.001))   
       
#ab effect (the indirect effect)) 
b.coeff <- c(.3923, .1890, .2033 )

test_that("total ab  on Y) -- slope)",
        expect_equivalent(mod6$ab,b.coeff[1],tolerance=.001)) 
        
test_that("each ab  on Y) -- slope)",
        expect_equivalent(mod6$all.ab,b.coeff[2:3],tolerance=.001))









#Moderated mediation is done for the Garcia (Garcia, 2010) data set.
# (see Hayes, 2013 for the protest data set
#n.iter set to 50 (instead of default of 5000) for speed of example
#no mediation, just an interaction
#note that we do not zero center to be compatible with Hayes p 226
mod7 <- mediate(liking ~  sexism * prot2 , zero=FALSE,data=Garcia, n.iter = 50, plot=FALSE)

beta=c(7.7062, -.4725, -3.7727,.8336)
t.coeff = c(7.3750, -2.3184, -3.0084, 3.4224)
test_that("Moderated regression with slopes and intercepts of",
    expect_equivalent(mod7$total.reg$beta,beta, tolerance=.001))
  
test_that("Moderated regression with t of",
    expect_equivalent(mod7$total.reg$t,t.coeff, tolerance=.001))







data(GSBE)   #The Garcia et al data set (aka GSBE)
mod11.4 <- mediate(liking ~  sexism * prot2 + (respappr), data=Garcia,
        n.iter = 50,zero=FALSE)   #to match Hayes
summary(mod11.4)
#to see this interaction graphically, run the examples in ?Garcia


