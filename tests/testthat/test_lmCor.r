#Tests for lmCor  -- does it match the output from lm?
#also a test that the graphics engine works.


mod <- (lm(rating ~ complaints + privileges, data = attitude))
mod1 <- lmCor(rating ~ complaints + privileges, data = attitude, std=FALSE, plot=FALSE) #do not standardize

test_that("lm and lmCor agree)",
     expect_equivalent(mod$coefficients,mod1$coefficients))

z.attitude <- data.frame(scale(attitude))  #standardize the data before doing lm
mod.z<-lm(rating ~ complaints + privileges, data = z.attitude)  #regressions on z scores
mod.z1 <- lmCor(rating ~ complaints + privileges, data = attitude, plot=FALSE)  #by default we standardize and 

test_that("lm and lmCor agree for standardized scores)",
     expect_equivalent(mod.z$coefficients,mod.z1$coefficients))


# the results are the same as the standardized lm


R <- cor(attitude) #find the correlations
#Do the regression on the correlations  
#Note that these match the regressions on the standard scores of the data
mod.r <- lmCor(rating ~ complaints + privileges, data =R, n.obs=30, plot=FALSE)


test_that("lmCor with raw data and  and lmCor for correlation matrices agree",
     expect_equivalent(mod.r$coefficients,mod.z1$coefficients[-1]))   #no intercept


#now, partial out learning and critical
mod.par <- lmCor(rating ~ complaints + privileges - learning - critical, data =R, n.obs=30)
#compare with the full regression:
mod.full <- lmCor(rating ~ complaints + privileges + learning + critical, data =R, n.obs=30)
test_that("lmCor with partial r match lmCor with full R agree",
     expect_equivalent(mod.full$coefficients[1:2],mod.par$coefficients))   #no intercept

