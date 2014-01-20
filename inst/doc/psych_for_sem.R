### R code from vignette source 'psych_for_sem.Rnw'

###################################################
### code chunk number 1: psych_for_sem.Rnw:128-136
###################################################
library(psych)
set.seed(42)
tau <- sim.congeneric(loads=c(.8,.8,.8,.8)) #population values
tau.samp <- sim.congeneric(loads=c(.8,.8,.8,.8),N=100) # sample correlation matrix for 100 cases
round(tau.samp,2)  
tau.samp <- sim.congeneric(loads=c(.8,.8,.8,.8),N=100, short=FALSE) 
tau.samp
dim(tau.samp$observed)


###################################################
### code chunk number 2: psych_for_sem.Rnw:145-147
###################################################
cong <- sim.congeneric(N=100)
round(cong,2)


###################################################
### code chunk number 3: psych_for_sem.Rnw:153-155
###################################################
#plot.new()
m1 <- structure.diagram(c("a","b","c","d"))


###################################################
### code chunk number 4: psych_for_sem.Rnw:170-177
###################################################
set.seed(42)
gload=matrix(c(.9,.8,.7),nrow=3)
fload <- matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,
rep(0,9),.7,.6,.4),   ncol=3)
fload #echo it to see the structureSw
bifact <- sim.hierarchical(gload=gload,fload=fload)
round(bifact,2)


###################################################
### code chunk number 5: psych_for_sem.Rnw:184-188
###################################################
op <- par(mfrow=c(1,2))
m.bi <- omega(bifact,title="A bifactor model")
m.hi <- omega(bifact,sl=FALSE,title="A hierarchical model")
op <- par(mfrow = c(1,1))


###################################################
### code chunk number 6: psych_for_sem.Rnw:206-209
###################################################
circ <- sim.circ(16)
f2 <- fa(circ,2)
plot(f2,title="16 simulated variables in a circumplex pattern")


###################################################
### code chunk number 7: psych_for_sem.Rnw:226-230
###################################################
set.seed(42)
fx <- c(.9,.8,.7,.6)
cong1 <- sim.structure(fx)
cong1


###################################################
### code chunk number 8: psych_for_sem.Rnw:236-240
###################################################
set.seed(42)
fx  <- matrix(c(.9,.8,.7,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)
three.fact <- sim.structure(fx)
three.fact


###################################################
### code chunk number 9: psych_for_sem.Rnw:245-247
###################################################
#plot.new()
three.fact.mod <-structure.diagram(fx)


###################################################
### code chunk number 10: psych_for_sem.Rnw:256-261
###################################################
Phi = matrix(c(1,.5,.3,.5,1,.2,.3,.2,1), ncol=3)
cor.f3 <- sim.structure(fx,Phi)
fx
Phi
cor.f3


###################################################
### code chunk number 11: psych_for_sem.Rnw:268-272
###################################################
fxs <- structure.list(9,list(F1=c(1,2,3),F2=c(4,5,6),F3=c(7,8,9)))
Phis <- phi.list(3,list(F1=c(2,3),F2=c(1,3),F3=c(1,2)))
fxs  #show the matrix
Phis #show this one as well


###################################################
### code chunk number 12: psych_for_sem.Rnw:278-280
###################################################
#plot.new()
corf3.mod <- structure.diagram(fxs,Phis)


###################################################
### code chunk number 13: psych_for_sem.Rnw:291-292
###################################################
f3.p <- Promax(fa(cor.f3$model,3))


###################################################
### code chunk number 14: psych_for_sem.Rnw:297-299
###################################################
#plot.new()
mod.f3p <- structure.diagram(f3.p,cut=.2)


###################################################
### code chunk number 15: psych_for_sem.Rnw:311-318
###################################################
set.seed(42)
fx  <- matrix(c(.9,.8,.7,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)
fy <- c(.6,.5,.4)
Phi <- matrix(c(1,.48,.32,.4,.48,1,.32,.3,.32,.32,1,.2,.4,.3,.2,1), ncol=4)
twelveV <- sim.structure(fx,Phi, fy)$model
colnames(twelveV) <-rownames(twelveV) <- c(paste("x",1:9,sep=""),paste("y",1:3,sep=""))
round(twelveV,2)


###################################################
### code chunk number 16: psych_for_sem.Rnw:324-327
###################################################
fxs <- structure.list(9,list(X1=c(1,2,3), X2 =c(4,5,6),X3 = c(7,8,9)))
phi <- phi.list(4,list(F1=c(4),F2=c(4),F3=c(4),F4=c(1,2,3)))
fyx <- structure.list(3,list(Y=c(1,2,3)),"Y")


###################################################
### code chunk number 17: psych_for_sem.Rnw:332-334
###################################################
#plot.new()
sg3 <- structure.diagram(fxs,phi,fyx)


###################################################
### code chunk number 18: psych_for_sem.Rnw:347-354
###################################################
fxh  <- structure.list(9,list(X1=c(1:3),X2=c(4:6),X3=c(7:9),g=NULL))
fy <- structure.list(3,list(Y=c(1,2,3)))
Phi <- diag(1,5,5)
Phi[4,c(1:3)] <- letters[1:3]
Phi[5,4] <- "r"
#plot.new()
hi.mod <-structure.diagram(fxh,Phi, fy)


###################################################
### code chunk number 19: psych_for_sem.Rnw:370-372
###################################################
principal(cong1$model)
fa(cong1$model)


###################################################
### code chunk number 20: psych_for_sem.Rnw:380-387
###################################################
pc3 <- principal(bifact,3)
pa3 <- fa(bifact,3,fm="pa")
ml3 <- fa(bifact,3,fm="ml")
pc3
pa3
ml3
factor.congruence(list(pc3,pa3,ml3))


###################################################
### code chunk number 21: psych_for_sem.Rnw:393-395
###################################################
ml3p <- Promax(ml3)
ml3p


###################################################
### code chunk number 22: psych_for_sem.Rnw:404-405
###################################################
om.bi <- omega(bifact)


###################################################
### code chunk number 23: psych_for_sem.Rnw:417-418
###################################################
om.hi <- omega(bifact,sl=FALSE)


###################################################
### code chunk number 24: psych_for_sem.Rnw:426-427
###################################################
om.bi


###################################################
### code chunk number 25: psych_for_sem.Rnw:435-436
###################################################
ic <- ICLUST(bifact,title="Hierarchical cluster analysis of bifact data")


###################################################
### code chunk number 26: psych_for_sem.Rnw:458-463
###################################################
library(sem)
 mod.tau <- structure.sem(c("a","a","a","a"),labels=paste("V",1:4,sep=""))
mod.tau   #show it
sem.tau <- sem(mod.tau,cong,100)
summary(sem.tau,digits=2)


###################################################
### code chunk number 27: psych_for_sem.Rnw:468-473
###################################################
mod.cong <- structure.sem(c("a","b","c","d"),labels=paste("V",1:4,sep=""))
mod.cong  #show the model
sem.cong <- sem(mod.cong,cong,100)
summary(sem.cong,digits=2)
anova(sem.cong,sem.tau) #test the difference between the two models


###################################################
### code chunk number 28: psych_for_sem.Rnw:484-489
###################################################
mod.one <- structure.sem(letters[1:9],labels=paste("V",1:9,sep=""))
mod.one  #show the model
sem.one <- sem(mod.one,bifact,100)
summary(sem.one,digits=2)
round(residuals(sem.one),2)


###################################################
### code chunk number 29: psych_for_sem.Rnw:495-499
###################################################
f1 <- factanal(covmat=bifact,factors=1)
mod.f1 <- structure.sem(f1)
sem.f1 <- sem(mod.f1,bifact,100)
sem.f1


###################################################
### code chunk number 30: psych_for_sem.Rnw:508-513
###################################################
f3 <-fa(bifact,3,rotate="varimax")
mod.f3 <- structure.sem(f3)
sem.f3 <- sem(mod.f3,bifact,100)
summary(sem.f3,digits=2)
round(residuals(sem.f3),2)


###################################################
### code chunk number 31: psych_for_sem.Rnw:522-525
###################################################
f3 <-fa(bifact,3)     #extract three factors and do an oblique rotation
mod.f3 <- structure.sem(f3) #create the sem model
mod.f3   #show it 


###################################################
### code chunk number 32: psych_for_sem.Rnw:533-535
###################################################
#plot.new()
mod.f3 <- structure.diagram(f3)


###################################################
### code chunk number 33: psych_for_sem.Rnw:552-553
###################################################
om.th.bi  <- omega(Thurstone,digits=3)


###################################################
### code chunk number 34: psych_for_sem.Rnw:562-565
###################################################
sem.bi <- sem(om.th.bi$model,Thurstone,213) #use the model created by omega
summary(sem.bi,digits=2)
stdCoef(sem.bi,digits=2)


###################################################
### code chunk number 35: psych_for_sem.Rnw:578-579
###################################################
om.hi <- omega(Thurstone,digits=3,sl=FALSE)


###################################################
### code chunk number 36: psych_for_sem.Rnw:587-591
###################################################
sem.hi <- sem(om.hi$model,Thurstone,213)
summary(sem.hi,digits=2)
stdCoef(sem.hi,digits=2)
anova(sem.hi,sem.bi)


###################################################
### code chunk number 37: psych_for_sem.Rnw:600-605
###################################################
om.holz.bi <- omega(Holzinger,4)
sem.holz.bi <- sem(om.holz.bi$model,Holzinger,355)
om.holz.hi <- omega(Holzinger,4,sl=FALSE)
sem.holz.hi <- sem(om.holz.hi$model,Holzinger,355)
anova(sem.holz.bi,sem.holz.hi)


###################################################
### code chunk number 38: psych_for_sem.Rnw:615-616
###################################################
om.sem <- omegaSem(Thurstone,n.obs=213)


###################################################
### code chunk number 39: psych_for_sem.Rnw:624-625
###################################################
omega.diagram(om.sem,sort=FALSE)


###################################################
### code chunk number 40: psych_for_sem.Rnw:640-641
###################################################
sessionInfo()


