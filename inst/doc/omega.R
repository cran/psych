### R code from vignette source 'omega.Rnw'

###################################################
### code chunk number 1: omega.Rnw:162-166
###################################################
library(psych)    #make the psych package active 
library(psychTools) #make psychTools active
om <- omega(Thurstone)  #do the analysis 
om  #show it 


###################################################
### code chunk number 2: Thurstone
###################################################
png('Thurstone.png')
omega.diagram(om)
dev.off()


###################################################
### code chunk number 3: omega
###################################################
omega(Thurstone)


###################################################
### code chunk number 4: anxiety
###################################################

anxiety <- sai[c("anxious", "jittery", "nervous" ,"tense", "upset","at.ease" ,  "calm" , 
   "confident", "content","relaxed")]
describe(anxiety)
lowerCor(anxiety)
om <- omega(anxiety,2)  #specify a two factor solution
summary(om)   #summarize the output


###################################################
### code chunk number 5: anxietyplot
###################################################
png('anxiety.png')
omega.diagram(om, main="Omega analysis of  two factors of anxiety")
dev.off()


###################################################
### code chunk number 6: direct
###################################################
om <- omegaDirect(Thurstone)  
om


###################################################
### code chunk number 7: drawdirect
###################################################
png('direct.png')
omega.diagram(om, main="Direct Schmid Leihman solution")
dev.off()


###################################################
### code chunk number 8: omega.Rnw:734-737
###################################################
om <- omega(holzinger.swineford[8:31],4)  #the exploratory solution

omegaSem(holzinger.swineford[8:31],4) #the confirmatory solution


###################################################
### code chunk number 9: holzinger
###################################################



###################################################
### code chunk number 10: omega.Rnw:831-834
###################################################
jen <- sim.hierarchical()  #use the default values
om <- omega(jen)
om


###################################################
### code chunk number 11: jensen
###################################################
png('jensen.png' )
omega.diagram(om)
dev.off()


###################################################
### code chunk number 12: Simulate1
###################################################
fx <- matrix(c(.7,.6,.5,.7,.6,.5,.8,.7,.6, 
    .6,.6,.6,rep(0,9),c(.6,.5,.6),rep(0,9),.6,.6,.6),ncol=4)
 simx <-sim.structure(fx)

om <- omega(simx$model)
dsl <- omegaDirect(simx$model)


###################################################
### code chunk number 13: Simulate.2
###################################################
 lowerMat(simx$model)
summary(om)
summary(dsl)
fa.congruence(list(om,dsl,fx))


###################################################
### code chunk number 14: omega.Rnw:1014-1015
###################################################
sessionInfo()


