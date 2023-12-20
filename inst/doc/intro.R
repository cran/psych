### R code from vignette source 'intro.Rnw'

###################################################
### code chunk number 1: intro.Rnw:483-484
###################################################
    options(width=100)


###################################################
### code chunk number 2: intro.Rnw:488-492
###################################################
library(psych)   #need to make psych active the first time you call it
library(psychTools)  #additional tools and data are here
data(sat.act) 
describe(sat.act)  #basic descriptive statistics


###################################################
### code chunk number 3: intro.Rnw:499-501
###################################################
 #basic descriptive statistics by a grouping variable.
describeBy(sat.act ~ gender,skew=FALSE,ranges=FALSE)


###################################################
### code chunk number 4: intro.Rnw:509-512
###################################################
sa.mat <-  describeBy(sat.act ~ gender + education,
 skew=FALSE,ranges=FALSE,mat=TRUE)
headTail(sa.mat)


###################################################
### code chunk number 5: outlier
###################################################
png( 'outlier.png' )
d2 <- outlier(sat.act,cex=.8)
dev.off()


###################################################
### code chunk number 6: intro.Rnw:545-549
###################################################
x <- matrix(1:120,ncol=10,byrow=TRUE)
colnames(x) <- paste('V',1:10,sep='')
new.x <- scrub(x,3:5,min=c(30,40,50),max=70,isvalue=45,newvalue=NA)
new.x


###################################################
### code chunk number 7: intro.Rnw:566-578
###################################################
x <- matrix(1:40,ncol=10,byrow=TRUE)
y <- matrix(1:20,ncol=4)
xy <- vJoin(x,y)
xy

XY <-  vJoin(x,y,cnames=FALSE)
XY
#match on ids and columns
x <- bfi[1:5,1:10]
y <- bfi[3:8,2:6] 
xy <- vJoin(x,y)
xy #the merged data 


###################################################
### code chunk number 8: pairspanels
###################################################
png( 'pairspanels.png' )
sat.d2 <- data.frame(sat.act,d2) #combine the d2 statistics from before with the sat.act data.frame
pairs.panels(sat.d2,bg=c("yellow","blue")[(d2 > 25)+1],pch=21,stars=TRUE)
dev.off()


###################################################
### code chunk number 9: affect
###################################################
png('affect.png')
pairs.panels(affect[14:17],bg=c("red","black","white","blue")[affect$Film],pch=21,
     main="Affect varies by movies ")
dev.off()


###################################################
### code chunk number 10: affect1
###################################################
keys <- list(
EA = c("active", "energetic", "vigorous", "wakeful", "wide.awake", "full.of.pep",
       "lively", "-sleepy", "-tired", "-drowsy"),
 TA =c("intense", "jittery", "fearful", "tense", "clutched.up", "-quiet", "-still", 
       "-placid", "-calm", "-at.rest") ,
PA =c("active", "excited", "strong", "inspired", "determined", "attentive", 
          "interested", "enthusiastic", "proud", "alert"),
NAf =c("jittery", "nervous", "scared", "afraid", "guilty", "ashamed", "distressed",  
         "upset", "hostile", "irritable" )) 
scores <- scoreItems(keys,psychTools::msq[,1:75])
#png('msq.png')
# pairs.panels(scores$scores,smoother=TRUE,
#  main ="Density distributions of four measures of affect" )

#dev.off()


###################################################
### code chunk number 11: violin
###################################################
png('violin.png')
data(sat.act)
violinBy(SATV+SATQ ~ gender, data=sat.act,grp.name=cs(Verbal.M,Verbal.F, Quan.M,Quant.F), 
     main="Density Plot by gender for SAT V and Q")
dev.off()


###################################################
### code chunk number 12: intro.Rnw:702-704
###################################################
data(epi.bfi)
error.bars.by(epi.bfi[,6:10],epi.bfi$epilie<4)


###################################################
### code chunk number 13: intro.Rnw:717-719
###################################################
error.bars.by(sat.act[5:6],sat.act$gender,bars=TRUE,
        labels=c("Male","Female"),ylab="SAT score",xlab="")


###################################################
### code chunk number 14: intro.Rnw:733-737
###################################################
T <- with(sat.act,table(gender,education))
rownames(T) <- c("M","F")
error.bars.tab(T,way="both",ylab="Proportion of Education Level",xlab="Level of Education",
main="Proportion of sample by education level")


###################################################
### code chunk number 15: intro.Rnw:756-766
###################################################
op <- par(mfrow=c(1,2))
  data(affect)
colors <- c("black","red","white","blue")
 films <- c("Sad","Horror","Neutral","Happy")
affect.stats <- errorCircles("EA2","TA2",data=affect[-c(1,20)],group="Film",labels=films,
xlab="Energetic Arousal",   ylab="Tense Arousal",ylim=c(10,22),xlim=c(8,20),pch=16,
cex=2,colors=colors, main =' Movies effect on arousal')
 errorCircles("PA2","NA2",data=affect.stats,labels=films,xlab="Positive Affect",
  ylab="Negative Affect", pch=16,cex=2,colors=colors,  main ="Movies effect on affect")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 16: bibars
###################################################
data(bfi) 
png( 'bibars.png' )
bi.bars(bfi,"age","gender",ylab="Age",main="Age by males and females")
dev.off()


###################################################
### code chunk number 17: histo
###################################################
png('histo.png')
densityBy(bfi,"age",grp="gender")
dev.off()


###################################################
### code chunk number 18: scatterhist
###################################################
data(GERAS)

png( 'scatterHist.png' )
psych::scatterHist(F ~ M + gender, data=GERAS.scales, cex.point=.3,smooth=FALSE, 
xlab="Masculine Scale",ylab="Feminine Scale",correl=FALSE, 
d.arrow=TRUE,col=c("red","blue"), bg=c("red","blue"), lwd=4, title="Combined  M and F 
scales",cex.cor=2,cex.arrow=1.25)
dev.off()


###################################################
### code chunk number 19: intro.Rnw:843-844
###################################################
lowerCor(sat.act)


###################################################
### code chunk number 20: intro.Rnw:851-857
###################################################
female <- subset(sat.act,sat.act$gender==2)
 male <- subset(sat.act,sat.act$gender==1)
lower <- lowerCor(male[-1])
upper <- lowerCor(female[-1])
both <- lowerUpper(lower,upper)
round(both,2)


###################################################
### code chunk number 21: intro.Rnw:863-865
###################################################
diffs <-  lowerUpper(lower,upper,diff=TRUE)
round(diffs,2)


###################################################
### code chunk number 22: corplot.png
###################################################
png('corplot.png')
corPlot(Thurstone,numbers=TRUE,upper=FALSE,diag=FALSE,cex=.7,
     main="9 cognitive variables from Thurstone")
dev.off()


###################################################
### code chunk number 23: circplot.png
###################################################
png('circplot.png')
circ <- sim.circ(24)
r.circ <- cor(circ)
corPlot(r.circ,main='24 variables in a circumplex')
dev.off()


###################################################
### code chunk number 24: spider.png
###################################################
png('spider.png')
op<- par(mfrow=c(2,2))
spider(y=c(1,6,12,18),x=1:24,data=r.circ,fill=TRUE,main="Spider plot of 24 circumplex variables")
op <- par(mfrow=c(1,1))
dev.off()


###################################################
### code chunk number 25: intro.Rnw:933-934
###################################################
corTest(sat.act)


###################################################
### code chunk number 26: intro.Rnw:945-946
###################################################
r.test(50,.3)


###################################################
### code chunk number 27: intro.Rnw:952-953
###################################################
r.test(30,.4,.6)


###################################################
### code chunk number 28: intro.Rnw:960-961
###################################################
r.test(103,.4,.5,.1)


###################################################
### code chunk number 29: intro.Rnw:967-968
###################################################
r.test(103,.5,.6,.7,.5,.5,.8)  #steiger Case B 


###################################################
### code chunk number 30: intro.Rnw:976-977
###################################################
cortest(sat.act)


###################################################
### code chunk number 31: intro.Rnw:991-991
###################################################



###################################################
### code chunk number 32: tetrar.png
###################################################
png('tetrar.png')
draw.tetra()
dev.off()


###################################################
### code chunk number 33: intro.Rnw:1067-1068
###################################################
lmCor(y = 5:9,x=1:4,data=Thurstone)


###################################################
### code chunk number 34: intro.Rnw:1075-1077
###################################################
sc <- lmCor(y = 5:9,x=3:4,data=Thurstone,z=1:2)
round(sc$residual,2)


###################################################
### code chunk number 35: intro.Rnw:1090-1104
###################################################
#data from Preacher and Hayes (2004)
sobel <- structure(list(SATIS = c(-0.59, 1.3, 0.02, 0.01, 0.79, -0.35, 
-0.03, 1.75, -0.8, -1.2, -1.27, 0.7, -1.59, 0.68, -0.39, 1.33, 
-1.59, 1.34, 0.1, 0.05, 0.66, 0.56, 0.85, 0.88, 0.14, -0.72, 
0.84, -1.13, -0.13, 0.2), THERAPY = structure(c(0, 1, 1, 0, 1, 
1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 
1, 1, 1, 0), value.labels = structure(c(1, 0), .Names = c("cognitive", 
"standard"))), ATTRIB = c(-1.17, 0.04, 0.58, -0.23, 0.62, -0.26, 
-0.28, 0.52, 0.34, -0.09, -1.09, 1.05, -1.84, -0.95, 0.15, 0.07, 
-0.1, 2.35, 0.75, 0.49, 0.67, 1.21, 0.31, 1.97, -0.94, 0.11, 
-0.54, -0.23, 0.05, -1.07)), .Names = c("SATIS", "THERAPY", "ATTRIB"
), row.names = c(NA, -30L), class = "data.frame", variable.labels = structure(c("Satisfaction", 
"Therapy", "Attributional Positivity"), .Names = c("SATIS", "THERAPY", 
"ATTRIB")))


###################################################
### code chunk number 36: intro.Rnw:1106-1107
###################################################
preacher <- mediate(SATIS ~ THERAPY + (ATTRIB),data=sobel)  #The example in Preacher and Hayes


###################################################
### code chunk number 37: mediate.png
###################################################
png('mediate.png')
mediate.diagram(preacher)
dev.off()


###################################################
### code chunk number 38: intro.Rnw:1130-1131
###################################################
preacher.lm <- lmCor(SATIS ~ THERAPY + ATTRIB,  data=sobel)  #The example in Preacher and Hayes


###################################################
### code chunk number 39: preacherlm.png
###################################################
png('preacherlm.png')
diagram(preacher.lm)
dev.off()


###################################################
### code chunk number 40: garcia.png
###################################################
png('garcia.png')
model <- mediate(respappr ~ prot2 * sexism +(sexism),data=Garcia,n.iter=50
  ,main="Moderated mediation (mean centered)")
  summary(model)
  dev.off()


###################################################
### code chunk number 41: intro.Rnw:1204-1215
###################################################

dancer  <- structure(list(TS = c(1, 7, 4.6, 1, 7, 7, 7, 7), TC = c(1, 1, 
5.6, 6.6, 4.9, 7, 1, 1), BS = c(1, 7, 7, 1, 7, 6.4, 7, 2.4), 
    BC = c(1, 1, 7, 5.9, 2.9, 3.8, 1, 1)), class = "data.frame", row.names = c(NA, 
-8L))
dancer   #show the data

model <- lmCor(TC + TS ~ BC + BS, data = dancer)
summary(model)  #show the summary statistics
cancorDiagram(model) #and the associated canonical figure



###################################################
### code chunk number 42: dancer
###################################################
png('dancerlm.png')
model <- lmCor(TC + TS ~ BC + BS, data = dancer)
  dev.off()


###################################################
### code chunk number 43: dancer
###################################################
png('dancer.png')
cancorDiagram(model)
  dev.off()


###################################################
### code chunk number 44: intro.Rnw:1264-1268
###################################################

C <- cov(sat.act,use="pairwise")
model1 <- lm(ACT~ gender + education + age, data=sat.act)
summary(model1)


###################################################
### code chunk number 45: intro.Rnw:1271-1273
###################################################
#compare with lmCor
lmCor(c(4:6),c(1:3),C, n.obs=700)


###################################################
### code chunk number 46: intro.Rnw:1394-1395
###################################################
sessionInfo()


