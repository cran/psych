### R code from vignette source 'overview.Rnw'

###################################################
### code chunk number 1: overview.Rnw:334-337
###################################################
library(psych)
data(sat.act) 
describe(sat.act)  #basic descriptive statistics


###################################################
### code chunk number 2: overview.Rnw:344-346
###################################################
 #basic descriptive statistics by a grouping variable.
describeBy(sat.act,sat.act$gender,skew=FALSE,ranges=FALSE)


###################################################
### code chunk number 3: overview.Rnw:354-357
###################################################
sa.mat <- describeBy(sat.act,list(sat.act$gender,sat.act$education),
 skew=FALSE,ranges=FALSE,mat=TRUE)
headTail(sa.mat)


###################################################
### code chunk number 4: outlier
###################################################
png( 'outlier.png' )
d2 <- outlier(sat.act,cex=.8)
dev.off()


###################################################
### code chunk number 5: overview.Rnw:388-392
###################################################
x <- matrix(1:120,ncol=10,byrow=TRUE)
colnames(x) <- paste('V',1:10,sep='')
new.x <- scrub(x,3:5,min=c(30,40,50),max=70,isvalue=45,newvalue=NA)
new.x


###################################################
### code chunk number 6: pairspanels
###################################################
png( 'pairspanels.png' )
sat.d2 <- data.frame(sat.act,d2) #combine the d2 statistics from before with the sat.act data.frame
pairs.panels(sat.d2,bg=c("yellow","blue")[(d2 > 25)+1],pch=21)
dev.off()


###################################################
### code chunk number 7: affect
###################################################
png('affect.png')
pairs.panels(affect[14:17],bg=c("red","black","white","blue")[affect$Film],pch=21,
     main="Affect varies by movies ")
dev.off()


###################################################
### code chunk number 8: overview.Rnw:452-454
###################################################
data(sat.act)
violinBy(sat.act[5:6],sat.act$gender,grp.name=c("M", "F"),main="Density Plot by gender for SAT V and Q")


###################################################
### code chunk number 9: overview.Rnw:479-481
###################################################
data(epi.bfi)
error.bars.by(epi.bfi[,6:10],epi.bfi$epilie<4)


###################################################
### code chunk number 10: overview.Rnw:494-496
###################################################
error.bars.by(sat.act[5:6],sat.act$gender,bars=TRUE,
        labels=c("Male","Female"),ylab="SAT score",xlab="")


###################################################
### code chunk number 11: overview.Rnw:512-522
###################################################
op <- par(mfrow=c(1,2))
  data(affect)
colors <- c("black","red","white","blue")
 films <- c("Sad","Horror","Neutral","Happy")
affect.stats <- errorCircles("EA2","TA2",data=affect[-c(1,20)],group="Film",labels=films,
xlab="Energetic Arousal",   ylab="Tense Arousal",ylim=c(10,22),xlim=c(8,20),pch=16,
cex=2,col=colors, main =' Movies effect on arousal')
 errorCircles("PA2","NA2",data=affect.stats,labels=films,xlab="Positive Affect",
  ylab="Negative Affect", pch=16,cex=2,col=colors,  main ="Movies effect on affect")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 12: overview.Rnw:536-538
###################################################
data(bfi)
with(bfi,{bi.bars(age,gender,ylab="Age",main="Age by males and females")})


###################################################
### code chunk number 13: overview.Rnw:552-553
###################################################
lowerCor(sat.act)


###################################################
### code chunk number 14: overview.Rnw:560-566
###################################################
female <- subset(sat.act,sat.act$gender==2)
 male <- subset(sat.act,sat.act$gender==1)
lower <- lowerCor(male[-1])
upper <- lowerCor(female[-1])
both <- lowerUpper(lower,upper)
round(both,2)


###################################################
### code chunk number 15: overview.Rnw:572-574
###################################################
diffs <-  lowerUpper(lower,upper,diff=TRUE)
round(diffs,2)


###################################################
### code chunk number 16: corplot.png
###################################################
png('corplot.png')
cor.plot(Thurstone,numbers=TRUE,main="9 cognitive variables from Thurstone")
dev.off()


###################################################
### code chunk number 17: circplot.png
###################################################
png('circplot.png')
circ <- sim.circ(24)
r.circ <- cor(circ)
cor.plot(r.circ,main='24 variables in a circumplex')
dev.off()


###################################################
### code chunk number 18: spider.png
###################################################
png('spider.png')
op<- par(mfrow=c(2,2))
spider(y=c(1,6,12,18),x=1:24,data=r.circ,fill=TRUE,main="Spider plot of 24 circumplex variables")
op <- par(mfrow=c(1,1))
dev.off()


###################################################
### code chunk number 19: overview.Rnw:641-642
###################################################
corr.test(sat.act)


###################################################
### code chunk number 20: overview.Rnw:653-654
###################################################
r.test(50,.3)


###################################################
### code chunk number 21: overview.Rnw:660-661
###################################################
r.test(30,.4,.6)


###################################################
### code chunk number 22: overview.Rnw:668-669
###################################################
r.test(103,.4,.5,.1)


###################################################
### code chunk number 23: overview.Rnw:675-676
###################################################
r.test(103,.5,.6,.7,.5,.5,.8)  #steiger Case B 


###################################################
### code chunk number 24: overview.Rnw:684-685
###################################################
cortest(sat.act)


###################################################
### code chunk number 25: overview.Rnw:696-697
###################################################
draw.tetra()


###################################################
### code chunk number 26: overview.Rnw:708-709
###################################################
draw.cor(expand=20,cuts=c(0,0))


###################################################
### code chunk number 27: overview.Rnw:729-730
###################################################
setCor(y = 5:9,x=1:4,data=Thurstone)


###################################################
### code chunk number 28: overview.Rnw:737-739
###################################################
sc <- setCor(y = 5:9,x=3:4,data=Thurstone,z=1:2)
round(sc$residual,2)


###################################################
### code chunk number 29: overview.Rnw:751-765
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
### code chunk number 30: overview.Rnw:767-768
###################################################
preacher <- mediate(1,2,3,sobel)  #The example in Preacher and Hayes


###################################################
### code chunk number 31: overview.Rnw:775-776
###################################################
mediate.diagram(preacher)


###################################################
### code chunk number 32: overview.Rnw:787-789
###################################################
preacher <- setCor(1,c(2,3),sobel,std=FALSE)
setCor.diagram(preacher)


###################################################
### code chunk number 33: overview.Rnw:856-857
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} 


###################################################
### code chunk number 34: overview.Rnw:865-868
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} else {
f3t <- fa(Thurstone,3,n.obs=213)
f3t }


###################################################
### code chunk number 35: overview.Rnw:889-893
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} else {
f3 <- fa(Thurstone,3,n.obs = 213,fm="pa")
f3o <- target.rot(f3)
f3o}


###################################################
### code chunk number 36: overview.Rnw:914-916
###################################################
f3w <- fa(Thurstone,3,n.obs = 213,fm="wls")
print(f3w,cut=0,digits=3)


###################################################
### code chunk number 37: overview.Rnw:928-929
###################################################
plot(f3t)


###################################################
### code chunk number 38: overview.Rnw:941-942
###################################################
fa.diagram(f3t)


###################################################
### code chunk number 39: overview.Rnw:961-963
###################################################
p3p <-principal(Thurstone,3,n.obs = 213,rotate="Promax")
p3p


###################################################
### code chunk number 40: overview.Rnw:982-984
###################################################
om.h <- omega(Thurstone,n.obs=213,sl=FALSE)
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 41: overview.Rnw:995-996
###################################################
om <- omega(Thurstone,n.obs=213)


###################################################
### code chunk number 42: overview.Rnw:1029-1031
###################################################
data(bfi)
ic <- iclust(bfi[1:25])


###################################################
### code chunk number 43: overview.Rnw:1043-1044
###################################################
summary(ic)  #show the results


###################################################
### code chunk number 44: overview.Rnw:1057-1059
###################################################
data(bfi)
r.poly <- polychoric(bfi[1:25]) #the ... indicate the progress of the function


###################################################
### code chunk number 45: overview.Rnw:1072-1074
###################################################
ic.poly <- iclust(r.poly$rho,title="ICLUST using polychoric correlations")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 46: overview.Rnw:1085-1087
###################################################
ic.poly <- iclust(r.poly$rho,5,title="ICLUST using polychoric correlations for nclusters=5")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 47: overview.Rnw:1098-1099
###################################################
ic.poly <- iclust(r.poly$rho,beta.size=3,title="ICLUST beta.size=3")


###################################################
### code chunk number 48: overview.Rnw:1111-1112
###################################################
print(ic,cut=.3)


###################################################
### code chunk number 49: overview.Rnw:1131-1133
###################################################
fa(bfi[1:10],2,n.iter=20)



###################################################
### code chunk number 50: overview.Rnw:1146-1148
###################################################
f4 <- fa(bfi[1:25],4,fm="pa")
factor.congruence(f4,ic)


###################################################
### code chunk number 51: overview.Rnw:1157-1158
###################################################
factor.congruence(list(f3t,f3o,om,p3p))


###################################################
### code chunk number 52: overview.Rnw:1202-1203
###################################################
vss <- vss(bfi[1:25],title="Very Simple Structure of a Big 5 inventory")


###################################################
### code chunk number 53: overview.Rnw:1211-1212
###################################################
vss


###################################################
### code chunk number 54: overview.Rnw:1222-1223
###################################################
fa.parallel(bfi[1:25],main="Parallel Analysis of a Big 5 inventory")


###################################################
### code chunk number 55: overview.Rnw:1241-1246
###################################################
v16 <- sim.item(16)
s <- c(1,3,5,7,9,11,13,15)
f2 <- fa(v16[,s],2)
fe <- fa.extension(cor(v16)[s,-s],f2)
fa.diagram(f2,fe=fe)


###################################################
### code chunk number 56: overview.Rnw:1299-1303
###################################################
set.seed(17)
r9 <- sim.hierarchical(n=500,raw=TRUE)$observed
round(cor(r9),2)
alpha(r9)


###################################################
### code chunk number 57: overview.Rnw:1310-1313
###################################################

alpha(attitude,keys=c("complaints","critical"))



###################################################
### code chunk number 58: overview.Rnw:1320-1322
###################################################

alpha(attitude)


###################################################
### code chunk number 59: overview.Rnw:1329-1331
###################################################
items <- sim.congeneric(N=500,short=FALSE,low=-2,high=2,categorical=TRUE) #500 responses to 4 discrete items
alpha(items$observed)  #item response analysis of congeneric measures


###################################################
### code chunk number 60: overview.Rnw:1384-1385
###################################################
om.9 <- omega(r9,title="9 simulated variables")


###################################################
### code chunk number 61: overview.Rnw:1396-1397
###################################################
om.9


###################################################
### code chunk number 62: overview.Rnw:1405-1406
###################################################
omegaSem(r9,n.obs=500)


###################################################
### code chunk number 63: overview.Rnw:1415-1416
###################################################
splitHalf(r9)


###################################################
### code chunk number 64: overview.Rnw:1438-1456
###################################################
#the old way is by location-- specify the total number of items
keys <- make.keys(nvars=28,list(Agree=c(-1,2:5),Conscientious=c(6:8,-9,-10),
 Extraversion=c(-11,-12,13:15),Neuroticism=c(16:20),
 Openness = c(21,-22,23,24,-25)),
 item.labels=colnames(bfi))
#the newer way is probably preferred -- specify the name of the data set to be scored.
     
keys <- make.keys(bfi,list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C2","-C4","-C5"),
      extraversion=c("-E1","-E2","E3","E4","E5"),neuroticism=c("N1","N2","N3","N4","N5"),
       openness = c("O1","-O2","O3","O4","-O5")) ) #specify the data file to be scored (bfi)
       
#These two approaches can be mixed if desired
 keys.list <- list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C2","-C4","-C5"),
           extraversion=c("-E1","-E2","E3","E4","E5"),
            neuroticism=c(16:20),openness = c(21,-22,23,24,-25)) 
keys <- make.keys(bfi,keys.list) #specify the data file to be scored (bfi)
#In any case, the resulting keys file is just a matrix of -1, 0 and 1s.
  keys


###################################################
### code chunk number 65: overview.Rnw:1463-1467
###################################################
 keys.1<- make.keys(10,list(Agree=c(-1,2:5),Conscientious=c(6:8,-9,-10)))
keys.2 <- make.keys(15,list(Extraversion=c(-1,-2,3:5),Neuroticism=c(6:10),
 Openness = c(11,-12,13,14,-15)))
 keys.25 <- superMatrix(list(keys.1,keys.2))


###################################################
### code chunk number 66: overview.Rnw:1477-1479
###################################################
 scores <- scoreItems(keys,bfi)
 scores


###################################################
### code chunk number 67: scores
###################################################
png('scores.png')
pairs.panels(scores$scores,pch='.',jiggle=TRUE)
dev.off()


###################################################
### code chunk number 68: overview.Rnw:1505-1508
###################################################
r.bfi <- cor(bfi,use="pairwise")
scales <- scoreItems(keys,r.bfi)
summary(scales)


###################################################
### code chunk number 69: overview.Rnw:1518-1524
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
score.multiple.choice(iq.keys,iqitems)
#just convert the items to true or false 
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)
describe(iq.tf)  #compare to previous results


###################################################
### code chunk number 70: overview.Rnw:1542-1548
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
scores <- score.multiple.choice(iq.keys,iqitems,score=TRUE,short=FALSE)
#note that for speed we can just do this on simple item counts rather than IRT based scores.
op <- par(mfrow=c(2,2))  #set this to see the output for multiple items
irt.responses(scores$scores,iqitems[1:4],breaks=11)


###################################################
### code chunk number 71: overview.Rnw:1560-1562
###################################################
 m <- colMeans(bfi,na.rm=TRUE)
  item.lookup(scales$item.corrected[,1:3],m,dictionary=bfi.dictionary[1:2])


###################################################
### code chunk number 72: overview.Rnw:1570-1572
###################################################
data(bfi)
bestScales(bfi,criteria=c("gender","education","age"),cut=.1,dictionary=bfi.dictionary[,1:3])


###################################################
### code chunk number 73: overview.Rnw:1596-1600
###################################################
set.seed(17)
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
test 


###################################################
### code chunk number 74: overview.Rnw:1607-1612
###################################################
op <- par(mfrow=c(3,1))
plot(test,type="ICC")
plot(test,type="IIC")
plot(test,type="test")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 75: overview.Rnw:1623-1626
###################################################
data(bfi)
e.irt <- irt.fa(bfi[11:15])
e.irt 


###################################################
### code chunk number 76: overview.Rnw:1633-1634
###################################################
e.info  <- plot(e.irt,type="IIC")


###################################################
### code chunk number 77: overview.Rnw:1645-1646
###################################################
print(e.info,sort=TRUE)


###################################################
### code chunk number 78: overview.Rnw:1675-1676
###################################################
iq.irt <- irt.fa(iq.tf)


###################################################
### code chunk number 79: overview.Rnw:1688-1689
###################################################
plot(iq.irt,type='test') 


###################################################
### code chunk number 80: overview.Rnw:1700-1701
###################################################
iq.irt 


###################################################
### code chunk number 81: overview.Rnw:1707-1708
###################################################
om <- omega(iq.irt$rho,4)


###################################################
### code chunk number 82: overview.Rnw:1722-1736
###################################################
v9 <- sim.irt(9,1000,-2.,2.,mod="normal") #dichotomous items
items <- v9$items
test <- irt.fa(items)
total <- rowSums(items)
ord <- order(total)
items <- items[ord,]
#now delete some of the data - note that they are ordered by score
items[1:333,5:9] <- NA
items[334:666,3:7] <- NA
items[667:1000,1:4] <- NA
scores <- score.irt(test,items)
unitweighted <- score.irt(items=items,keys=rep(1,9))
scores.df <- data.frame(true=v9$theta[ord],scores,unitweighted)
colnames(scores.df) <- c("True theta","irt theta","total","fit","rasch","total","fit")


###################################################
### code chunk number 83: overview.Rnw:1745-1747
###################################################
 pairs.panels(scores.df,pch='.',gap=0)
 title('Comparing true theta for IRT, Rasch and  classically based scoring',line=3)


###################################################
### code chunk number 84: overview.Rnw:1810-1814
###################################################

C <- cov(sat.act,use="pairwise")
model1 <- lm(ACT~ gender + education + age, data=sat.act)
summary(model1)


###################################################
### code chunk number 85: overview.Rnw:1817-1819
###################################################
#compare with sector
setCor(c(4:6),c(1:3),C, n.obs=700)


###################################################
### code chunk number 86: overview.Rnw:1902-1926
###################################################
xlim=c(0,10)
ylim=c(0,10)
plot(NA,xlim=xlim,ylim=ylim,main="Demontration of dia functions",axes=FALSE,xlab="",ylab="")
ul <- dia.rect(1,9,labels="upper left",xlim=xlim,ylim=ylim)
ll <- dia.rect(1,3,labels="lower left",xlim=xlim,ylim=ylim)
lr <- dia.ellipse(9,3,"lower right",xlim=xlim,ylim=ylim)
ur <- dia.ellipse(9,9,"upper right",xlim=xlim,ylim=ylim)
ml <- dia.ellipse(3,6,"middle left",xlim=xlim,ylim=ylim)
mr <- dia.ellipse(7,6,"middle right",xlim=xlim,ylim=ylim)
bl <- dia.ellipse(1,1,"bottom left",xlim=xlim,ylim=ylim)
br <- dia.rect(9,1,"bottom right",xlim=xlim,ylim=ylim)
dia.arrow(from=lr,to=ul,labels="right to left")
dia.arrow(from=ul,to=ur,labels="left to right")
dia.curved.arrow(from=lr,to=ll$right,labels ="right to left")
dia.curved.arrow(to=ur,from=ul$right,labels ="left to right")
dia.curve(ll$top,ul$bottom,"double",-1)  #for rectangles, specify where to point 
dia.curved.arrow(mr,ur,"up")  #but for ellipses, just point to it.
dia.curve(ml,mr,"across")
dia.arrow(ur,lr,"top down")
dia.curved.arrow(br$top,lr$bottom,"up")
dia.curved.arrow(bl,br,"left to right")
dia.arrow(bl,ll$bottom)
dia.curved.arrow(ml,ll$top,scale=-1)
dia.curved.arrow(mr,lr$top)


###################################################
### code chunk number 87: overview.Rnw:2037-2038
###################################################
sessionInfo()


