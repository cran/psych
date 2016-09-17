### R code from vignette source 'overview.Rnw'

###################################################
### code chunk number 1: overview.Rnw:338-341
###################################################
library(psych)
data(sat.act) 
describe(sat.act)  #basic descriptive statistics


###################################################
### code chunk number 2: overview.Rnw:348-350
###################################################
 #basic descriptive statistics by a grouping variable.
describeBy(sat.act,sat.act$gender,skew=FALSE,ranges=FALSE)


###################################################
### code chunk number 3: overview.Rnw:358-361
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
### code chunk number 5: overview.Rnw:392-396
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
### code chunk number 8: overview.Rnw:456-458
###################################################
data(sat.act)
violinBy(sat.act[5:6],sat.act$gender,grp.name=c("M", "F"),main="Density Plot by gender for SAT V and Q")


###################################################
### code chunk number 9: overview.Rnw:484-486
###################################################
data(epi.bfi)
error.bars.by(epi.bfi[,6:10],epi.bfi$epilie<4)


###################################################
### code chunk number 10: overview.Rnw:499-501
###################################################
error.bars.by(sat.act[5:6],sat.act$gender,bars=TRUE,
        labels=c("Male","Female"),ylab="SAT score",xlab="")


###################################################
### code chunk number 11: overview.Rnw:515-519
###################################################
T <- with(sat.act,table(gender,education))
rownames(T) <- c("M","F")
error.bars.tab(T,way="both",ylab="Proportion of Education Level",xlab="Level of Education",
main="Proportion of sample by education level")


###################################################
### code chunk number 12: overview.Rnw:538-548
###################################################
op <- par(mfrow=c(1,2))
  data(affect)
colors <- c("black","red","white","blue")
 films <- c("Sad","Horror","Neutral","Happy")
affect.stats <- errorCircles("EA2","TA2",data=affect[-c(1,20)],group="Film",labels=films,
xlab="Energetic Arousal",   ylab="Tense Arousal",ylim=c(10,22),xlim=c(8,20),pch=16,
cex=2,colorsl=colors, main =' Movies effect on arousal')
 errorCircles("PA2","NA2",data=affect.stats,labels=films,xlab="Positive Affect",
  ylab="Negative Affect", pch=16,cex=2,colors=colors,  main ="Movies effect on affect")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 13: overview.Rnw:562-564
###################################################
data(bfi)
with(bfi,{bi.bars(age,gender,ylab="Age",main="Age by males and females")})


###################################################
### code chunk number 14: overview.Rnw:580-581
###################################################
lowerCor(sat.act)


###################################################
### code chunk number 15: overview.Rnw:588-594
###################################################
female <- subset(sat.act,sat.act$gender==2)
 male <- subset(sat.act,sat.act$gender==1)
lower <- lowerCor(male[-1])
upper <- lowerCor(female[-1])
both <- lowerUpper(lower,upper)
round(both,2)


###################################################
### code chunk number 16: overview.Rnw:600-602
###################################################
diffs <-  lowerUpper(lower,upper,diff=TRUE)
round(diffs,2)


###################################################
### code chunk number 17: corplot.png
###################################################
png('corplot.png')
corPlot(Thurstone,numbers=TRUE,upper=FALSE,diag=FALSE,main="9 cognitive variables from Thurstone")
dev.off()


###################################################
### code chunk number 18: circplot.png
###################################################
png('circplot.png')
circ <- sim.circ(24)
r.circ <- cor(circ)
corPlot(r.circ,main='24 variables in a circumplex')
dev.off()


###################################################
### code chunk number 19: spider.png
###################################################
png('spider.png')
op<- par(mfrow=c(2,2))
spider(y=c(1,6,12,18),x=1:24,data=r.circ,fill=TRUE,main="Spider plot of 24 circumplex variables")
op <- par(mfrow=c(1,1))
dev.off()


###################################################
### code chunk number 20: overview.Rnw:669-670
###################################################
corr.test(sat.act)


###################################################
### code chunk number 21: overview.Rnw:681-682
###################################################
r.test(50,.3)


###################################################
### code chunk number 22: overview.Rnw:688-689
###################################################
r.test(30,.4,.6)


###################################################
### code chunk number 23: overview.Rnw:696-697
###################################################
r.test(103,.4,.5,.1)


###################################################
### code chunk number 24: overview.Rnw:703-704
###################################################
r.test(103,.5,.6,.7,.5,.5,.8)  #steiger Case B 


###################################################
### code chunk number 25: overview.Rnw:712-713
###################################################
cortest(sat.act)


###################################################
### code chunk number 26: overview.Rnw:724-725
###################################################
draw.tetra()


###################################################
### code chunk number 27: overview.Rnw:736-737
###################################################
draw.cor(expand=20,cuts=c(0,0))


###################################################
### code chunk number 28: overview.Rnw:757-758
###################################################
setCor(y = 5:9,x=1:4,data=Thurstone)


###################################################
### code chunk number 29: overview.Rnw:765-767
###################################################
sc <- setCor(y = 5:9,x=3:4,data=Thurstone,z=1:2)
round(sc$residual,2)


###################################################
### code chunk number 30: overview.Rnw:779-793
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
### code chunk number 31: overview.Rnw:795-796
###################################################
preacher <- mediate(1,2,3,sobel)  #The example in Preacher and Hayes


###################################################
### code chunk number 32: overview.Rnw:803-804
###################################################
mediate.diagram(preacher)


###################################################
### code chunk number 33: overview.Rnw:815-817
###################################################
preacher <- setCor(1,c(2,3),sobel,std=FALSE)
setCor.diagram(preacher)


###################################################
### code chunk number 34: overview.Rnw:884-885
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} 


###################################################
### code chunk number 35: overview.Rnw:893-896
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} else {
f3t <- fa(Thurstone,3,n.obs=213)
f3t }


###################################################
### code chunk number 36: overview.Rnw:917-921
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} else {
f3 <- fa(Thurstone,3,n.obs = 213,fm="pa")
f3o <- target.rot(f3)
f3o}


###################################################
### code chunk number 37: overview.Rnw:942-944
###################################################
f3w <- fa(Thurstone,3,n.obs = 213,fm="wls")
print(f3w,cut=0,digits=3)


###################################################
### code chunk number 38: overview.Rnw:956-957
###################################################
plot(f3t)


###################################################
### code chunk number 39: overview.Rnw:969-970
###################################################
fa.diagram(f3t)


###################################################
### code chunk number 40: overview.Rnw:989-991
###################################################
p3p <-principal(Thurstone,3,n.obs = 213,rotate="Promax")
p3p


###################################################
### code chunk number 41: overview.Rnw:1010-1012
###################################################
om.h <- omega(Thurstone,n.obs=213,sl=FALSE)
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 42: overview.Rnw:1023-1024
###################################################
om <- omega(Thurstone,n.obs=213)


###################################################
### code chunk number 43: overview.Rnw:1057-1059
###################################################
data(bfi)
ic <- iclust(bfi[1:25])


###################################################
### code chunk number 44: overview.Rnw:1071-1072
###################################################
summary(ic)  #show the results


###################################################
### code chunk number 45: overview.Rnw:1085-1087
###################################################
data(bfi)
r.poly <- polychoric(bfi[1:25]) #the ... indicate the progress of the function


###################################################
### code chunk number 46: overview.Rnw:1100-1102
###################################################
ic.poly <- iclust(r.poly$rho,title="ICLUST using polychoric correlations")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 47: overview.Rnw:1113-1115
###################################################
ic.poly <- iclust(r.poly$rho,5,title="ICLUST using polychoric correlations for nclusters=5")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 48: overview.Rnw:1126-1127
###################################################
ic.poly <- iclust(r.poly$rho,beta.size=3,title="ICLUST beta.size=3")


###################################################
### code chunk number 49: overview.Rnw:1139-1140
###################################################
print(ic,cut=.3)


###################################################
### code chunk number 50: overview.Rnw:1159-1161
###################################################
fa(bfi[1:10],2,n.iter=20)



###################################################
### code chunk number 51: overview.Rnw:1174-1176
###################################################
f4 <- fa(bfi[1:25],4,fm="pa")
factor.congruence(f4,ic)


###################################################
### code chunk number 52: overview.Rnw:1185-1186
###################################################
factor.congruence(list(f3t,f3o,om,p3p))


###################################################
### code chunk number 53: overview.Rnw:1230-1231
###################################################
vss <- vss(bfi[1:25],title="Very Simple Structure of a Big 5 inventory")


###################################################
### code chunk number 54: overview.Rnw:1239-1240
###################################################
vss


###################################################
### code chunk number 55: overview.Rnw:1250-1251
###################################################
fa.parallel(bfi[1:25],main="Parallel Analysis of a Big 5 inventory")


###################################################
### code chunk number 56: overview.Rnw:1269-1274
###################################################
v16 <- sim.item(16)
s <- c(1,3,5,7,9,11,13,15)
f2 <- fa(v16[,s],2)
fe <- fa.extension(cor(v16)[s,-s],f2)
fa.diagram(f2,fe=fe)


###################################################
### code chunk number 57: overview.Rnw:1290-1297
###################################################
fx <-matrix(c( .9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)  
fy <- matrix(c(.6,.5,.4),ncol=1)
rownames(fx) <- c("V","Q","A","nach","Anx")
rownames(fy)<- c("gpa","Pre","MA")
Phi <-matrix( c(1,0,.7,.0,1,.7,.7,.7,1),ncol=3)
gre.gpa <- sim.structural(fx,Phi,fy)
print(gre.gpa)


###################################################
### code chunk number 58: overview.Rnw:1303-1305
###################################################
esem.example <- esem(gre.gpa$model,varsX=1:5,varsY=6:8,nfX=2,nfY=1,n.obs=1000,plot=FALSE)
esem.example


###################################################
### code chunk number 59: overview.Rnw:1310-1311
###################################################
 esem.diagram(esem.example)


###################################################
### code chunk number 60: overview.Rnw:1363-1367
###################################################
set.seed(17)
r9 <- sim.hierarchical(n=500,raw=TRUE)$observed
round(cor(r9),2)
alpha(r9)


###################################################
### code chunk number 61: overview.Rnw:1374-1377
###################################################

alpha(attitude,keys=c("complaints","critical"))



###################################################
### code chunk number 62: overview.Rnw:1384-1386
###################################################

alpha(attitude)


###################################################
### code chunk number 63: overview.Rnw:1393-1395
###################################################
items <- sim.congeneric(N=500,short=FALSE,low=-2,high=2,categorical=TRUE) #500 responses to 4 discrete items
alpha(items$observed)  #item response analysis of congeneric measures


###################################################
### code chunk number 64: overview.Rnw:1448-1449
###################################################
om.9 <- omega(r9,title="9 simulated variables")


###################################################
### code chunk number 65: overview.Rnw:1460-1461
###################################################
om.9


###################################################
### code chunk number 66: overview.Rnw:1469-1470
###################################################
omegaSem(r9,n.obs=500)


###################################################
### code chunk number 67: overview.Rnw:1479-1480
###################################################
splitHalf(r9)


###################################################
### code chunk number 68: overview.Rnw:1494-1514
###################################################
 #the newer way is probably preferred 
     
keys.list  <- list(agree=c("-A1","A2","A3","A4","A5"),
      conscientious=c("C1","C2","C2","-C4","-C5"),
      extraversion=c("-E1","-E2","E3","E4","E5"),
      neuroticism=c("N1","N2","N3","N4","N5"),
      openness = c("O1","-O2","O3","O4","-O5")) 
      
      
#this can also be done  by location--
keys.list  <- list(Agree=c(-1,2:5),Conscientious=c(6:8,-9,-10),
 Extraversion=c(-11,-12,13:15),Neuroticism=c(16:20),
 Openness = c(21,-22,23,24,-25))
          
#These two approaches can be mixed if desired
 keys.list <- list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C2","-C4","-C5"),
           extraversion=c("-E1","-E2","E3","E4","E5"),
            neuroticism=c(16:20),openness = c(21,-22,23,24,-25)) 

  keys.list


###################################################
### code chunk number 69: overview.Rnw:1552-1554
###################################################
 scores <- scoreItems(keys.list,bfi)
 scores


###################################################
### code chunk number 70: scores
###################################################
png('scores.png')
pairs.panels(scores$scores,pch='.',jiggle=TRUE)
dev.off()


###################################################
### code chunk number 71: overview.Rnw:1580-1583
###################################################
r.bfi <- cor(bfi,use="pairwise")
scales <- scoreItems(keys.list,r.bfi)
summary(scales)


###################################################
### code chunk number 72: overview.Rnw:1593-1599
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
score.multiple.choice(iq.keys,iqitems)
#just convert the items to true or false 
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)
describe(iq.tf)  #compare to previous results


###################################################
### code chunk number 73: overview.Rnw:1617-1623
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
scores <- score.multiple.choice(iq.keys,iqitems,score=TRUE,short=FALSE)
#note that for speed we can just do this on simple item counts rather than IRT based scores.
op <- par(mfrow=c(2,2))  #set this to see the output for multiple items
irt.responses(scores$scores,iqitems[1:4],breaks=11)


###################################################
### code chunk number 74: overview.Rnw:1635-1637
###################################################
 m <- colMeans(bfi,na.rm=TRUE)
  item.lookup(scales$item.corrected[,1:3],m,dictionary=bfi.dictionary[1:2])


###################################################
### code chunk number 75: overview.Rnw:1645-1647
###################################################
data(bfi)
bestScales(bfi,criteria=c("gender","education","age"),cut=.1,dictionary=bfi.dictionary[,1:3])


###################################################
### code chunk number 76: overview.Rnw:1671-1675
###################################################
set.seed(17)
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
test 


###################################################
### code chunk number 77: overview.Rnw:1682-1687
###################################################
op <- par(mfrow=c(3,1))
plot(test,type="ICC")
plot(test,type="IIC")
plot(test,type="test")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 78: overview.Rnw:1698-1701
###################################################
data(bfi)
e.irt <- irt.fa(bfi[11:15])
e.irt 


###################################################
### code chunk number 79: overview.Rnw:1708-1709
###################################################
e.info  <- plot(e.irt,type="IIC")


###################################################
### code chunk number 80: overview.Rnw:1720-1721
###################################################
print(e.info,sort=TRUE)


###################################################
### code chunk number 81: overview.Rnw:1750-1751
###################################################
iq.irt <- irt.fa(iq.tf)


###################################################
### code chunk number 82: overview.Rnw:1763-1764
###################################################
plot(iq.irt,type='test') 


###################################################
### code chunk number 83: overview.Rnw:1775-1776
###################################################
iq.irt 


###################################################
### code chunk number 84: overview.Rnw:1782-1783
###################################################
om <- omega(iq.irt$rho,4)


###################################################
### code chunk number 85: overview.Rnw:1797-1811
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
scores <- scoreIrt(test,items)
unitweighted <- scoreIrt(items=items,keys=rep(1,9))
scores.df <- data.frame(true=v9$theta[ord],scores,unitweighted)
colnames(scores.df) <- c("True theta","irt theta","total","fit","rasch","total","fit")


###################################################
### code chunk number 86: overview.Rnw:1820-1822
###################################################
 pairs.panels(scores.df,pch='.',gap=0)
 title('Comparing true theta for IRT, Rasch and  classically based scoring',line=3)


###################################################
### code chunk number 87: overview.Rnw:1834-1850
###################################################
keys.list  <- list(agree=c("-A1","A2","A3","A4","A5"),
      conscientious=c("C1","C2","C3","-C4","-C5"),
      extraversion=c("-E1","-E2","E3","E4","E5"),
      neuroticism=c("N1","N2","N3","N4","N5"),
      openness = c("O1","-O2","O3","O4","-O5")) 
    
   item.list <- list(agree=c("A1","A2","A3","A4","A5"),
      conscientious=c("C1","C2","C3","C4","C5"),
      extraversion=c("E1","E2","E3","E4","E5"),
      neuroticism=c("N1","N2","N3","N4","N5"),
      openness = c("O1","O2","O3","O4","O5")) 

bfi.1pl <- scoreIrt.1pl(keys.list,bfi)  #the one parameter solution
bfi.2pl <- scoreIrt.2pl(item.list,bfi)  #the two parameter solution
bfi.ctt <- scoreFast(keys.list,bfi) # fast scoring function



###################################################
### code chunk number 88: overview.Rnw:1855-1859
###################################################
#compare the solutions using the cor2 function
 cor2(bfi.1pl,bfi.ctt)
cor2(bfi.2pl,bfi.ctt)
 cor2(bfi.2pl,bfi.1pl)


###################################################
### code chunk number 89: overview.Rnw:1916-1920
###################################################

C <- cov(sat.act,use="pairwise")
model1 <- lm(ACT~ gender + education + age, data=sat.act)
summary(model1)


###################################################
### code chunk number 90: overview.Rnw:1923-1925
###################################################
#compare with sector
setCor(c(4:6),c(1:3),C, n.obs=700)


###################################################
### code chunk number 91: overview.Rnw:2008-2032
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
### code chunk number 92: overview.Rnw:2144-2145
###################################################
sessionInfo()


