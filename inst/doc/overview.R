### R code from vignette source 'overview.Rnw'

###################################################
### code chunk number 1: overview.Rnw:325-328
###################################################
library(psych)
data(sat.act) 
describe(sat.act)  #basic descriptive statistics


###################################################
### code chunk number 2: overview.Rnw:335-337
###################################################
 #basic descriptive statistics by a grouping variable.
describeBy(sat.act,sat.act$gender,skew=FALSE,ranges=FALSE)


###################################################
### code chunk number 3: overview.Rnw:345-348
###################################################
sa.mat <- describeBy(sat.act,list(sat.act$gender,sat.act$education),
 skew=FALSE,ranges=FALSE,mat=TRUE)
headTail(sa.mat)


###################################################
### code chunk number 4: overview.Rnw:359-363
###################################################
x <- matrix(1:120,ncol=10,byrow=TRUE)
colnames(x) <- paste('V',1:10,sep='')
new.x <- scrub(x,3:5,min=c(30,40,50),max=70,isvalue=45,newvalue=NA)
new.x


###################################################
### code chunk number 5: pairspanels
###################################################
png( 'pairspanels.png' )
pairs.panels(sat.act,pch='.')
dev.off()


###################################################
### code chunk number 6: affect
###################################################
png('affect.png')
pairs.panels(affect[14:17],bg=c("red","black","white","blue")[affect$Film],pch=21,
     main="Affect varies by movies ")
dev.off()


###################################################
### code chunk number 7: overview.Rnw:430-432
###################################################
data(epi.bfi)
error.bars.by(epi.bfi[,6:10],epi.bfi$epilie<4)


###################################################
### code chunk number 8: overview.Rnw:445-447
###################################################
error.bars.by(sat.act[5:6],sat.act$gender,bars=TRUE,
        labels=c("Male","Female"),ylab="SAT score",xlab="")


###################################################
### code chunk number 9: overview.Rnw:463-472
###################################################
op <- par(mfrow=c(1,2))
  data(affect)
colors <- c("black","red","white","blue")
 films <- c("Sad","Horror","Neutral","Happy")
affect.stats <- errorCircles("EA2","TA2",data=affect,group="Film",labels=films,xlab="Energetic Arousal",ylab="Tense Arousal",ylim=c(10,22),xlim=c(8,20),pch=16,cex=2,col=colors,
  main =' Movies effect on arousal')
 errorCircles("PA2","NA2",data=affect.stats,labels=films,xlab="Positive Affect",ylab="Negative Affect",pch=16,cex=2,col=colors,
   main ="Movies effect on affect")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 10: overview.Rnw:486-488
###################################################
data(bfi)
with(bfi,{bi.bars(age,gender,ylab="Age",main="Age by males and females")})


###################################################
### code chunk number 11: overview.Rnw:502-503
###################################################
lowerCor(sat.act)


###################################################
### code chunk number 12: overview.Rnw:510-516
###################################################
female <- subset(sat.act,sat.act$gender==2)
 male <- subset(sat.act,sat.act$gender==1)
lower <- lowerCor(male[-1])
upper <- lowerCor(female[-1])
both <- lowerUpper(lower,upper)
round(both,2)


###################################################
### code chunk number 13: overview.Rnw:522-524
###################################################
diffs <-  lowerUpper(lower,upper,diff=TRUE)
round(diffs,2)


###################################################
### code chunk number 14: corplot.png
###################################################
png('corplot.png')
cor.plot(Thurstone,numbers=TRUE,main="9 cognitive variables from Thurstone")
dev.off()


###################################################
### code chunk number 15: circplot.png
###################################################
png('circplot.png')
circ <- sim.circ(24)
r.circ <- cor(circ)
cor.plot(r.circ,main='24 variables in a circumplex')
dev.off()


###################################################
### code chunk number 16: spider.png
###################################################
png('spider.png')
op<- par(mfrow=c(2,2))
spider(y=c(1,6,12,18),x=1:24,data=r.circ,fill=TRUE,main="Spider plot of 24 circumplex variables")
op <- par(mfrow=c(1,1))
dev.off()


###################################################
### code chunk number 17: overview.Rnw:591-592
###################################################
corr.test(sat.act)


###################################################
### code chunk number 18: overview.Rnw:603-604
###################################################
r.test(50,.3)


###################################################
### code chunk number 19: overview.Rnw:610-611
###################################################
r.test(30,.4,.6)


###################################################
### code chunk number 20: overview.Rnw:618-619
###################################################
r.test(103,.4,.5,.1)


###################################################
### code chunk number 21: overview.Rnw:625-626
###################################################
r.test(103,.5,.6,.7,.5,.5,.8)  #steiger Case B 


###################################################
### code chunk number 22: overview.Rnw:634-635
###################################################
cortest(sat.act)


###################################################
### code chunk number 23: overview.Rnw:646-647
###################################################
draw.tetra()


###################################################
### code chunk number 24: overview.Rnw:666-667
###################################################
set.cor(y = 5:9,x=1:4,data=Thurstone)


###################################################
### code chunk number 25: overview.Rnw:674-676
###################################################
sc <- set.cor(y = 5:9,x=3:4,data=Thurstone,z=1:2)
round(sc$residual,2)


###################################################
### code chunk number 26: overview.Rnw:742-744
###################################################
f3t <- fa(Thurstone,3,n.obs=213)
f3t 


###################################################
### code chunk number 27: overview.Rnw:765-768
###################################################
f3 <- fa(Thurstone,3,n.obs = 213,fm="pa")
f3o <- target.rot(f3)
f3o


###################################################
### code chunk number 28: overview.Rnw:789-791
###################################################
f3w <- fa(Thurstone,3,n.obs = 213,fm="wls")
print(f3w,cut=0,digits=3)


###################################################
### code chunk number 29: overview.Rnw:803-804
###################################################
plot(f3t)


###################################################
### code chunk number 30: overview.Rnw:816-817
###################################################
fa.diagram(f3t)


###################################################
### code chunk number 31: overview.Rnw:836-838
###################################################
p3p <-principal(Thurstone,3,n.obs = 213,rotate="Promax")
p3p


###################################################
### code chunk number 32: overview.Rnw:857-859
###################################################
om.h <- omega(Thurstone,n.obs=213,sl=FALSE)
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 33: overview.Rnw:870-871
###################################################
om <- omega(Thurstone,n.obs=213)


###################################################
### code chunk number 34: overview.Rnw:904-906
###################################################
data(bfi)
ic <- iclust(bfi[1:25])


###################################################
### code chunk number 35: overview.Rnw:918-919
###################################################
summary(ic)  #show the results


###################################################
### code chunk number 36: overview.Rnw:932-934
###################################################
data(bfi)
r.poly <- polychoric(bfi[1:25]) #the ... indicate the progress of the function


###################################################
### code chunk number 37: overview.Rnw:947-949
###################################################
ic.poly <- iclust(r.poly$rho,title="ICLUST using polychoric correlations")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 38: overview.Rnw:960-962
###################################################
ic.poly <- iclust(r.poly$rho,5,title="ICLUST using polychoric correlations for nclusters=5")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 39: overview.Rnw:973-974
###################################################
ic.poly <- iclust(r.poly$rho,beta.size=3,title="ICLUST beta.size=3")


###################################################
### code chunk number 40: overview.Rnw:986-987
###################################################
print(ic,cut=.3)


###################################################
### code chunk number 41: overview.Rnw:1006-1008
###################################################
fa(bfi[1:10],2,n.iter=20)



###################################################
### code chunk number 42: overview.Rnw:1021-1023
###################################################
f4 <- fa(bfi[1:25],4,fm="pa")
factor.congruence(f4,ic)


###################################################
### code chunk number 43: overview.Rnw:1032-1033
###################################################
factor.congruence(list(f3t,f3o,om,p3p))


###################################################
### code chunk number 44: overview.Rnw:1077-1078
###################################################
vss <- vss(bfi[1:25],title="Very Simple Structure of a Big 5 inventory")


###################################################
### code chunk number 45: overview.Rnw:1086-1087
###################################################
vss


###################################################
### code chunk number 46: overview.Rnw:1097-1098
###################################################
fa.parallel(bfi[1:25],main="Parallel Analysis of a Big 5 inventory")


###################################################
### code chunk number 47: overview.Rnw:1116-1121
###################################################
v16 <- sim.item(16)
s <- c(1,3,5,7,9,11,13,15)
f2 <- fa(v16[,s],2)
fe <- fa.extension(cor(v16)[s,-s],f2)
fa.diagram(f2,fe=fe)


###################################################
### code chunk number 48: overview.Rnw:1174-1178
###################################################
set.seed(17)
r9 <- sim.hierarchical(n=500,raw=TRUE)$observed
round(cor(r9),2)
alpha(r9)


###################################################
### code chunk number 49: overview.Rnw:1185-1188
###################################################

alpha(attitude,keys=c("complaints","critical"))



###################################################
### code chunk number 50: overview.Rnw:1195-1197
###################################################

alpha(attitude)


###################################################
### code chunk number 51: overview.Rnw:1204-1206
###################################################
items <- sim.congeneric(N=500,short=FALSE,low=-2,high=2,categorical=TRUE) #500 responses to 4 discrete items
alpha(items$observed)  #item response analysis of congeneric measures


###################################################
### code chunk number 52: overview.Rnw:1259-1260
###################################################
om.9 <- omega(r9,title="9 simulated variables")


###################################################
### code chunk number 53: overview.Rnw:1271-1272
###################################################
om.9


###################################################
### code chunk number 54: overview.Rnw:1280-1281
###################################################
omegaSem(r9,n.obs=500)


###################################################
### code chunk number 55: overview.Rnw:1290-1291
###################################################
guttman(r9)


###################################################
### code chunk number 56: overview.Rnw:1313-1331
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
### code chunk number 57: overview.Rnw:1338-1342
###################################################
 keys.1<- make.keys(10,list(Agree=c(-1,2:5),Conscientious=c(6:8,-9,-10)))
keys.2 <- make.keys(15,list(Extraversion=c(-1,-2,3:5),Neuroticism=c(6:10),
 Openness = c(11,-12,13,14,-15)))
 keys.25 <- superMatrix(list(keys.1,keys.2))


###################################################
### code chunk number 58: overview.Rnw:1352-1354
###################################################
 scores <- score.items(keys,bfi)
 scores


###################################################
### code chunk number 59: scores
###################################################
png('scores.png')
pairs.panels(scores$scores,pch='.',jiggle=TRUE)
dev.off()


###################################################
### code chunk number 60: overview.Rnw:1380-1383
###################################################
r.bfi <- cor(bfi,use="pairwise")
scales <- cluster.cor(keys,r.bfi)
summary(scales)


###################################################
### code chunk number 61: overview.Rnw:1393-1399
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
score.multiple.choice(iq.keys,iqitems)
#just convert the items to true or false 
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)
describe(iq.tf)  #compare to previous results


###################################################
### code chunk number 62: overview.Rnw:1417-1423
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
scores <- score.multiple.choice(iq.keys,iqitems,score=TRUE,short=FALSE)
#note that for speed we can just do this on simple item counts rather than IRT based scores.
op <- par(mfrow=c(2,2))  #set this to see the output for multiple items
irt.responses(scores$scores,iqitems[1:4],breaks=11)


###################################################
### code chunk number 63: overview.Rnw:1449-1453
###################################################
set.seed(17)
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
test 


###################################################
### code chunk number 64: overview.Rnw:1460-1465
###################################################
op <- par(mfrow=c(3,1))
plot(test,type="ICC")
plot(test,type="IIC")
plot(test,type="test")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 65: overview.Rnw:1476-1479
###################################################
data(bfi)
e.irt <- irt.fa(bfi[11:15])
e.irt 


###################################################
### code chunk number 66: overview.Rnw:1486-1487
###################################################
e.info  <- plot(e.irt,type="IIC")


###################################################
### code chunk number 67: overview.Rnw:1498-1499
###################################################
print(e.info,sort=TRUE)


###################################################
### code chunk number 68: overview.Rnw:1513-1514
###################################################
iq.irt <- irt.fa(iq.tf)


###################################################
### code chunk number 69: overview.Rnw:1526-1527
###################################################
plot(iq.irt,type='test') 


###################################################
### code chunk number 70: overview.Rnw:1538-1539
###################################################
iq.irt 


###################################################
### code chunk number 71: overview.Rnw:1545-1546
###################################################
om <- omega(iq.irt$rho,4)


###################################################
### code chunk number 72: overview.Rnw:1560-1574
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
### code chunk number 73: overview.Rnw:1583-1585
###################################################
 pairs.panels(scores.df,pch='.',gap=0)
 title('Comparing true theta for IRT, Rasch and  classically based scoring',line=3)


###################################################
### code chunk number 74: overview.Rnw:1634-1638
###################################################

C <- cov(sat.act,use="pairwise")
model1 <- lm(ACT~ gender + education + age, data=sat.act)
summary(model1)


###################################################
### code chunk number 75: overview.Rnw:1641-1643
###################################################
#compare with mat.regress
set.cor(c(4:6),c(1:3),C, n.obs=700)


###################################################
### code chunk number 76: overview.Rnw:1726-1750
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
### code chunk number 77: overview.Rnw:1861-1862
###################################################
sessionInfo()


