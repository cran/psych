### R code from vignette source 'overview.Rnw'

###################################################
### code chunk number 1: overview.Rnw:448-449
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} 


###################################################
### code chunk number 2: overview.Rnw:457-462
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} else {
library(psych)
library(psychTools)
f3t <- fa(Thurstone,3,n.obs=213)
f3t }


###################################################
### code chunk number 3: overview.Rnw:483-487
###################################################
if(!require('GPArotation')) {stop('GPArotation must be installed to do rotations')} else {
f3 <- fa(Thurstone,3,n.obs = 213,fm="pa")
f3o <- target.rot(f3)
f3o}


###################################################
### code chunk number 4: overview.Rnw:510-512
###################################################
f3w <- fa(Thurstone,3,n.obs = 213,fm="wls")
print(f3w,cut=0,digits=3)


###################################################
### code chunk number 5: overview.Rnw:525-526
###################################################
plot(f3t)


###################################################
### code chunk number 6: overview.Rnw:538-539
###################################################
fa.diagram(f3t)


###################################################
### code chunk number 7: overview.Rnw:558-560
###################################################
p3p <-principal(Thurstone,3,n.obs = 213,rotate="Promax")
p3p


###################################################
### code chunk number 8: overview.Rnw:579-581
###################################################
om.h <- omega(Thurstone,n.obs=213,sl=FALSE)
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 9: overview.Rnw:592-593
###################################################
om <- omega(Thurstone,n.obs=213)


###################################################
### code chunk number 10: overview.Rnw:626-628
###################################################
data(bfi)
ic <- iclust(bfi[1:25])


###################################################
### code chunk number 11: overview.Rnw:640-641
###################################################
summary(ic)  #show the results


###################################################
### code chunk number 12: overview.Rnw:654-656
###################################################
data(bfi)
r.poly <- polychoric(bfi[1:25],correct=0) #the ... indicate the progress of the function


###################################################
### code chunk number 13: overview.Rnw:668-670
###################################################
ic.poly <- iclust(r.poly$rho,title="ICLUST using polychoric correlations")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 14: overview.Rnw:681-683
###################################################
ic.poly <- iclust(r.poly$rho,5,title="ICLUST using polychoric correlations for nclusters=5")
iclust.diagram(ic.poly) 


###################################################
### code chunk number 15: overview.Rnw:694-695
###################################################
ic.poly <- iclust(r.poly$rho,beta.size=3,title="ICLUST beta.size=3")


###################################################
### code chunk number 16: overview.Rnw:707-708
###################################################
print(ic,cut=.3)


###################################################
### code chunk number 17: overview.Rnw:731-733
###################################################
fa(bfi[1:10],2,n.iter=20)



###################################################
### code chunk number 18: overview.Rnw:746-748
###################################################
f4 <- fa(bfi[1:25],4,fm="pa")
factor.congruence(f4,ic)


###################################################
### code chunk number 19: overview.Rnw:757-758
###################################################
factor.congruence(list(f3t,f3o,om,p3p))


###################################################
### code chunk number 20: overview.Rnw:802-803
###################################################
vss <- vss(bfi[1:25],title="Very Simple Structure of a Big 5 inventory")


###################################################
### code chunk number 21: overview.Rnw:811-812
###################################################
vss


###################################################
### code chunk number 22: overview.Rnw:822-823
###################################################
fa.parallel(bfi[1:25],main="Parallel Analysis of a Big 5 inventory")


###################################################
### code chunk number 23: overview.Rnw:843-848
###################################################
v16 <- sim.item(16)
s <- c(1,3,5,7,9,11,13,15)
f2 <- fa(v16[,s],2)
fe <- fa.extension(cor(v16)[s,-s],f2)
fa.diagram(f2,fe=fe)


###################################################
### code chunk number 24: overview.Rnw:864-871
###################################################
fx <-matrix(c( .9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)  
fy <- matrix(c(.6,.5,.4),ncol=1)
rownames(fx) <- c("V","Q","A","nach","Anx")
rownames(fy)<- c("gpa","Pre","MA")
Phi <-matrix( c(1,0,.7,.0,1,.7,.7,.7,1),ncol=3)
gre.gpa <- sim.structural(fx,Phi,fy)
print(gre.gpa)


###################################################
### code chunk number 25: overview.Rnw:877-879
###################################################
esem.example <- esem(gre.gpa$model,varsX=1:5,varsY=6:8,nfX=2,nfY=1,n.obs=1000,plot=FALSE)
esem.example


###################################################
### code chunk number 26: overview.Rnw:884-885
###################################################
 esem.diagram(esem.example)


###################################################
### code chunk number 27: overview.Rnw:938-942
###################################################
set.seed(17)
r9 <- sim.hierarchical(n=500,raw=TRUE)$observed
round(cor(r9),2)
alpha(r9)


###################################################
### code chunk number 28: overview.Rnw:949-952
###################################################

alpha(attitude,keys=c("complaints","critical"))



###################################################
### code chunk number 29: overview.Rnw:959-961
###################################################

alpha(attitude)


###################################################
### code chunk number 30: overview.Rnw:968-970
###################################################
items <- sim.congeneric(N=500,short=FALSE,low=-2,high=2,categorical=TRUE) #500 responses to 4 discrete items
alpha(items$observed)  #item response analysis of congeneric measures


###################################################
### code chunk number 31: overview.Rnw:1023-1024
###################################################
om.9 <- omega(r9,title="9 simulated variables")


###################################################
### code chunk number 32: overview.Rnw:1035-1036
###################################################
om.9


###################################################
### code chunk number 33: overview.Rnw:1044-1045
###################################################
omegaSem(r9,n.obs=500,lavaan=TRUE)


###################################################
### code chunk number 34: overview.Rnw:1054-1055
###################################################
splitHalf(r9)


###################################################
### code chunk number 35: overview.Rnw:1069-1089
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
 keys.list <- list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C3","-C4","-C5"),
           extraversion=c("-E1","-E2","E3","E4","E5"),
            neuroticism=c(16:20),openness = c(21,-22,23,24,-25)) 

  keys.list


###################################################
### code chunk number 36: overview.Rnw:1111-1113
###################################################
 scores <- scoreItems(keys.list,bfi)
 scores


###################################################
### code chunk number 37: scores
###################################################
png('scores.png')
pairs.panels(scores$scores,pch='.',jiggle=TRUE)
dev.off()


###################################################
### code chunk number 38: overview.Rnw:1139-1142
###################################################
r.bfi <- cor(bfi,use="pairwise")
scales <- scoreItems(keys.list,r.bfi)
summary(scales)


###################################################
### code chunk number 39: overview.Rnw:1152-1158
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
score.multiple.choice(iq.keys,iqitems)
#just convert the items to true or false 
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)
describe(iq.tf)  #compare to previous results


###################################################
### code chunk number 40: overview.Rnw:1176-1182
###################################################
data(iqitems)
iq.keys <- c(4,4,4, 6,6,3,4,4,  5,2,2,4,  3,2,6,7)
scores <- score.multiple.choice(iq.keys,iqitems,score=TRUE,short=FALSE)
#note that for speed we can just do this on simple item counts rather than IRT based scores.
op <- par(mfrow=c(2,2))  #set this to see the output for multiple items
irt.responses(scores$scores,iqitems[1:4],breaks=11)


###################################################
### code chunk number 41: overview.Rnw:1194-1196
###################################################
 m <- colMeans(bfi,na.rm=TRUE)
  item.lookup(scales$item.corrected[,1:3],m,dictionary=bfi.dictionary[1:2])


###################################################
### code chunk number 42: overview.Rnw:1204-1206
###################################################
data(bfi)
bestScales(bfi,criteria=c("gender","education","age"),cut=.1,dictionary=bfi.dictionary[,1:3])


###################################################
### code chunk number 43: overview.Rnw:1230-1234
###################################################
set.seed(17)
d9 <- sim.irt(9,1000,-2.0,2.0,mod="normal") #dichotomous items
test <- irt.fa(d9$items,correct=0)
test 


###################################################
### code chunk number 44: overview.Rnw:1241-1246
###################################################
op <- par(mfrow=c(3,1))
plot(test,type="ICC")
plot(test,type="IIC")
plot(test,type="test")
op <- par(mfrow=c(1,1))


###################################################
### code chunk number 45: overview.Rnw:1257-1260
###################################################
data(bfi)
e.irt <- irt.fa(bfi[11:15])
e.irt 


###################################################
### code chunk number 46: overview.Rnw:1267-1268
###################################################
e.info  <- plot(e.irt,type="IIC")


###################################################
### code chunk number 47: overview.Rnw:1279-1280
###################################################
print(e.info,sort=TRUE)


###################################################
### code chunk number 48: overview.Rnw:1309-1310
###################################################
iq.irt <- irt.fa(ability)


###################################################
### code chunk number 49: overview.Rnw:1322-1323
###################################################
plot(iq.irt,type='test') 


###################################################
### code chunk number 50: overview.Rnw:1334-1335
###################################################
iq.irt 


###################################################
### code chunk number 51: overview.Rnw:1341-1342
###################################################
om <- omega(iq.irt$rho,4)


###################################################
### code chunk number 52: overview.Rnw:1356-1370
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
### code chunk number 53: overview.Rnw:1379-1381
###################################################
 pairs.panels(scores.df,pch='.',gap=0)
 title('Comparing true theta for IRT, Rasch and  classically based scoring',line=3)


###################################################
### code chunk number 54: overview.Rnw:1393-1409
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
### code chunk number 55: overview.Rnw:1414-1418
###################################################
#compare the solutions using the cor2 function
 cor2(bfi.1pl,bfi.ctt)
cor2(bfi.2pl,bfi.ctt)
 cor2(bfi.2pl,bfi.1pl)


###################################################
### code chunk number 56: overview.Rnw:1482-1486
###################################################

C <- cov(sat.act,use="pairwise")
model1 <- lm(ACT~ gender + education + age, data=sat.act)
summary(model1)


###################################################
### code chunk number 57: overview.Rnw:1489-1491
###################################################
#compare with setCor
setCor(gender + education + age ~ ACT + SATV + SATQ, data = C, n.obs=700)


###################################################
### code chunk number 58: overview.Rnw:1574-1598
###################################################
xlim=c(0,10)
ylim=c(0,10)
plot(NA,xlim=xlim,ylim=ylim,main="Demonstration of dia functions",axes=FALSE,xlab="",ylab="")
ul <- dia.rect(1,9,labels="upper left",xlim=xlim,ylim=ylim)
ll <- dia.rect(1,3,labels="lower left",xlim=xlim,ylim=ylim)
lr <- dia.ellipse(9,3,"lower right",xlim=xlim,ylim=ylim,e.size=.09)
ur <- dia.ellipse(7,9,"upper right",xlim=xlim,ylim=ylim,e.size=.1)
ml <- dia.ellipse(3,6,"middle left",xlim=xlim,ylim=ylim,e.size=.1)
mr <- dia.ellipse(7,6,"middle right",xlim=xlim,ylim=ylim,e.size=.08)
bl <- dia.ellipse(1,1,"bottom left",xlim=xlim,ylim=ylim,e.size=.08)
br <- dia.rect(9,1,"bottom right",xlim=xlim,ylim=ylim)
dia.arrow(from=lr,to=ul,labels="right to left")
dia.arrow(from=ul,to=ur,labels="left to right")
dia.curved.arrow(from=lr,to=ll$right,labels ="right to left")
dia.curved.arrow(to=ur,from=ul$right,labels ="left to right")
dia.curve(ll$top,ul$bottom,"double",-1)  #for rectangles, specify where to point 
dia.curved.arrow(mr,ur,"up")  #but for ellipses, just point to it.
dia.curve(ml,mr,"across")
dia.curved.arrow(ur,lr,"top down")
dia.curved.arrow(br$top,lr$bottom,"up")
dia.curved.arrow(bl,br,"left to right")
dia.arrow(bl$top,ll$bottom)
dia.curved.arrow(ml,ll$top,scale=-1)
dia.curved.arrow(mr,lr$top)


###################################################
### code chunk number 59: overview.Rnw:1719-1720
###################################################
sessionInfo()


