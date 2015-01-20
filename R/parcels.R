"sim.parcels" <-
function(x,n.samp=100,reps=10,delta=.5,full=FALSE,congruence=FALSE,max=TRUE) {
result <- matrix(NA,nrow=reps,ncol=8)
n.obs <- nrow(x)
nvar <- ncol(x)

for (i in 1:reps){

samp <- sample(n.obs,n.samp)
r <- cor(x[samp,])
#The normal factor case
f2 <- fa(r,2,rotate="geominQ",delta=delta)
clust <- factor2cluster(f2)
count <- sum(clust[1:(nvar/2),1]) + sum(clust[(nvar/2+1):nvar,2]) 

#the pairwise parcel case
if(full) {
keys <- parcels(r,2,max=max,congruence=congruence)
keys <- keysort(keys)
} else {
keys1 <- parcels(r[1:(nvar/2),1:(nvar/2)],2,max=max,congruence=congruence)
keys2 <- parcels(r[(nvar/2+1):nvar,(nvar/2+1):nvar],2,max=max,congruence=congruence)
keys <- super.matrix(keys1,keys2)
}

x.p <- score.items(keys,x[samp,])
r2 <- cor(x.p$scores)
f2.p <- fa(r2,2,rotate="geominQ",delta=delta)
clust <- factor2cluster(f2.p)
nrclust <- nrow(clust)
count2 <- sum(clust[1:(nrclust/2),1]) + sum(clust[(nrclust/2+1):nrclust,2]) 

Roe <- cor(x.p$scores,x[samp,])
fe2 <- fa.extension(Roe,f2.p)
clustex <- factor2cluster(fe2)
nrclustex <- nrow(clustex)
count2e <- sum(clustex[1:(nrclustex/2),1]) + sum(clustex[(nrclustex/2+1):nrclustex,2]) 



#the tri parcel solution
if(full) {
keys <- parcels(r,3,max=max,congruence=congruence)
keys <- keysort(keys) } else {
keys3.1 <- parcels(r[1:(nvar/2),1:(nvar/2)],3,max=max,congruence=congruence)
keys3.2 <- parcels(r[(nvar/2+1):nvar,(nvar/2+1):nvar],3,max=max,congruence=congruence)
keys <- super.matrix(keys3.1,keys3.2)
}

x.p3 <- scoreItems(keys,x[samp,],)
f2.p3 <- fa(x.p3$scores,2,rotate="geominQ",delta=delta)
clust3 <- factor2cluster(f2.p3)
nrclust3 <- nrow(clust3)
count3 <- sum(clust3[1:(nrclust3/2),1]) + sum(clust3[(nrclust3/2+1):nrclust3,2]) 

Roe <- cor(x.p3$scores,x[samp,])
fe3 <- fa.extension(Roe,f2.p3)
clustex3 <- factor2cluster(fe3)
nrclustex3 <- nrow(clustex3)
count3e <- sum(clustex3[1:(nrclustex3/2),1]) + sum(clustex3[(nrclustex3/2+1):nrclustex3,2]) 


 
result[i,1] <- f2$Phi[1,2]
result[i,2] <- f2.p$Phi[1,2]
result[i,3] <- f2.p3$Phi[1,2]
result[i,4] <- min(count,24-count)
result[i,5] <- min(count2,nrclust-count2)
result[i,6] <- min(count3,nrclust3-count3)
result[i,7] <- min(count2e,nrclustex-count2e)
result[i,8] <- min(count3e,nrclustex3-count3e)
}
colnames(result) <- c("fa Phi","P2 Phi","P3 Phi","error","error2","error3","error2e","error3e")
return(result)
}


"keysort" <-
function(keys) {
items <- 1:nrow(keys)
weights <- items %*% abs(keys)
ord <- order(weights)
keys[] <- keys[,ord]
}


"parcels" <- function
(x,size=3,max=TRUE,flip=TRUE,congruence = FALSE) {

if(nrow(x) !=ncol(x)) {x <- cor(x,use="pairwise")}
if(congruence) { x <- factor.congruence(x,x) }
if(max) {diag(x) <- 0
nvar <- nrow(x)
row.range <- apply(x,1,range,na.rm=TRUE)     
row.max <- pmax(abs(row.range[1,]),abs(row.range[2,]))  #find the largest absolute similarity
diag(x) <- row.max
sd.inv <- diag(1/sqrt(row.max))
similar <- sd.inv %*% x %*% sd.inv
} 
if (size==2) {key <- parcels2(x,flip=flip)} else {key <- parcels3(x,flip=flip)}
rownames(key) <- colnames(x)
colnames(key) <- paste("P",1:ncol(key),sep="")
return(key)}

#this just works for parcels of size 2
"parcels2" <- 
function(x,size=2,flip) {
nvar <- nrow(x)
key <- matrix(0,nvar,nvar/size)
similar <- x
similar <- similar * lower.tri(similar)
for(i in 1:(nvar/size)) {  
max.cell <- which.max(abs(similar))             #global maximum
max.col <- trunc(max.cell/nrow(similar))+1    #is in which row and column?
max.row <- max.cell - (max.col-1)*nrow(similar) #need to fix the case of first column
if (max.row < 1) {max.row <- nrow(similar)
 max.col <- max.col-1 }
 key[max.col,i] <- 1
 if(flip && (similar[max.row,max.col] < 0)) {key[max.row,i] <- -1} else {key[max.row,i] <- 1}
 similar[max.row,] <- similar[,max.row] <- NA
 similar[,max.col] <- similar[max.col,] <- NA
 }
return(key)
}

"parcels3" <-
function(x,flip) {
nvar <- nrow(x)
nkeys <- floor(nvar/3)
keys <- matrix(0,nvar,nkeys)
pointers <- 1:nvar
point <- pointers
for(i in 1:nkeys) {

best <- trimax(abs(x[point,point]))[1:3]
items <- point[best]

keys[items,i] <- 1
if(flip) {if(x[items[1],items[2]] < 0 ) keys[items[2],i] <- -1
          if(x[items[1],items[ 3]] < 0 ) keys[items[3],i] <- -1 
          #if(x[items[2],items[ 3]] < 0 )  keys[items[3],i] <- -keys[items[3],i]   
        } 
point <- point[-best]
 }  
keys
 }
    


"trimax" <- 
function(x) {
  nvar <- nrow(x)
  cs1 <- cumsum(1:nvar)
  cs2 <- cumsum(cs1)
simil <- rep(NA,length=nvar*(nvar-1)* (nvar-2)/6)
ind=1
for(i in 3:nvar) {
  for (j in 2:(i-1)) {
    for (k in 1:(j-1)) {simil[ind] <- x[i,j] + x[i,k] + x[j,k]
    ind <- ind+1}
     }
     } 
 maxi <- which.max(simil)
 if(maxi ==1) {
    m1 <- 3
    m2 <- 2
    m3 <- 1 }  else {
 	m1 <- min(which(cs2 >= maxi)) +2 
	maxi2 <- (maxi - cs2[m1-3]) 
 	if( maxi2 <  2)  {m2 <- 2 
                   m3 <- 1} else {
               		 m2 <- min(which(cs1 >= maxi2))+1
               		 if((maxi2- cs1[m2-2])  > 0 ) { m3 <- maxi2  - cs1[m2-2]
                        } else {m3 <- m2 -1}
    }}
results <- c(m1,m2,m3,maxi)  #maxi is included for debugging
return(results) 
}


#Testing parcels
if(FALSE) {
ind=1
nvar<- 24
cs <- cumsum(cumsum(1:nvar))
print(cs)
for(i in 3:nvar) {
  for (j in 2:(i-1)) {
  
    for (k in 1:(j-1)) {
    m1 <- min(which(cs >= ind)) +2
    newind <- ind - cs[m1-3]
    m2 <- min(which(cs > newind)) +1
   
    cat(ind, i,j,k , "m1 = ",m1,"newind = ",newind, "m2=",m2,", \n")
    
    ind <- ind+1}
     }
     }
}
