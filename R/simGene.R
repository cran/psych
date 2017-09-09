"simGene" <- function(ng=10,traits=1,n.obs=1000,dom=TRUE) {
 X <- array(sample(2,ng*n.obs*traits*3,replace=TRUE),dim=c(n.obs,ng,traits,3))
 MZ <- DZ <- array(NA,dim=c(n.obs,ng,traits))
  MZt <- DZt <- matrix(NA,n.obs,traits)
for(t in 1:traits) {
if(dom) { MZ[,1:ng,t] <- X[,1:ng,t,1] * X[,1:ng,t,2] #the allele values are mulitplied
 DZ[,1:ng,t] <- X[,1:ng,t,1] * X[,1:ng,t,3]} else {
 MZ[,1:ng,t] <- X[,1:ng,t,1] + X[,1:ng,t,2]    #the allele values are added
 DZ[,1:ng,t] <- X[,1:ng,t,1] + X[,1:ng,t,3]}

 
 MZt[,t] <- rowMeans(MZ[,,t])      #the trait values
 DZt[,t] <- rowMeans(DZ[,,t])
 }
  X.df <- data.frame(genes=X[,1:ng,1:traits,sample(2,1,replace=TRUE)],MZ=MZt,DZ=DZt)
 return(X.df)}
 
test.simGene <- function(x, ng=10) {
 t1 <-rowMeans(x[1:ng])
 t2 <- rowMeans(x[(ng+1):(ng*2)])
 t11 <-rowMeans(x[1:(ng/2)])
 t12 <-rowMeans(x[(ng/2 +1):ng])
 t21 <-rowMeans(x[(ng+1):(ng/2 + ng)])
 t22 <- rowMeans(x[(ng/2 + ng+1):(ng*2)])

 scores <- data.frame(t1=t1,t2=t2,t11=t11,t12=t12,t21 = t21,t22=t22,traits=x[(ng*2 +1):(ng*2+4)])
 }