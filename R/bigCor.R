#Function to find very big correlation matrices by stitching together smaller ones
#see the help file for various timings
bigCor <- function(x,size=NULL,use="pairwise",cor="pearson"){
  nvar <- NCOL(x)
 if (is.null(size)) size <- 20  #added round 6/21/22  switched to 20  2/22/23
 n.steps <- ceiling(nvar/size)
#first makes sure the data are all numeric  
x <- char2numeric(x)     #added 2/21/23

small.r <- matrix(NA,size,size)
 
short <- function(i) {
res <- list()
 loweri <- (i-1) * size + 1
 upperi <- min(i * size,nvar)
  lowerj <- min((i-1) * size + 1,nvar)
 upperj <- min(i * size,nvar)
#small.r <- cor(x[,loweri:upperi],use=use)
 for (j in i:n.steps) {
 lowerj <- min((j-1) * size + 1,nvar)
 upperj <- min(j * size,nvar)
 if(cor=="pearson") {
 small.r <- cor(x[,loweri:upperi],x[lowerj:upperj],use=use)} else {
  small.r <- t(polychoric(x[loweri:upperi],x[lowerj:upperj])$rho)}  #fixed for polychoric
 res[[j]] <- t(small.r)
 }
return(res)

}

result <- mcmapply(short,c(1:n.steps))

#now, stitch them together
R <- matrix(NA,nvar,nvar)
k <- 1

for(i in 1:NROW(result)) {
 nx <- dim(result[[i,1]])[1]
 
  #temp <- result[[i]]
 temp <- unlist(result[i,])
 
  R[k:(k+ nx-1),1:(k+ nx-1) ] <- temp
  R[1:(k+nx -1),k:(k+ nx-1)] <- t( R[k:(k+ nx-1),1:(k+ nx-1) ])
  k <- k +nx 
  }
colnames(R) <- rownames(R) <- colnames(x)

return(R)
}

#result <- mapply(short,c(1:n.iter),MoreArgs=list(x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=frac,min.item=min.item,max.item=max.item))

#result <- mcmapply(short,c(1:n.iter),MoreArgs=list(x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=frac,min.item=min.item,max.item=max.item))

#options command: options("mc.cores"=4) will set the number of cores to 4.

pairwiseCountBig <- function(x,size=NULL){
  nvar <- NCOL(x)
 if (is.null(size)) size <- nvar/4
 n.steps <- ceiling(nvar/size)

small.r <- matrix(NA,size,size)
 
short <- function(i) {
res <- list()
 loweri <- (i-1) * size + 1
 upperi <- min(i * size,nvar)
  lowerj <- min((i-1) * size + 1,nvar)
 upperj <- min(i * size,nvar)

 for (j in i:n.steps) {
 lowerj <- min((j-1) * size + 1,nvar)
 upperj <- min(j * size,nvar)
 small.r <- pairwiseCount(x[,loweri:upperi],x[lowerj:upperj])
 res[[j]] <- t(small.r)
 }
return(res)

}
result <- mcmapply(short,c(1:n.steps))

#now, stitch them together
R <- matrix(NA,nvar,nvar)
k <- 1
for(i in 1:NROW(result)) {
  nx <- dim(result[[i,1]])[1]
  
  temp <- unlist(result[i,])
  R[k:(k+ nx-1),1:(k+ nx-1) ] <- temp
  R[1:(k+nx -1),k:(k+ nx-1)] <- t( R[k:(k+ nx-1),1:(k+ nx-1) ])
  k <- k +nx 
  }
colnames(R) <- rownames(R) <- colnames(x)

return(R)
}
