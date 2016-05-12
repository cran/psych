"phi2poly" <-
function(ph,cp,cc,n=NULL,correct=TRUE) {
	#.Deprecated(phi2tetra, msg='phi2poly is deprecated, please use phi2tetra')
     #ph is the phi coefficient
     #cp is the selection ratio of the predictor
     #cc is the success rate of the criterion
     r.marg<-rep(0,2)
     c.marg<- rep(0,2)
     p<-array(rep(0,4),dim=c(2,2))
     r.marg[1]<- cp
     r.marg[2]<- 1 -cp 
     c.marg[1]<- cc
     c.marg[2]<- 1-cc
     
     p[1,1]<- r.marg[1]*c.marg[1]+ ph*sqrt(prod(r.marg,c.marg))
     p[2,2]<- r.marg[2]*c.marg[2]+ ph*sqrt(prod(r.marg,c.marg))
     p[1,2]<- r.marg[1]*c.marg[2]- ph*sqrt(prod(r.marg,c.marg))
     p[2,1]<- r.marg[2]*c.marg[1]- ph*sqrt(prod(r.marg,c.marg))
     if(!is.null(n)) p <- p*n
     result<-tetrachoric(p,correct=correct )$rho 
     return(result)}

"phi2tet" <-
function(ph,cp,cc,n=NULL,correct=TRUE) {
	if(is.null(n)) n <- 1 
     #ph is the phi coefficient
     #cp is the selection ratio of the predictor
     #cc is the success rate of the criterion
     r.marg<-rep(0,2)
     c.marg<- rep(0,2)
     p<-array(rep(0,4),dim=c(2,2))
     r.marg[1]<- cp/n
     r.marg[2]<- 1 -cp/n 
     c.marg[1]<- cc/n
     c.marg[2]<- 1-cc/n
     
     p[1,1]<- r.marg[1]*c.marg[1]+ ph*sqrt(prod(r.marg,c.marg))
     p[2,2]<- r.marg[2]*c.marg[2]+ ph*sqrt(prod(r.marg,c.marg))
     p[1,2]<- r.marg[1]*c.marg[2]- ph*sqrt(prod(r.marg,c.marg))
     p[2,1]<- r.marg[2]*c.marg[1]- ph*sqrt(prod(r.marg,c.marg))
     if(!is.null(n)) p <- p*n
     result<-tetrachoric(p,correct=correct )$rho 
     return(result)}

"phi2tetra" <-
function(ph,m,n=NULL,correct=TRUE) {
if(!is.matrix(ph) && !is.data.frame(ph)) {result <- phi2tet(ph,m[1],m[2],n=n,correct=correct) } else {
nvar <- nrow(ph)
if(nvar !=ncol(ph)) {stop('Matrix must be square')}
if (length(m) !=nvar) {stop("length of m must match the number of variables")}
result <- as.matrix(ph)
for(i in 2:nvar) {
  for (j in 1:(i-1)) { 
   result[i,j] <- result[j,i] <- phi2tet(ph[i,j],m[i],m[j],n=n,correct=correct) 
  }
}
}
return(result) }
	
	
   