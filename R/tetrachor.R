#adapted from John Fox's Polychor
#this does all the work
"tetra" <- 
function(x,y=NULL,correct=TRUE) {
 binBvn <- function (rho,rc,cc)    #adapted from John Fox's polychor
{ row.cuts <- c(-Inf, rc, Inf)
    col.cuts <- c(-Inf, cc, Inf)
    P <- matrix(0, 2,2)
    R <- matrix(c(1, rho, rho, 1), 2, 2)
    for (i in 1:2) {
        for (j in 1:2) {
            P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), 
                upper = c(row.cuts[i + 1], col.cuts[j + 1]), 
                corr = R)
        }}
    P   #the estimated 2 x 2 predicted by rho, rc, cc
}
 f <- function(rho,cc,rc) { 
      P <- binBvn(rho, cc, rc) 
       -sum(tab * log(P)) }  #the ML criterion to be minimized
      
 if(is.null(y)) {tab <- x} else {tab <- table(x,y)
 if((sum(tab) > 1) && (min(tab) == 0) && correct) {
    warning("A cell entry of 0 was replaced with .5.  Check your data!")
    tab[tab==0] <-.5  #correction for continuity
    }
    }
  tot <- sum(tab)
  tab <- tab/tot
  rc <- qnorm(colSums(tab))[1]
  cc <- qnorm(rowSums(tab))[1]
  rho <- optimize(f,interval=c(-1,1),rc=rc,cc=cc)
  result <- list(rho=rho$minimum,tau=c(cc,rc),objective=rho$objective)
  return(result)
  }
  
 #repeatedly do the analysis to form a matrix of output 
"tetra.mat" <- 
function(x,correct=TRUE) {nvar <- dim(x)[2]
x <- x -min(x,na.rm=TRUE) #in case the numbers are not 0,1
if(max(x,na.rm=TRUE) > 1) {stop("Tetrachoric correlations require dictomous data")}
tau <- -qnorm(colMeans(x,na.rm=TRUE))
mat <- matrix(0,nvar,nvar)
colnames(mat) <- rownames(mat) <- colnames(x)
names(tau) <- colnames(x)
for (i in 2:nvar) {
  for (j in 1:(i-1)) {
  tetra <-  tetra(x[,i],x[,j],correct=correct)
   mat[i,j] <- mat[j,i] <- tetra$rho
   }
   }
   diag(mat) <- 1
  result <- list(rho = mat,tau = tau)
  return(result) 
  }
  
 #convert comorbidity type numbers to a table
 pqr <- function(q1,q2=NULL,p=NULL) {
    if(length(q1) > 1) {
       q2 <- q1[2]
       p <- q1[3]
       q1 <- q1[1]}
   tab <- matrix(0,2,2)
   tab[1,1] <- p
   tab[2,1] <- q1-p
   tab[1,2] <- q2-p
   tab[2,2] <- 1-q1 - tab[1,2]
   return(tab)}
   
 #the public function
 "tetrachor" <- 
 function(x,correct=TRUE) {
 cl <- match.call() 
 if (!is.matrix(x) && !is.data.frame(x)) {
  if (length(x) ==4) {x <- matrix(x,2,2) } else {
  if(length(x) ==3 ) {x <- pqr(x) } else {
  stop("Data must be either a 1 x 4 vector, a 2 x 2 matrix, or a data.frame/matrix of data")} 
  }}
  nvar <- dim(x)[2]
 if (dim(x)[1] == nvar) {result <- tetra(x,correct=correct)} else {
 result <- tetra.mat(x,correct=correct)}
 
 result$Call <- cl
 class(result) <- c("psych","tetra")
 return(result) 
 }