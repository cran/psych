"corr.test" <-
function(x,y=NULL,use="pairwise",method="pearson",adjust="holm",alpha=.05,ci=TRUE){
cl <- match.call()
if(is.null(y)) {r <- cor(x,use=use,method=method)
 sym <- TRUE
n <- t(!is.na(x)) %*% (!is.na(x))
} else {r <- cor(x,y,use=use,method=method)
  sym=FALSE
n <- t(!is.na(x)) %*% (!is.na(y))}
if((use=="complete") | (min(n) == max(n))) n <- min(n)
t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
se <- sqrt((1-r*r)/(n-2))

nvar <- ncol(r)

p[p>1] <- 1
if (adjust !="none") {
  if (is.null(y)) {lp <- upper.tri(p)  #the case of a symmetric matrix
     pa <- p[lp]
     pa <- p.adjust(pa,adjust)
     p[upper.tri(p,diag=FALSE)] <- pa
  } else {
  p[] <- p.adjust(p,adjust)  #the case of an asymmetric matrix 
} }
#find confidence intervals
  z <- fisherz(r[lower.tri(r)])
 if(ci) { 
   if (min(n) < 4) {
      warning("Number of subjects must be greater than 3 to find confidence intervals.")
   }
     
   if(adjust!="holm") {dif.corrected <- qnorm(1-alpha/(nvar*(nvar-1))) } else {   # 1- alpha/2  /nvar *(nvar-1) /2)
      ord <- order(abs(z),decreasing=FALSE)  #to find the HOlm correction, we need to order the size of the correlations
      dif.corrected <- qnorm(1-alpha/(2*order(ord))) } #holm

    alpha <- 1-alpha/2  #the raw alpha level for confidence intervals
    dif <- qnorm(alpha)
    if(sym) {
    if(is.matrix(n)) {
  	 sef <- 1/sqrt(n[lower.tri(n)] - 3)
     } else { sef <- 1/sqrt(n - 3)}
    lower <- fisherz2r(z - dif * sef)
    upper <- fisherz2r(z + dif * sef)
   
    lower.corrected <- fisherz2r(z - dif.corrected * sef)
    upper.corrected <- fisherz2r(z + dif.corrected * sef)
     ci <- data.frame(lower=lower,r=r[lower.tri(r)],upper=upper,p=p[lower.tri(p)])
     ci.adj <- data.frame(lower.adj=lower.corrected,upper.adj=upper.corrected)
  

     

      cnR <- abbreviate(colnames(r),minlength=5) 
      
       k <- 1
     for(i in 1:(nvar-1)) {for (j in (i+1):nvar) {
      rownames(ci)[k] <- paste(cnR[i],cnR[j],sep="-")
      k<- k +1 }}
      
    } else { #non symmetric case 
    n.x <- NCOL(x)
    n.y <- NCOL(y)
     z <- fisherz(r)
     if(adjust != "holm") {dif.corrected <- qnorm(1-(1-alpha)/(n.x * n.y)) #we have already adjust alpha by 2
        } else {ord <- order(abs(z),decreasing=FALSE)  #to find the HOlm correction, we need to order the size of the correlations
           dif.corrected <- qnorm(1-(1-alpha)/(order(ord)))
     }
     sef <- 1/sqrt(n - 3)
     lower <- as.vector(fisherz2r(z - dif * sef))
     upper <- as.vector(fisherz2r(z + dif * sef))
     lower.corrected <- fisherz2r(z - dif.corrected * sef)
    upper.corrected <- fisherz2r(z + dif.corrected * sef)
    
    
  ci <- data.frame(lower=lower,r=as.vector(r),upper=upper,p=as.vector(p))
  ci.adj <- data.frame(lower.adj=as.vector(lower.corrected),r=as.vector(r),upper.adj= as.vector(upper.corrected))
  cnR <- abbreviate(rownames(r),minlength=5) 
  cnC <- abbreviate(colnames(r),minlength=5)
  k <- 1
      for(i in 1:ncol(y)) {for (j in 1:ncol(x)) {
      rownames(ci)[k] <- paste(cnR[j],cnC[i],sep="-")
      k<- k +1 }}
    }
} else {ci <- NULL}
result <- list(r = r,n=n,t=t,p=p,se=se,sef=sef, adjust=adjust,sym =sym,ci=ci,ci.adj=ci.adj, Call=cl)
class(result) <- c("psych", "corr.test")
return(result)
}
#modified 1/4/14 to report sample size once if they are all equal
#modified 3/12/14 to report confidence intervals (suggested by Alexander Weiss)
#modified 3/27/14 to correct bug detected by Clemens Fell
#modified 3/27/14 to correct bug reported by Louis-Charles Vannier
#modified 2/21/15 to make confidence intervals an option (incredible decrease in speed if doing cis)
#modified 8/24/17 to include Bonferoni adjusted confidence intervals


"corr.p" <-
function(r,n,adjust="holm",alpha=.05) {
cl <- match.call()
if(missing(n)) stop("The number of subjects must be specified")
sym <- FALSE
t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
p[p>1] <- 1
if (adjust !="none") {
if(isSymmetric(unclass(p))) {sym <- TRUE
 lp <- upper.tri(p)  #the case of a symmetric matrix
     pa <- p[lp]
     pa <- p.adjust(pa,adjust)
     p[upper.tri(p,diag=FALSE)] <- pa
  } else {
  p[] <- p.adjust(p ,adjust)  #the case of an asymmetric matrix
  sym <- FALSE}
} 
nvar <- ncol(r)
if(sym) {z <- fisherz(r[lower.tri(r)])} else {z <- fisherz(r)
 n.x <- NCOL(r)
 n.y <- NROW(r)
 if(adjust != "holm") {dif.corrected <- qnorm((1-alpha/2)/(n.x * n.y)) # adjust alpha by 2
        } else {ord <- order(abs(z),decreasing=FALSE)  #to find the Holm correction, we need to order the size of the correlations
           dif.corrected <- qnorm(1-alpha/(2*(order(ord))))
     }}
 if (min(n) < 4) {
      warning("Number of subjects must be greater than 3 to find confidence intervals.")
   }
   if(sym & is.matrix(n)) {
   se <- 1/sqrt(n[lower.tri(n)] - 3) } else { se <- 1/sqrt(n - 3)}
if(sym) {    dif.corrected <- qnorm(1-alpha/(nvar*(nvar-1)))  } # 1- alpha/2  /nvar *(nvar-1) /2
    alpha <- 1-alpha/2
    dif <- qnorm(alpha)
    lower <- fisherz2r(z - dif * se)
    upper <- fisherz2r(z + dif * se)
      lower.corrected <- fisherz2r(z - dif.corrected * se)
    upper.corrected <- fisherz2r(z + dif.corrected * se)
  if(sym) {ci <- data.frame(lower=lower,r=r[lower.tri(r)],upper=upper,p=p[lower.tri(p)])
          ci.adj <- data.frame(lower.adj = as.vector(lower.corrected),r=r[lower.tri(r)],upper.adj=as.vector(upper.corrected))} else {
   ci <- data.frame(lower=as.vector(lower),r=as.vector(r),upper=as.vector(upper),p=as.vector(p))
      ci.adj <- data.frame(lower.adj =as.vector( lower.corrected),r=as.vector(r),upper.adj= as.vector(upper.corrected))}
      cnR <- abbreviate(colnames(r),minlength=5)  
      rnR <- abbreviate(rownames(r),minlength=5) 
     if(sym) {k <- 1
      for(i in 2:nvar) {for (j in 1:(i-1)) {
      rownames(ci)[k] <- paste(cnR[j],cnR[i],sep="-")
      k<- k +1 } }
    
      
      
      } else {k <- 1
      for(i in 1:ncol(r)) {for (j in 1:nrow(r)) {
      rownames(ci)[k] <- paste(rnR[j],cnR[i],sep="-")
      k<- k +1 }}
      }
result <- list(r = r,n=n,t=t,p=p,sym=sym,adjust=adjust,ci=ci,ci.adj = ci.adj,Call=cl)
class(result) <- c("psych", "corr.p")
return(result)
}
#revised March 28, 2014 to be compatible with corr.test
#revised August 28, 2017 to include holm and bonferroini adjusted confidence intervals

#could be replaced with the following
corr.test1 <- function(x,y=NULL,use="pairwise",method="pearson",adjust="holm",alpha=.05){
cl <- match.call()
if(is.null(y)) {r <- cor(x,use=use,method=method)
 sym <- TRUE
n <- t(!is.na(x)) %*% (!is.na(x))
} else {r <- cor(x,y,use=use,method=method)
  sym=FALSE
n <- t(!is.na(x)) %*% (!is.na(y))}
if((use=="complete") | (min(n) == max(n))) n <- min(n)

result <- corr.p(r,n,adjust=adjust,alpha=alpha)
result$Call<- cl
class(result) <- c("psych", "corr.test")
return(result)
 }
 
 
 
 
 


