"cohen.d" <- function(x,group,alpha=.05,std=TRUE,sort=NULL,dictionary=NULL,MD=TRUE) {
cl <- match.call()
if ((length(group) ==1) && ( group %in% colnames(x) )) {group <- which(colnames(x) %in% group)
  group.in <- TRUE} else {group.in <- FALSE}
 stats <- statsBy(x,group)
 S <- stats$rwg

 S.inv <- Pinv(S)   #the pseudo inverse because it is possible this is not PSD   added December 15, 2019
d <- stats$mean[2,] - stats$mean[1,]
sd.p <- sqrt((( (stats$n[1,]-1) * stats$sd[1,]^2) + (stats$n[2,]-1) * stats$sd[2,]^2)/(stats$n[1,]+stats$n[2,])) #if we subtract 2 from n, we get Hedges g
sd.ph <- sqrt((((stats$n[1,]-1) * stats$sd[1,]^2) + (stats$n[2,]-1) * stats$sd[2,]^2)/(stats$n[1,]+stats$n[2,]-2)) #if we subtract 2 from n, we get Hedges g
n <- stats$n[1,]+ stats$n[2,]
cohen.d <- d/sd.p

hedges.g <- d/sd.ph
d <- cohen.d    #basically use this in the Mahalanobis distance
names(cohen.d) <- colnames(x)
if(!group.in) {group <- which(colnames(stats$n)=="group")}
n <- n[-group]
d <- d[-group]
n1 <- stats$n[1,-group]
n2 <- stats$n[2,-group]
p1 <- n1/n
p2 <- n2/n
cohen.d <- cohen.d[-group]
hedges.g <- hedges.g[-group]

r <- cohen.d/sqrt(cohen.d^2 + 1/(  p1*p2)) #} else {r <- cohen.d/sqrt(cohen.d^2 + 1/( p1*p2) )} #for unequal n otherwise this is just 4
t <- d2t(cohen.d,n)
p <- ( 1-pt(abs(t),n-2)) * 2
if(MD) {D <- sqrt(t(d) %*% S.inv %*% d)   #convert to D units from D2 units
wt.d <- t(d) %*% S.inv *d    #what is this?
D <- as.vector(D)} else {D <- NA}
cohen.d.conf <- cohen.d.ci(cohen.d,n1=n1,n2=n2,alpha=alpha)
if(!is.null(dictionary)) {dict <- dictionary[match(rownames(cohen.d.conf),rownames(dictionary)),,drop=FALSE]
   cohen.d.conf <- cbind(cohen.d.conf,dict)} else {dict=NULL}
if(!is.null(sort)) {if(sort %in%( c("decreasing","descending","TRUE"))) {ord <- order(cohen.d.conf["effect"],decreasing=TRUE)} else {ord <- order(cohen.d.conf["effect"],decreasing=FALSE)}
       cohen.d.conf <- cohen.d.conf[ord,]
       dict <- dict[ord,,drop=FALSE]
      }
se <- (cohen.d.conf[,3] - cohen.d.conf[,1])/2  #average upper - lower 
result <- list(cohen.d = cohen.d.conf,hedges.g = hedges.g, M.dist = D, r=r,t=t,n=n,p=p, wt.d =wt.d,descriptive=stats,se=se,dict=dict,Call=cl)
class(result) <- c("psych","cohen.d")
return(result)
}

"d2t" <- function(d,n=NULL,n2=NULL,n1=NULL) {if(is.null(n1)) {t <- d*sqrt(n)/2} else 
    if(is.null(n2)) {t <- d*sqrt(n1) } else {
    t <- d /sqrt(1/n1 + 1/n2)}
       return(t)}
   
"t2d" <- function(t,n=NULL,n2=NULL, n1=NULL) {if(is.null(n1)) { d <- 2*t/sqrt(n)} else {
   if(is.null(n2)) { d <- t/sqrt(n1)} else {
     d <- t * sqrt(1/n1 + 1/n2)}}
  return(d)}
   
"d.ci" <- "cohen.d.ci" <- function(d,n=NULL,n2=NULL,n1=NULL,alpha=.05) { t <- d2t(d=d,n=n,n2=n2,n1=n1)
    tail <- 1- alpha/2 
    ci <- matrix(NA,ncol=3,nrow=length(d))
    for(i in 1:length(d)) {
    nmax <- pmax(c(n,n1+1,n1+n2))

 upper <- try(t2d( uniroot(function(x) {suppressWarnings(pt(q=t[i],df=nmax[i]-2,ncp=x)) - alpha/2}, c(min(-5,-abs(t[i])*10),max(5,abs(t[i])*10)))$root,n=n[i],n2=n2[i],n1=n1[i]),silent=TRUE)
  
     if(inherits( upper, "try-error")) {ci[i,3] <- NA} else {ci[i,3] <- upper}
    ci[i,2] <- d[i]
  lower.ci  <- try(t2d(uniroot(function(x) {suppressWarnings(pt(q=t[i],df=nmax[i]-2,ncp=x)) - tail}, c(min(-5,-abs(t[i])*10),max(5,abs(t[i]) *10)))$root,n=n[i],n2=n2[i],n1=n1[i]),silent=TRUE)
    if(inherits( lower.ci,"try-error")) {ci[i,1] <- NA} else {ci[i,1] <- lower.ci}
   }
   colnames(ci) <- c("lower","effect","upper")
   rownames(ci) <- names(d)
    return(ci)
     }
     
"m2t" <- function(m1,m2,s1,s2,n1=NULL,n2=NULL,n=NULL,pooled=TRUE ) { 
     if(!is.null(n) ) { 
        t <- (m1-m2)/sqrt((s1^2 + s2^2)/(n/2))
        d <- 2*t/sqrt(n)
        df <- n-2} else {
    if(pooled) {vp <- ((n1-1) * s1^2 +  (n2-1)* s2^2)/(n1+n2 -2 )
            se = sqrt(vp*(1/n1 + 1/n2))} else {se = sqrt(s1^2/n1 + s2^2/n2)}
      t <- (m1-m2)/se
      df=n1 +n2 -2
      if(!pooled) {df = (s1^2/n1 + s2^2/n2)^2/(s1^4/(n1^2 *(n1-1)) + s2^4/(n2^2 * (n2-1)))}
        d <- t * sqrt(1/n1 + 1/n2)}
      p <- 2* pt(abs(t),df,lower.tail=FALSE)
      result <- list(t=t,df=df,p= p,d=d)
     cat("\n t = ",t, "df =", df, " with probability = 
",p,"\n")
invisible(result) #return the values as well
     }
     

  
 
 "cohen.d.by" <- 
function(x,group,group2,alpha=.05,MD=TRUE)  {


  group1 <- group
  group1name <- group
  group2name <- group2
   group2 <- which(colnames(x) %in% group2)
   group1 <- which(colnames(x) %in% group)
    categories <- names(table(x[group2]))
    result <- list()
    for(i in 1:length(categories)) {
      group <- subset(x,x[group2]==categories[i])
      group <- group[-group2]
      result[[i]] <- cohen.d(group,group1name,MD=MD)     
      }
      names(result) <- paste0(group1name,"for",group2name ,categories)
      class(result) <- class(result) <- c("psych","cohen.d.by")
      return(result)
   }
    
    
    "print.cohen.d" <- function(x,digits=2) {cat("Call: ")
           print(x$Call)
            cat("Cohen d statistic of difference between two means\n")
            if(NCOL(x$cohen.d) ==  3)  {print(round(x$cohen.d,digits=digits))} else {print( data.frame(round(x$cohen.d[1:3],digits=digits),x$cohen.d[4:NCOL(x$cohen.d)]))}
          
            cat("\nMultivariate (Mahalanobis) distance between groups\n")
            print(x$M.dist,digits=digits) 
            cat("r equivalent of difference between two means\n")
            print(round(x$r,digits=digits))
            }
 
    "print.cohen.d.by" <- function(x,digits=2) {cat("Call: ")
            print(x$Call)
            ncases <- length(x)
            for (i in (1:ncases)) {cat("\n Group levels = ",names(x[i]),"\n")
               cat("Cohen d statistic of difference between two means\n")
            print(x[[i]]$cohen.d,digits=digits)
            cat("\nMultivariate (Mahalanobis) distance between groups\n")
            print(x[[i]]$M.dist,digits=digits) 
            cat("r equivalent dof difference between two means\n")
            print(x[[i]]$r,digits=digits)
             }
            }
            

#Following Algina 2015
"d.robust" <- function(x,group,trim=.2) {

 valid <- function(x) {
        sum(!is.na(x))
    }
  nvar <- NCOL(x)

  means <- list()
  vars <- list()
  Sw <- d.robust <-  rep(NA,nvar)
  n.by.grp <- list()
  if(nvar ==1) {
     means[1] <- by(x,group,function(x) mean(x,trim = trim, na.rm=TRUE))
     vars[1] <- by(x,group,function(x) winsor.var(x,trim=trim,na.rm=TRUE))
     } else {  cn <- colnames(x)
      for (i in 1:nvar) {
      n.by.grp[[cn[i]]] <-  by(x[,i],x[group],valid)
       means[[cn[i]]] <- by(x[,i],x[group],function(x) mean(x,trim = trim, na.rm=TRUE))
       vars[[cn[i]]] <- by(x[,i],x[group],function(x) winsor.var(x,trim = trim, na.rm=TRUE))
     }
 }
mean.by.grp <- matrix(unlist(means),ncol=2,byrow=TRUE)
vars.by.grp <- matrix(unlist(vars),ncol=2,byrow=TRUE)
n.by.grp  <- matrix(unlist(n.by.grp),ncol=2,byrow=TRUE)
rownames(mean.by.grp) <- cn
rownames(vars.by.grp) <- cn
colnames(mean.by.grp) <-colnames(vars.by.grp) <- paste0("Grp",1:2)
for(i in 1:nvar) {
Sw[i] = sqrt((vars.by.grp[i,1] * (n.by.grp[i,1]-1) + vars.by.grp[i,2] * (n.by.grp[i,2]-1))/(n.by.grp[i,1] + n.by.grp[i,2]-2))
d.robust[i] <- .642 * (mean.by.grp[i,2] - mean.by.grp[i,1])/Sw[i]
names(d.robust) <- cn
}

 result <- list(means=mean.by.grp,vars=vars.by.grp,Sw,d.robust) 
 return(result)   
}

#find the resampled M.dist November 3, 2018
cohen.d.expected <- function(x,group,n.rep=10 ) {
  summary <- list()
  n.obs <- nrow(x)
  observed <- cohen.d(x=x,group=group)$M.dist
 ind <- 1:n.obs
  for(i in 1:n.rep){
  samp <- sample(ind,n.obs,replace=FALSE)  #this is a random permutation of the order variable
  x[,group] <- x[samp,group]
  summary[[i]] <- cohen.d(x,group)$M.dist

  }
  result <- unlist(summary)
  mean.boot <- mean(result)
  sd.boot <- sd(result)
  result <-list(observed=observed,mean = mean.boot,sd=sd.boot,trials =result)
  return(result)
  }    
 

  