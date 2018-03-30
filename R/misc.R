#A number of useful helper functions
#added January, 2012

#a dummy function to allow the help to find misc
"psych.misc" <-
function() {}

"lowerMat" <- 
function(R,digits=2) {
   lowleft <- lower.tri(R,diag=TRUE)
   nvar <- ncol(R)
	nc <- digits+3
	width <- getOption("width")
	k1 <- width/(nc+2)
   if(is.null(colnames(R))) {colnames(R) <- paste("C",1:nvar,sep="")}
   if(is.null(rownames(R))) {rownames(R) <- paste("R",1:nvar,sep="")}
	colnames(R) <- abbreviate(colnames(R),minlength=digits+3)
	
	nvar <- ncol(R)
	nc <- digits+3
	#k1 <- width/(nc+2)
	
	if(k1 * nvar < width) {k1 <- nvar}  
	k1 <- floor(k1)
	fx <- format(round(R,digits=digits))
	 if(nrow(R) == ncol(R) ) {fx[!lowleft] <- ""}
	 for(k in seq(0,nvar,k1)) { if(k<nvar) {
	print(fx[(k+1):nvar,(k+1):min((k1+k),nvar)],quote=FALSE)}
	}
	invisible(R[lower.tri(R,diag=FALSE)])}
	
"lowerCor" <- 
function(x,digits=2,use="pairwise",method="pearson") {
   R <- cor(x,use=use,method=method)
   lowerMat(R,digits)
   invisible(R)
   }

	
#adapted from utils::txtProgressBar
#modified August 10, 2012 to print just 100 times. 
"progressBar" <- 
function(value,max,label=NULL) {
if(class(stdout())[1]=="terminal") { pc <- round(100 * value/max)  #only print to the screen, not to a file
if(ceiling(100 * value/max)==floor(100 * value/max)) {
width <- 100
char="."
nw <- nchar(char, "w")
nb <- round(width * value/max )
       
        #cat(paste(c("\r  ",label," |", rep.int(" ", nw * width + 6)), collapse = ""))
        cat(paste(c("\r  ",label," |", rep.int(char, nb), rep.int(" ", nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""))
        }
            }
 #flush.console()
 flush(stdout())
 }


"reflect" <- 
function(f,flip=NULL) {
rnames <- colnames(f$loadings)
rnames[flip] <- paste(rnames[flip],'(R)',sep="")
flipper <- rep(1,ncol(f$loadings))
flipper[flip] <- -1
flip <- diag(flipper)
colnames(flip) <- rownames(flip) <- rnames
f$loadings <- f$loadings %*% flip
if(!is.null(f$weights)) f$weights <- f$weights %*% flip
if(!is.null(f$scores)) f$scores <- f$scores %*% flip
if(!is.null(f$Phi)) f$Phi <- flip %*% f$Phi %*% t(flip)
return(f)
}
#and futher patched May 10, 2017 to get the right denominator (I was dividing by n-1 - lag instead of n-lag
#patched March 1, 2017 to get the df right for lags 
#developed January 10, 2012 to find the mean square of successive differences
#see Von Neuman et al. 1941
#added the minus 1 to the case of a single vector.  Sept 21, 2017
#the matrix version was correct, just not the single case version
"mssd" <- function(x,group=NULL,lag=1,na.rm=TRUE) {
if(is.null(group)) {
if(is.vector(x)) { result <- sum(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/(sum(!is.na(x)) -lag )} else {
x <- as.matrix(x)
if(NCOL(x) == 1) {
  result <- sum(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/(sum(!is.na(x))-lag)
  } else {
n <- colSums(!is.na(x))  -lag
result <- colSums(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/n}
} } else {
x <- as.matrix(x) #added 26/5/14
if(NROW(group) != NROW(x)) group <- x[,group]  #added 26/5/14
nvar <- ncol(x)
cname <- colnames(x)
temp <- by(x,group, mssd,na.rm=na.rm,lag=lag)
 rownn <- lapply(temp,is.null)
               if(sum(as.integer(rownn)) > 0) {
               rown <-  names(temp)[-which(rownn==TRUE)] } else {rown <- names(temp) } 
   
result <- t(matrix(unlist(temp),nrow=nvar))
colnames(result) <- cname
rownames(result) <- rown
}
return(result)}


"rmssd" <- function(x,group=NULL,lag=1,na.rm=TRUE) {
return(sqrt(mssd(x,group=group,lag=lag,na.rm=na.rm))) }



##### Added March 1, 2017
"autoR" <- function(x,group=NULL,lag=1,na.rm=TRUE,use="pairwise") {

if(is.null(group)) {
n.obs <- NROW(x)
if(is.vector(x)) {
x <- as.vector(scale(x,scale=FALSE))
 mssd <- sum(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/(sum(!is.na(x))-lag)
 v1 <- sd(x[1:(n.obs-lag)],na.rm=na.rm)^2
 v2 <- sd(x[(lag+1):n.obs],na.rm=na.rm)^2
# r <- -(mssd - v1 - v2)/(2*sqrt(v1*v2))
r <- cor(x[1:(n.obs-lag)],x[(lag+1):n.obs],use=use)
 result <- list(autoR = r,rssd=sqrt(mssd))  #fixed May 10 ,2017 to correct autorR-
 } else {
x <- as.matrix(x)
n <- colSums(!is.na(x))
mssd <- colSums(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/(n-lag)
v1 <- apply(x[1:(n.obs-lag),],2,sd, na.rm=na.rm)^2
v2 <- apply(x[(lag+1):n.obs,],2, sd,na.rm=na.rm)^2
# r <- -(mssd - v1 - v2)/(2*sqrt(v1*v2))
r <- diag(cor(x[1:(n.obs-lag),],x[(lag+1):n.obs,],use=use))
 result <- list(autoR = r,rssd=sqrt(mssd))
}
} else {
 cl <- match.call()
x <- as.matrix(x) #added 26/5/14
if(NROW(group) != NROW(x)) group <- x[,group]  #added 26/5/14
nvar <- ncol(x)
cname <- colnames(x)
temp <- by(x,group, autoR,na.rm=na.rm,lag=lag)


 rownn <- lapply(temp,is.null)
               if(sum(as.integer(rownn)) > 0) {
               rown <-  names(temp)[-which(rownn==TRUE)] } else {rown <- names(temp) } 

tm <- t(matrix(unlist(temp),nrow=nvar*2))
autoR <- tm[,1:nvar]
rmssd  <- tm[,(nvar+1):(nvar*2)]

 colnames(autoR) <- colnames(rmssd) <- cname
 rownames(autoR) <- rownames(rmssd) <- rown

result <- list(autoR = autoR,rmssd=rmssd,Call=cl)
}
class(result) <- c("psych","autoR")
return(result)}



#####
sim.mssd <- function(n,r,g=.1) {
rw <- rnorm(n,sqrt(1/(1-r^2)))
x  <- xg <- rep(0,n)
for(i in 2:n) {x[i] <- r*x[i-1] + rw[i]
              xg[i] <- x[i] + g*i }
rx <- sample(x,n,replace=FALSE)
x2 <- x*2
rx2 <- rx*2
x.df <- data.frame(x,rx,x2,rx2,xg)
return(x.df)}


"build.html.help" <- function(p="psych",fn = "/Volumes/Test/psych/man",out="/Volumes/Test/help/") {        
    
     db <- list.files(fn)
     for (f in db) {tools::Rd2HTML(paste(fn,db[f]),out=paste(out,db[f]),package=p)
     } 
     }
        

    
#shannon complexity index
"shannon" <-  
   function(x,correct=FALSE,base=2) {if(is.null(dim(x))) {
        t <- table(x)
        s <- sum(t)
        p <- t/s
        H <- -sum(p * log(p,base))     
       if(correct) {
           Hmax <- -log(1/length(p),base)
           H <- H/Hmax}
   } else {  H <- apply(x,2,function(x) shannon(x, correct=correct, base=base))}      
     return(H)
}


test.all <- function(p) {
 library(p,character.only=TRUE)
  ob <- paste("package",p,sep=":")
  ol <- objects(ob)
  nf <- length(ol)
  for(i in 1:nf) {
    fn <- as.character(ol[[i]])
    example(topic=fn,package=p,character.only=TRUE)
    }
 detach(ob,character.only=TRUE)
}



  
  
  #lookup which x's are found in y[c1],return matches for y[]
"lookup" <- 
function(x,y,criteria=NULL) {
if (is.null(criteria)) {temp <- match(x,rownames(y))} else {
     temp <- match(x,y[,criteria])}
 y <- (y[temp[!is.na(temp)],,drop=FALSE])
return(y)}
 
 #use lookup to take fa/ic output and show the results 
 #modified July 4, 2017 to allow for omega output as well
"fa.lookup"  <-
   function(f,dictionary,digits=2) {
    f <- fa.sort(f)
   if(length(class(f)) > 1){ value <- class(f)[2] } else {value <- "none"}
   
   switch(value,
    
    omega = {f <- f$schmid$sl
            h2 <- NULL},
    fa    = {
             h2 <- f$communality
              com <- f$complexity
              f <- f$loading        },
    principal= {
                h2 <- f$communality
              com <- f$complexity
              f <- f$loading},
    iclust = {f <- f$loadings
              h2 <- NULL},
    none = {f <- f
            h2 <- NULL})
    
    

   ord <- rownames(f)
     
   contents <- lookup(rownames(f),dictionary)
   if(!is.null(h2)) {results <- data.frame(round(unclass(f),digits=digits),com=round(com,digits=digits),h2=round(h2,digits=digits))} else {
   results <- data.frame(round(unclass(f),digits=digits))}
   results <- merge(results,contents,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(results) <- results[,"Row.names"]
   results <- results[ord,-1]  #now put it back into the correct order
return(results)}
  



 
  "fa.organize" <- 
function(fa.results,o=NULL,i=NULL,cn=NULL) {
  if(!is.null(o)) {fa.results$loadings <- fa.results$loadings[,o]
       fa.results$Structure <- fa.results$Structure[,o]
       fa.results$weights <- fa.results$weights[,o]
       fa.results$valid <- fa.results$valid[o]
       fa.results$score.cor <- fa.results$score.cor[o,o]
       fa.results$r.scores <- fa.results$r.scores[o,o]
       fa.results$R2 <- fa.results$R2[o]
  if(!is.null(cn)) {colnames(fa.results$loadings) <- cn}
 fa.results$Phi <- fa.results$Phi[o,o]}
  if(!is.null(i)) {fa.results$loadings <- fa.results$loadings[i,]
   fa.results$Structure <- fa.results$Structure[i,]
       fa.results$weights <- fa.results$weights[i,]
       fa.results$complexity=fa.results$complexity[i]
       fa.results$uniquenesses <- fa.results$uniquenesses[i]}
  return(fa.results)
  }
  
  
  #fixed 13/6/14 to solve the missing data problem
  "con2cat" <- function(old,cuts=c(0,1,2,3),where) {
  new <- old
  nc <- length(cuts)
  if(missing(where)) where <- 1:ncol(old)
  lw <- length(where)
  if(is.matrix(cuts)) {mcuts <- cuts} else {mcuts <- matrix(rep(cuts,lw),nrow=lw,byrow=TRUE)}
  vwhere <- as.vector(where)
       for (w in 1:lw) {where <- vwhere[w]
       cuts <- mcuts[w,]
       nc <- length(cuts)
  if(nc < 2 ) {new[(!is.na(new[,where]) & ( old[,where] <= cuts)),where] <- 0
               new[(!is.na(new[,where]) & ( old[,where] > cuts)),where] <- 1}   else {
  
       
    
     new[(!is.na(new[,where]) & ( old[,where] <= cuts[1])),where] <- 0                
      for(i in(2:nc)) { 
           new[(!is.na(new[,where]) & ( old[,where] > cuts[i-1] )),where] <- i-1
          # & (new[(!is.na(new[,where]) & ( old[,where] > cuts[i-1] )),where]),where]  <- i-1
     }   
     new[(!is.na(new[,where]) & ( old[,where] > cuts[nc])),where]  <- nc }
     }
   new}
                         
   "keys.lookup" <- function(keys.list,dictionary) {
      if(is.list(keys.list)) { items <-  sub("-","",unlist(keys.list))
   f <- make.keys(items,keys.list)}
    keys.list <- fa.sort(f)
     contents <- lookup(rownames(f), y=dictionary)
    rownames(contents)[rowSums(f) <0 ] <- paste0(rownames(contents)[rowSums(f)<0],"-")
     return(contents)
    }
  
  "item.lookup" <- 
function (f,m, dictionary,cut=.3, digits = 2) {
    f <- fa.sort(f)
    if (!(is.matrix(f) || is.data.frame(f))) {
        h2 <- f$communality
        com <- f$complexity
        ord <- rownames(f$loadings)
        nfact <- ncol(f$loadings)
        f <- f$loadings
    }
    else {
        h2 <- NULL
        com <- NULL
        ord <- rownames(f)
        nfact <- ncol(f)
    }
    means <- m[ord]
    f <- data.frame(unclass(f),means=means)

    contents <- lookup(rownames(f), y=dictionary)
    if (!is.null(h2)) {
        results <- data.frame(round(f, digits = digits), 
            com = round(com, digits = digits), h2 = round(h2, 
                digits = digits))
    }
    else {
        results <- data.frame(round(f, digits = digits))
    }
    results <- merge(results, contents, by = "row.names", all.x = TRUE, 
        sort = FALSE)
    rownames(results) <- results[, "Row.names"]
    results <- results[ord, -1]
    res <- results[0,]  #make an empty data frame of the structure of results
    for (i in 1:nfact) { temp <-results[abs(results[,i]) > cut,]
     ord <- order(temp[,"means"])
     res <- rbind(res,temp[ord,])
     }
    return(res)
}



"falsePositive" <- function(sexy=.1,alpha=.05,power=.8) {
pf <- alpha * (sexy)
vp <- power * (1-sexy)
pf/(pf+vp)}



"bullseye" <- function(x,y,n) {
for(i in 1:n) {dia.ellipse(x,y,e.size=i)}}


"rel.val" <- function(n,sdx=.2,bars=TRUE,arrow.len=.05) {
if(n>20) {pc <- "."} else {pc <- 16}

plot(NA,xlim=c(0,10),ylim=c(0,10),axes=FALSE,xlab="",ylab="",main="Reliability and Validity as target shooting")
#Reliable and valid
x=3
y=2
bullseye(x,y,4)
x1 <- x + rnorm(n,0,sdx)
y1 <- y + rnorm(n,0,sdx)
xm <- mean(x1)
ym <- mean(y1)
points(x1,y1,pch=pc)
points(xm,ym,pch=20,col="red")
if(bars) error.crosses(x1,y1,add=TRUE,arrow.len=arrow.len,labels="")
text(x,y-2,"Reliable and valid")

#unReliable and  invalid
x=7
y=7
bullseye(x,y,4)
x1 <- x + rnorm(n,1,1)
y1 <- y + rnorm(n,1,1)
xm <- mean(x1)
ym <- mean(y1)
points(x1,y1,pch=pc)
points(xm,ym,pch=20,col="red")
if(bars) error.crosses(x1,y1,add=TRUE,arrow.len=arrow.len,labels="")
text(x,y-2,"Unreliable and Invalid")

#reliable and invalid
x=7
y=2
bullseye(x,y,4)
x1 <- x + rnorm(n,1,sdx)
y1 <- y + rnorm(n,1,sdx)
xm <- mean(x1)
ym <- mean(y1)
points(x1,y1,pch=pc)
points(xm,ym,pch=20,col="red")
if(bars)error.crosses(x1,y1,add=TRUE,arrow.len=arrow.len,labels="")

text(x,y-2,"Reliable and Invalid")


#unreliable, but valid
x=3
y=7

bullseye(x,y,4)
x1 <- x + rnorm(n,0,1)
y1 <- y + rnorm(n,0,1)
xm <- mean(x1)
ym <- mean(y1)
points(x1,y1,pch=pc)
points(xm,ym,pch=20,col="red")
if(bars) error.crosses(x1,y1,add=TRUE,arrow.len=arrow.len,labels="")
text(x,y-2,"Unreliable but Valid")
}
#rel.val(10,.5)


# "cor2" <- function(x,y,digits=2,use="pairwise",method="pearson") {
# R <- cor(x,y,use=use,method=method)
# print(round(R,digits))
# invisible(R)}

"cor2" <- function(x,y=NULL,digits=2,use="pairwise",method="pearson") {
multi <- FALSE
if(is.list(x) && is.null(y)) {multi <- TRUE
 n <- length(x)
xi <- x[[1]]
 for (i in 2:n) {xi <- cbind(xi,x[[i]])} 
R <- cor(xi,use=use,method=method) }else {
R <- cor(x,y,use=use,method=method)}
if(multi) {lowerMat(R,digits) } else {print(round(R,digits))}
invisible(R)}

levels2numeric <- function(x) {
n.var <- ncol(x)
for(item in 1:n.var) {
if (is.factor(x[,item])) x[,item] <- as.numeric(x[,item])}
invisible(x)
}


signifNum <-  function(x,digits=2) {
if(!is.null(ncol(x))) {sign <- rep(1,prod(dim(x)))} else {
      sign <- rep(1,length(x))}
sign[which(x < 0)] <- -1
base <- trunc(log10(sign*x))
mantissa <-  x/10^base
pretty <- round(mantissa,digits=digits-1) * 10^base 
pretty[which ((sign * x) == 0,arr.ind=TRUE)] <- 0 #fix the ones that are -Inf 
pretty}


#October 25, 2016
"char2numeric" <- function(x) {
 nvar <- NCOL(x)
 for(i in 1:nvar) {   
        if(!is.numeric(x[[i]] ))  {
                                  if(is.factor(unlist(x[[i]])) | is.character(unlist(x[[i]]))) {  x[[i]] <- as.numeric(x[[i]]) 
                          } else {x[[i]] <- NA} }
              } 
              invisible(x)} 

#this just shows if it is a matrix
"isCorrelation" <-  function(x) {value <- FALSE
  if(NROW(x) == NCOL(x)) {
  if( is.data.frame(x)) {if(isSymmetric(unname(as.matrix(x)))) { value <- TRUE}} else {if(isSymmetric(unname(x))) {value <- TRUE}}}
return(value)}


    
    
 

