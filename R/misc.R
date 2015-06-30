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
	}}
	
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
if(class(stdout())[1]=="terminal")  #only print to the screen, not to a file
 { pc <- round(100 * value/max)
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


#developed January 10, 2012 to find the mean square of successive differences
#see Von Neuman et al. 1941
"mssd" <- function(x,group=NULL,lag=1,na.rm=TRUE) {
if(is.null(group)) {
if(is.vector(x)) { result <- sum(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/(length(x)-1)} else {
x <- as.matrix(x)
n <- colSums(!is.na(x)) -1
result <- colSums(diff(x,lag=lag,na.rm=na.rm)^2,na.rm=na.rm)/n}
} else {
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
return(sqrt(mssd(x,lag=lag,group=group,na.rm=na.rm))) }


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



  "bestItems" <- 
function(x,criteria=1,cut=.3, abs=TRUE, dictionary=NULL,cor=TRUE,digits=2) {
if((nrow(x) !=ncol(x)) && cor) {x <- cor(x,use="pairwise")} #convert to correlation if necessary
if(abs) {ord <- order(abs(x[,criteria]),decreasing=TRUE)
  value <- x[ord,criteria,drop=FALSE]
  value <- value[(abs(value) >cut),,drop=FALSE]
  } else {ord <- order(x[,criteria],decreasing=TRUE)
  value <- x[ord,criteria]
  value <- value[abs(value) > cut] }
value <- round(data.frame(value),digits)
if((!is.null(dictionary)) && !is.factor(dictionary)) {temp <- lookup(rownames(value),dictionary)
   value <- merge(value,temp,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(value) <- value[,"Row.names"]
   value <- value[-1]
  if(abs) {ord <- order(abs(value[,criteria]),decreasing=TRUE) } else {ord <- order(value[,criteria],decreasing=TRUE)}
   value <- value[ord,] 
   }
return(value)
}
  
  
  #lookup which x's are found in y[c1],return matches for y[]
 "lookup" <- 
function(x,y,criteria=NULL) {
if (is.null(criteria)) {temp <- match(x,rownames(y))} else {
     temp <- match(x,y[,criteria])}
 y <- (y[temp[!is.na(temp)],,drop=FALSE])
  return(y)}
 
 #use lookup to take fa/ic output and show the results 
"fa.lookup"  <-
   function(f,dictionary,digits=2) {
   f <- fa.sort(f)
  
   if(!(is.matrix(f) || is.data.frame(f))) {h2 <- f$communality
    com <- f$complexity
     ord <- rownames(f$loadings) #keep this order
    f <- f$loadings} else {h2<- NULL
       com <- NULL
       ord <- rownames(f)}
    
   contents <- lookup(rownames(f),dictionary)
   if(!is.null(h2)) {results <- data.frame(round(unclass(f),digits=digits),com=round(com,digits=digits),h2=round(h2,digits=digits))} else {
   results <- data.frame(round(unclass(f),digits=digits))}
   results <- merge(results,contents,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(results) <- results[,"Row.names"]
   results <- results[ord,-1]  #now put it back into the correct order
   return(results)}
  


 #created 20/2/14
 #find the scales based upon the items that most correlate with a criteria
 #pure dust bowl empiricism
 #modified 13/3/15 to handle the problem of missing item labels
"bestScales" <- 
 function(x,criteria,cut=.1,n.item =10, overlap=FALSE,dictionary=NULL,digits=2) {

#first, declare a function to identify the bad items and drop them from the keys
 findBad <- function(key,r) { 
ss <- abs(key) > 0 
rss <- r[ss,ss]
if(any(is.na(rss))){ #some of these are bad
n.bad <-  apply(rss,1,function(x) sum(is.na(x)))
key[names(which.max(n.bad))] <- 0
findBad(key,r)}
return(key)
}

short <- function(key,r) { 
 kn <- names(key[abs(key[,1]) >0,1])
 if(is.null(kn)) kn <- names(which(abs(key[,1]) >0))
 cn <- colnames(key)
 ord <- order(abs(r[kn,cn]),decreasing=TRUE)
 kn <- kn[ord]
 result <- r[kn,cn,drop=FALSE]
 return(result)
}
#begin the main function
 nvar <- ncol(x)
 if(nrow(x) != nvar) {r <- cor(x,use="pairwise")} else {r <- x}  #convert data to a correlation matrix  #don't actually need to have  a square matrix
 ny <- length(criteria)
 nc <- length(cut)
 ni <- length(n.item)
 ord.name <- NULL
if(length(cut) == 1)  cut <- rep(cut,ny)
if(length(n.item) == 1) n.item <- rep(n.item,ny)
 if(ny > 1 ) {ord <- apply(abs(r[,criteria]),2,order,decreasing=TRUE) 
     for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[criteria[i],ord[n.item[i]+1,i]])) 
     ord.name <- c(ord.name, rownames(r)[ord[1:n.item[i],i]] )
    }
     } else {
         ord <- order(abs(r[criteria,]),decreasing=TRUE)
         for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[ord[n.item[i]+1],criteria])) }
        }

 key <- matrix(0,ncol=ny,nrow=nvar)
 key[t(r[criteria,] >= cut)] <- 1
 key[t(r[criteria,] <= -cut)] <- -1
 rownames(key)  <- colnames(r)
 colnames(key)  <- criteria
 c <- key  #this just gets it to be a matrix of the right size and names
 if(!overlap)  {key[criteria,criteria] <- 0} else {for(i in 1:ny) key[criteria[i],criteria[i]] <- 0}
 #colnames(key) <- paste(criteria,"S",sep=".")
 colnames(key) <- criteria
 
if(any(is.na(r))) {#Are there any bad values
  for(i in 1:ny) {#key[,i] <- findBad(key[,i],r)  #Drop the bad items from any scoring key
  c[,i] <- colSums(key[,i] * r,na.rm=TRUE)}    #replace matrix addition with a colSums
  c <- t(c)
} else {#otherwise, don't bother


 c<- t(key) %*% r    #we can do the matrix multiply because there are no bad data         
 }

 C <- c %*% key
 if(ny < 2) {re <- c[,criteria]/sqrt(C) } else {
 re <- diag(c[,criteria])/sqrt(diag(C))}
ni <- colSums(abs(key))
R <- cov2cor(C)

short.key <- list()
value <- list()


for(i in 1:ny) {short.key[[criteria[i]]] <- short(key[,i,drop=FALSE],r) 

if(!is.null(dictionary)) {if(!is.factor(dictionary)) {temp <- lookup(rownames(short.key[[criteria[i]]]),dictionary)

  value[[criteria[[i]]]] <- merge(short.key[[i]],temp,by="row.names",all.x=TRUE,sort=FALSE)
 ord <- order(abs(value[[criteria[[i]]]][[criteria[[i]]]]),decreasing=TRUE)
  value[[criteria[[i]]]] <- value[[criteria[[i]]]][ord,]
 } 
 }}
results <- list(r=re,n.items=ni,R=R,cut=cut,short.key=short.key,value=value,key=key,ordered=ord.name)
class(results) <- c("psych","bestScales")
return(results)
}
 
 
 
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
                         
   
  
  "item.lookup" <- 
function (f,m, dictionary,cut=.3, digits = 2) 
{
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
