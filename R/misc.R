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

  
  "best.items" <- 
function(x,c1=1,cut=.3, abs=TRUE, contents=NULL,digits=2) {
if(abs) {ord <- order(abs(x[,c1]),decreasing=TRUE)
  value <- x[ord,c1]
  value <- value[abs(value) >cut]
  } else {ord <- order(x[,c1],decreasing=TRUE)
  value <- x[ord,c1]
  value <- value[abs(value) >cut] }
value <- round(data.frame(value),digits)
if(!is.null(contents)) {if(!is.factor(contents)) {temp <- lookup(rownames(value),contents,c1=1)
  if(nrow(value) > nrow(temp))  value <- value[rownames(value) %in% contents[,1], ]
  value <- data.frame(value,temp)}}
return(value)
}
  
  
  #lookup which x's are found in y[c1],return matches for y[]
 "lookup" <- 
function(x,y,c1=1) {
temp <- match(x,y[,c1])
 y <- (y[temp[!is.na(temp)],])
  return(y)}


