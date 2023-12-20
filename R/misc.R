#A number of useful helper functions
#added January, 2012
#most are public, some are local just for me

#a dummy function to allow the help to find misc
"psych.misc" <-
function() {}

"lowerMat" <- 
function(R,digits=2,minlength=5) {
   lowleft <- lower.tri(R,diag=TRUE)
   nvar <- ncol(R)
	nc <- digits+3
	width <- getOption("width")
	k1 <- width/(nc+2)
   if(is.null(colnames(R))) {colnames(R) <- paste("C",1:nvar,sep="")}
   if(is.null(rownames(R))) {rownames(R) <- paste("R",1:nvar,sep="")}
	colnames(R) <- abbreviate(colnames(R),minlength=minlength)
	
	nvar <- ncol(R)
	nc <- digits+3
	#k1 <- width/(nc+2)
	
	if(k1 * nvar < width) {k1 <- nvar}  
	k1 <- floor(k1)
	if(!is.character(R) ) {
	fx <- format(round(R,digits=digits))} else {fx <- format(R)}
	 if(nrow(R) == ncol(R) ) {fx[!lowleft] <- ""}
	 for(k in seq(0,nvar,k1)) { if(k<nvar) {
	print(fx[(k+1):nvar,(k+1):min((k1+k),nvar)],quote=FALSE)}
	}
	invisible(R[lower.tri(R,diag=FALSE)])}
	
"lowerCor" <- 
function(x,digits=2,use="pairwise",method="pearson",minlength=5,cor="cor",show=TRUE) {
nvar <- NCOL(x)
x <- char2numeric(x)    #added 1/2/21
if(cor %in% c("tetrachoric","tetra","binary")) cor<- "tetrachoric"
if(cor %in% c("polychoric","poly")) cor<- "polychoric"
switch(cor, cor = {
            R <- cor(x, use = use, method = method)
        }, tetrachoric = {
            R <- tetrachoric(x)$rho
        }, polychoric = {
            R <- polychoric(x)$rho
        }, cov = {
            R <- cov(x, use = use)
        }, mixed = {
            R <- mixedCor(x)$rho
        })

  # R <- cor(x,use=use,method=method)
   if(show) lowerMat(R,digits,minlength=minlength)
   invisible(R)
   }

	
#adapted from utils::txtProgressBar
#modified August 10, 2012 to print just 100 times. 
"progressBar" <- 
function(value,max,label=NULL) {
if(inherits(stdout()[1],"terminal")) { pc <- round(100 * value/max)  #only print to the screen, not to a file
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
if(is.null(names(f))) {temp <- f
    f <-list()
    f$loadings <- temp}
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
v1 <- apply(x[1:(n.obs-lag),,drop=FALSE],2,sd, na.rm=na.rm)^2
v2 <- apply(x[(lag+1):n.obs,,drop=FALSE],2, sd,na.rm=na.rm)^2
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

#lookup a set of items from a bigger set
lookupItem <- function(x,y) {
n.look <- NROW(x)
possible <- list()
y <- as.character(y)
x <- as.character(x)
for(i in 1:n.look){
 temp <- grep(x[i],y)
 if(length(temp)>0) possible[i]<- temp

}
return(possible)
}
  
  
  #lookup which x's are found in y[c1],return matches for y[]
"lookup" <- 
function(x,y,criteria=NULL) {
if (is.null(criteria)) {temp <- match(x,rownames(y))} else {
     temp <- match(x,y[,criteria])}
     if(any(!is.na(temp))) {
 y <- (y[temp[!is.na(temp)],,drop=FALSE]) } else {y <- NA}
return(y)}
 
 #use lookup to take fa/ic output and show the results 
 #modified July 4, 2017 to allow for omega output as well
 #modified June 23, 2018 to limit to top n items per factor and abs(loading) > cut
 #modified April 30, 2020 to include the ability to handle bassAckward output
"fa.lookup"  <-
   function(f,dictionary=NULL,digits=2,cut=.0,n=NULL,sort=TRUE) {
    omega <- bassAck <- none <- lavaan <- NA    #weird requirement to avoid being global
   if(length(class(f)) > 1){ obnames <- cs(omega, fa, principal, iclust,bassAck, none)
     value <- inherits(f, obnames, which=TRUE)
			   if (any(value > 1)) {value <- obnames[which(value >0)]} else {value <- "none"}}
   if(sort & (value !="bassAck")) {f <- fa.sort(f)}
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
    bassAck = {nlevels <- length(f$bass.ack)	
              h2 <- NULL
              ord1 <- rownames(fa.sort(f$bass.ack[[nlevels-1]]))
               f <- f$bass.ack[[nlevels]]
               f <- f[,ord1]
               f <- fa.sort(f)
               colnames(f) <- paste0("F",1:ncol(f))           
    },
    none = {f <- f
            h2 <- NULL})
    
    n.fact <- NCOL(f)

   ord <- rownames(f)
   old.names <- ord
   ord <- sub("-","",ord)
   rownames(f) <- ord
   if(!is.null(dictionary))  { 
   contents <- lookup(rownames(f),dictionary)} else {message("fa.lookup requires a dictionary, otherwise just use fa.sort")}
   if(!is.null(h2)) {results <- data.frame(round(unclass(f),digits=digits),h2=round(h2,digits=digits),com=round(com,digits=digits))} else {
   results <- data.frame(round(unclass(f),digits=digits))}
   
   results <- merge(results,contents,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(results) <- results[,"Row.names"]
   results <- results[ord,-1]  #now put it back into the correct order
   rownames(results) <- old.names
   if(!is.null(n)) {
     rn <-rownames(results)
     results <- cbind(results,rn)
  	 f2c <- table(apply(abs(results[1:n.fact]),1,which.max))   #which column is the maximum value
   	 k <- 1
   	 j <- 1
   	 for(i in 1:n.fact) {
   		 results[k:(k+min(n,f2c[i])),] <-  results[j:(j+ min(n,f2c[i])),]
   		 k <- (k+min(n,f2c[i]))  
   		 j <- j + f2c[i]   }
    	results <- results[1:(k-1),]
    	rownames(results) <- results[,"rn"]
    	results <- results[,-NCOL(results)]
    }
    
if(cut > 0) {
   r.max <- apply(abs(results[,1:n.fact]),1,max)
   results <- results[abs(r.max) > cut,]
   }
return(results)}
  



 #revised 07/07/18 to add the cluster option
 #revised 12/07/18 to allow for simple matrix input
 #read a matrix, return a matrix
 #read a list, return a list
 #modified 9/2/20 for the case of not full rank cluster scores
  "fa.organize" <- 
function(fa.results,o=NULL,i=NULL,cn=NULL,echelon=TRUE,flip=TRUE) {
if(is.null(names(fa.results)) )  {temp <- fa.results   #the matrix form          
                 if(flip) {
                 total.load <-colSums(temp)
                 flipper <- sign(total.load)
                 flipper[flipper==0] <-1 
                 temp <- t( t(temp) * flipper ) }
                 if(!is.null(o)) {temp <- temp[,o]}
                 if(!is.null(i)) {temp <-temp[i,]}
                 fa.results <- temp 
                 nf <- ncol(temp)} else { # the list form 
    nf <- ncol(fa.results$loadings)
     if(echelon & is.null(o) ) {temp <- apply(abs(  fa.results$loadings),1,which.max)
     nf <- ncol(fa.results$loadings)
     nvar <- nrow(fa.results$loadings)
  o <- rep(NA,nf)
  k <- 1
  o[k] <- temp[k]
  for (ki in 2:nvar) {if (!(temp[ki] %in% o[1:k])) {o[k+1] <- temp[ki]
    k <- k + 1}  }
    }
    
    #now, a kludge for the case where there are some empty columns
   
    n.good <- sum(1:nf %in% o)
    if (n.good < nf) { tt <- 1:nf
      o[(n.good + 1):nf] <- tt[!1:nf %in% o]
       }
  if(flip) {
        total.load <- colSums(fa.results$loadings)
       flipper <- sign(total.load)
        flipper[flipper==0] <-1 } else { flipper <- rep(1,NCOL(fa.results$loadings)) }
   fa.results$loadings <- t(t(fa.results$loadings) * flipper)                      


  if(!is.null(o)) {fa.results$loadings <- fa.results$loadings[,o]
      flipper <- flipper[o] 
       fa.results$Structure <- t(t(fa.results$Structure[,o]) * flipper)
       fa.results$valid <- t(t(fa.results$valid[o])*flipper)
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
       }
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
   # f <- fa.sort(f)
    none<- NULL   #A strange requirement of R 4.0
     if(length(class(f)) > 1){ obnames <- cs(omega, fa, principal, iclust, none)
     value <- inherits(f, obnames, which=TRUE)
			   if (any(value > 1)) {value <- obnames[which(value >0)]} else {value <- "none"}
			   f <- fa.sort(f) } else {value <- "none"}
  old.names <- NULL
   switch(value,
    
    omega = {f <- f$schmid$sl
            h2 <- NULL
            old.names <- rownames(f)
            rownames(f) <- sub("-","",old.names)},
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
    means <- m[ord]   #incorrectly added a comma to allow it sort dataframes 6/22/21
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
    if(!is.null(old.names)) rownames(results) <- old.names
    res <- results
  #   res <- results[0,]  #make an empty data frame of the structure of results
#     for (i in 1:nfact) { temp <- results[abs(results[,i]) > cut,]
#      ord <- order(temp[,"means"])
#      res <- rbind(res,temp[ord,])
#      }
    return(res)
}

lmCorLookup <- setCorLookup <- function(x,dictionary=NULL,cut=0,digits=2,p=.05) {
coef <- x$coefficients
probs <- x$Probability
labels <- dictionary[rownames(coef),,drop=FALSE]

coef[probs > p] <- NA
result <- list()
nvar <- NCOL(coef)
for (i in 1:nvar) {
ord <- order(abs(coef[,i]),decreasing=TRUE)
temp <- cbind(coef=round(coef[ord,i],digits),labels[ord,])
result[[i]] <- data.frame(temp[!is.na(temp[,1]),])
}
names(result) <- colnames(coef)
result
}


"lookupItems" <- function(content=NULL,dictionary=NULL,search=c("Item","Content","item")) {
  location <-which(colnames(dictionary) %in% search)
  value <- grep(content,dictionary[,location])
  dictionary[value,,drop=FALSE]
}


"falsePositive" <- function(sexy=.1,alpha=.05,power=.8) {
pf <- alpha * (sexy)
vp <- power * (1-sexy)
pf/(pf+vp)}

"build.html.help" <- function(p="psych",fn = "/Volumes/Test/psych/man",out="/Volumes/Test/help/") {        
    
     db <- list.files(fn)
     for (f in db) {tools::Rd2HTML(paste(fn,db[f]),out=paste(out,db[f]),package=p)
     } 
     }
        


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

#finally added the char2numeric so it will not choke on character variables   1/3/21

#added the cor option 11/18/23
"cor2" <- function(x,y=NULL,digits=2,use="pairwise",method="pearson",cor="cor") {
multi <- FALSE
if(is.list(x) && is.null(y)) {multi <- TRUE
 n <- length(x)
xi <- x[[1]]
 for (i in 2:n) {xi <- cbind(xi,x[[i]])}

switch(cor,
cor= {R <- cor(xi,use=use,method=method) },
tetrachoric = {R <- tetrachoric(xi)$rho},
polychoric= {R <- polychoric(xi)$rho},
cov = {R <- cov(xi, use=use)},
mixed = {R <- mixedCor(xi)$rho}
)
}else {
x <- char2numeric(x)
y <- char2numeric(y)

switch(cor,

cor= {R <- cor(x,y,use=use,method=method)},
tetrachoric = {R <- tetrachoric(x,y)$rho},
polychoric= {R <- polychoric(x,y)$rho},
cov = {R <- cov(x,y, use=use)},
mixed = {R <- mixedCor(cbind(x,y))$rho}

)
}
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
#June 18 2020  added  as.factor to convert character strings that are not stored as factors
#this has a downsize that it converts numbers stored as characters to factors (see nchar2numeric)
#added the flag option and the change the colname option 1/2/21

"isAllNumeric" <- function(x){
      nvar <- NCOL(x)
      bad <- rep(0,nvar)
 for(i in 1:nvar) { 
     if(!is.numeric(x[[i]]) ){if(is.factor(unlist(x[[i]]))| is.character(unlist(x[[i]]))    
     ) bad [i] <- 1
     }}
     if(prod(bad)==0) {result <- TRUE} else {result<- bad}
     invisible(result)
     }

"char2numeric" <- function(x,flag=TRUE) {

 nvar <- NCOL(x)
 for(i in 1:nvar) {   
        if(!is.numeric(x[[i]] ))  {
                                  if(is.factor(unlist(x[[i]])) | is.character(unlist(x[[i]]))) {  x[[i]] <- as.numeric(as.factor(x[[i]])) 
                          } else {x[[i]] <- NA}
                         if(flag) colnames(x)[i] <- paste0(colnames(x)[i],"*")
                          } 
                          
              } 
              invisible(x)} 
              
              
#added June 25, 2020 to handle the case of numeric data stored as characters  
#added flag and colnames option  1/2/21             
"nchar2numeric" <- function(x,flag=TRUE) {
 nvar <- NCOL(x)
 for(i in 1:nvar) {   
        if(!is.numeric(x[[i]] ))  { 
                if(is.factor(unlist(x[[i]])) | is.character(unlist(x[[i]]))) {
                                  if(is.factor(unlist(x[[i]]))) {  x[[i]] <- as.numeric(as.factor(x[[i]])) } else {x[[i]] <- as.numeric(x[[i]])}
                          } else {x[[i]] <- NA}
                          if(flag) colnames(x)[i] <- paste0(colnames(x)[i],"*")
                            }
              } 
              invisible(x)} 


#this just shows if it is a matrix is symmetric and has diagonals of 1
#Added the unclass to  handle a problem with class partial.r  4/10/21
"isCorrelation" <-  function(x,na.rm=FALSE) {value <- FALSE
  if(NROW(x) == NCOL(x)) {
  if( is.data.frame(x)) {if(isSymmetric(unclass(unname(as.matrix(x))))) { value <- TRUE}} else {if(isSymmetric(unclass(unname(x)))) {value <- TRUE}}
  value <- value && isTRUE(all.equal(prod(diag(as.matrix(x))),1) )
  if(!value  && (na.rm) && any(is.na(diag(as.matrix(x))))) stop("Although the matrix is symmetric, one of the elements of the  diagonal is NA. Check your data.")
  value <- value && isTRUE((min(x,na.rm=TRUE)>= -1) & (max(x,na.rm=TRUE) <= 1))  
  }
  return(value)}
  
  #this just shows if it is a symmetric matrix
"isCovariance" <-  function(x) {value <- FALSE
  if(NROW(x) == NCOL(x)) {
  if( is.data.frame(x)) {if(isSymmetric(unclass(unname(as.matrix(x))))) { value <- TRUE}} else {if(isSymmetric(unclass(unname(x)))) {value <- TRUE}}}
 # value <- value && isTRUE(all.equal(prod(diag(as.matrix(x))),1) )  #don't check for diagonal of 1
  return(value)}
  
  


#cs is taken from Hmisc:::Cs
cs <- function(...) {as.character(sys.call())[-1]}
#acs is modified to produce a single string
acs <- function(...) {gsub(",","",toString(sys.call()[-1]))}

fromTo <- function(data,from,to=NULL) {cn <- colnames(data)
    if(is.null(to)) {to <- from[2]
       from <- from[1]}
  from <- which(cn == as.character(from))
 to =  which(cn == as.character(to))
  select <- from:to
  return(data[select])}
  
  

#flip -- not public but useful trick
flip <- function(R,key) {#reverse correlations with keys  < 0
	R[key<0,] <- -R[key < 0,]
	R[,key<0] <- -R[,key < 0]
  return(R)}
    
#do matrix multiplication with missing values
 matMult <- function(x,y) {  #matrix multiplication without matrices!  
  nvar1 <- NCOL(x)
  nvar2 <- NROW(y)
  M <- matrix(NA,ncol =NCOL(y),nrow=NROW(x))
    if(NCOL(x) != NROW(y)) {stop("incompatible dimensions")}
   for(i in 1:NROW(x)) {
    for (j in 1:NCOL(y))  {
    M[i,j] <-  sum(x[i,] * y[,j], na.rm=TRUE)
    } 
    }
    return(M)
    }
     
   
   
   factorScoresSapa  <- function(weights,items) {
  nf <- NCOL(weights)
  scores<- matrix(NA,nrow=NROW(items),ncol=nf)
  n.obs <- NROW(items)
  n.items <- NCOL(items)
  titem <- t(items)
  for (factors in (1:nf)) {
           
    scores[,factors] <-t(colMeans( weights[,factors] * titem,  na.rm=TRUE))
    }
    scores<- scale(scores)  #return standardized values
    return(scores)
  }

countBad <- function(r.list) {
n <- length(r.list)
count <- 0
for (i in 1:n) {
 ev <-eigen(r.list[[i]])
 bad <- any(ev$values < 0)
 count <- count + bad}
 return(count)
 }

#Create SAPA data by sampling
SAPAfy <- function(x,y)
{
  n_vars <- ncol(x)
  n_questions_unanswered <- round(n_vars - y)  #  in case someone inputs a double.
  
  for (xxx in 1:nrow(x)) 
  {
    x[xxx, sample.int(n_vars, n_questions_unanswered)  ] <- NA 
  }
  return(x)
}


