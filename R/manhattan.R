#Developed 12/1/18 to be part of bestScales approach
"manhattan" <- function(x,criteria=NULL,keys=NULL,abs=TRUE,ylab=NULL,labels=NULL,log.p=FALSE,ci=.05,pch=21,main="Manhattan Plot of",adjust="holm",ylim=NULL, ...) {
 
 if(is.null(ylab) ) {if(log.p) {ylab = "- log(10) p of " } else {ylab <- "Correlations with " } 
                    }
 n.obs <- nrow(x)
 nvar <- ncol(x)
 
 if(!is.null(dim(criteria))| (length(criteria) == NROW(x)))  { x <- cbind(x,criteria)
    criteria <- colnames(criteria)} else {
    if(length(criteria) > 1 ) {criteria <- (criteria) } else {criteria <- "criteria"}
    }
    
 if(!is.null(keys) ) {key.ord <- selectFromKeys(keys)
 num.keys <- length(keys)
 for (i in 1:num.keys) {
 select <- sub("-", "", unlist(keys[i]))
 keys[[i]]<- select}
 if(is.null(labels)) labels <- names(keys)
} else {key.ord <- colnames(x)[1:nvar] }
 if(is.null(dim(criteria)) ) {x <- x[c(key.ord,criteria)] } else {x <- x[key.ord]}
 nvar <- length(key.ord)
# if(!is.null(dim(criteria))| (length(criteria) == NROW(x)))  { x <- cbind(x,criteria)
 #   if(length(criteria) > 1 ){criteria <- names(criteria) } else {criteria <- "criteria"}
 #   }
 n.crit <- length(criteria) 
  
 
 if(isCorrelation(x)) {r <- x      #  case 1
    raw <- FALSE} else {  #case 2
    y <- x[,criteria]
    R <- corr.test(x,y,adjust=adjust)
    r <- R$r
    colnames(r) <- criteria
     x <- as.matrix(x)
     raw <- TRUE}
  n.crit <- length(criteria)
  if(abs) r <- abs(r)
  if(log.p) r <-  - log10(R$p)
  if(is.null(ylim)) ylim <- c(min(r[1:nvar,]),max(r[1:nvar,]))
  
 for(i in 1:n.crit) {
 if(is.null(keys)) {
 plot(r[1:nvar,i],main=paste(main,criteria[i]),ylab=paste(ylab,criteria[i]), ylim=ylim,pch=21,xlab="",xaxt="n",...)
 axis(1, at =1:nvar, labels = labels,las=2)
                     } else {
                     
                
  plot(NA,xlim=c(1,num.keys),ylim=ylim,ylab=paste(ylab,criteria[i]),main= paste(main,criteria[i]),xlab="",xaxt="n",...)

axis(1, at =1:num.keys, labels = labels,las=2)
    for (k in 1:num.keys) {
    xloc <- rep(k,length(keys[[k]])) 
    points(xloc,r[keys[[k]],i],pch=pch,bg=c("black","blue","red")[k %% 3 + 1],...)
   }
   }
if(!is.null(ci))  {  dif.corrected <- qnorm(1 - ci/(nvar * n.crit))
                 se <-  1/sqrt(n.obs - 3)
                 upper <- fisherz2r(dif.corrected * se)
                 abline(h=upper,lty="dashed")
                 if(!abs) {abline(h=0)
                 abline(h=-upper,lty="dashed")}


        }
       
    }
invisible(r)
}