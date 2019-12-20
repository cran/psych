#Developed 12/1/18 to be part of bestScales approach
"manhattan" <- function(x,criteria=NULL,keys=NULL,raw=TRUE,n.obs=NULL,abs=TRUE,ylab=NULL,labels=NULL,log.p=FALSE,ci=.05,pch=21,main="Manhattan Plot of",adjust="holm",ylim=NULL,digits=2,dictionary=NULL, ...) {
 
 if(is.null(ylab) ) {if(log.p) {ylab = "- log(10) p of " } else {ylab <- "Correlations with " } 
                    }
                    
 #There are two basic cases:
 #1) a correlation matrix is provided, and the number of variables is nrow(x)   (raw == FALSE)
 #2) a data matrix is provided and the number of variables is ncol(x)   (raw == TRUE, the default)  
  pt.col <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 n.col <- 8
  nvar <- ncol(x)
 if(!is.null(dim(criteria))| (length(criteria) == NROW(x)))  { x <- cbind(x,criteria)
    criteria <- colnames(criteria)} else {
    if(length(criteria) > 1 ) {criteria <- (criteria) } else {criteria <- "criteria"}
    }
                  
if(raw) {n.obs <- nrow(x)
 #nvar <- ncol(x)
 #pt.col <- rainbow(n.col)  #here is one way of doing it
 #these next color are from GGplot2 and are said to color blind friendly


   
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
  }
 
 if(isCorrelation(x) | !raw) {r <- x      #  case 1
    raw <- FALSE
    nvar <- NROW(x)
    if(log.p) {if(is.null(n.obs)) {n.obs <- 1000
       message("\nn.obs not specified, arbitrarily set to 1000")}
     }
     if(!is.null(keys) ) {key.ord <- selectFromKeys(keys)
 num.keys <- length(keys)
 for (i in 1:num.keys) {
 select <- sub("-", "", unlist(keys[i]))
 keys[[i]]<- select}
 if(is.null(labels)) {labels <- names(keys)} else {if(labels=="none") {labels =rep("",num.keys)}}
} else {key.ord <- rownames(x)[1:nvar] }
     } else {  #case 2
    y <- x[,criteria]
    R <- corr.test(x,y,adjust=adjust,ci=FALSE)  #can we make this find polychorics and biserials
    r <- R$r
    colnames(r) <- criteria
     x <- as.matrix(x)
     raw <- TRUE}
  n.crit <- length(criteria)
  if(abs) r <- abs(r)
  if(log.p) {if(!raw) {r <- -log10(corr.p(r,n.obs,adjust=adjust,ci=FALSE)$p)}  else { r <-  - log10(R$p)}
   temp.p <- r
   temp.p[is.infinite(r)] <- NA
   min.p <- min(temp.p,na.rm=TRUE)
   r[is.infinite(r)] <-min.p }
  if(is.null(ylim)) ylim <- c(min(r[1:nvar,],na.rm=TRUE),max(r[1:nvar,],na.rm=TRUE))
  
  
 for(i in 1:n.crit) {
 if(is.null(keys)) {
 plot(r[1:nvar,i],main=paste(main,criteria[i]),ylab=paste(ylab,criteria[i]), ylim=ylim,pch=21,xlab="",xaxt="n",...)
 axis(1, at =1:nvar, labels = labels,las=2)
                     } else {
                     
                
plot(NA,xlim=c(1,num.keys),ylim=ylim,ylab=paste(ylab,criteria[i]),main= paste(main,criteria[i]),xlab="",xaxt="n",...)

axis(1, at =1:num.keys, labels = labels,las=2)
    for (k in 1:num.keys) {
    xloc <- rep(k,length(keys[[k]])) 
   # points(xloc,r[keys[[k]],i],pch=pch,bg=c("black","blue","red")[k %% 3 + 1],...)
   points(xloc,r[keys[[k]],i],pch=pch,bg=pt.col[k %% n.col + 1],...)
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
if(raw) {r <- r[1:nvar,]}   else {r <- r[1:nvar,criteria] }
 if(!is.null(dictionary)) {
   contents <- lookup(rownames(r),dictionary) 
   results <- data.frame(round(r,digits=digits))
   results <- merge(results,contents,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(results) <- results[,"Row.names"]
   results <- results[-1]  #now put it back into the correct order 
   r <- results    
    }
   
invisible(r)
}