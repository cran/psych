#Developed 06/24/21 to find mean item validities
# x is the set of items
# criteria is a single or vector of variables to correlate with the items
# keys is a keys.list of scoring keys used to find the scales

validityItem <- item.validity <- function(x,criteria,keys)  {
#first find the correlation of the items with the criteria

n.keys <- length(keys)
nvar <- NCOL(criteria)
key.n <- names(keys)
select <- lapply(keys,selectFromKeys)
valid <- matrix(NA,nrow=n.keys,ncol=nvar)
for(i in 1:n.keys) {
   list.i <-  unlist(keys[[i]])
        pos <- rep(1,length(keys[[i]]))
		neg <- grep("-",list.i)
		if(!is.null(neg)) pos [neg] <- -1
		
if(NCOL(criteria)==1) {valid[i,1] <- mean(cor(x[,select[[i]]],criteria,use="pairwise")*pos,na.rm=TRUE)   
                       } else {
for (j in 1:NCOL(criteria))  {	
 valid[i,j] <- mean(cor(x[,select[[i]]],criteria[j],use="pairwise")*pos,na.rm=TRUE)   
 } 
  }
  }
rownames(valid)<- key.n
colnames(valid) <- colnames(criteria)
return(valid)
}


#this combines reliability and predicted validity to produce interesting results
predicted.validity <- function(x,criteria,keys,scale.rel=NULL,item.val =NULL) {
 cl <- match.call()
n.keys <- length(keys)
asymp <- predicted <- matrix(NA,ncol=NCOL(criteria),nrow=n.keys)
#If we provide these, we don't need to find them
if(is.null(item.val)) item.val <-  validityItem(x,criteria,keys)
if(is.null(scale.rel)) scale.rel <- reliability(keys,x)

for(i in 1:n.keys) {
n.item <- scale.rel$result.df[i,"n.items"]
predicted[i,] <- item.val[i,] * n.item /sqrt(n.item + n.item * (n.item -1) * scale.rel$result.df[i,"mean.r"])
asymp[i, ] <- item.val[i,] /sqrt(scale.rel$result.df[i,"mean.r"])
}
rownames(predicted) <- rownames(asymp) <- names(keys)
colnames(predicted) <- colnames(asymp) <-colnames(item.val)
result <- list(predicted=predicted,item.validities=item.val,scale.reliabilities = scale.rel, asymptotic = asymp, Call=cl)
class(result) <- c("psych","validity")
return(result)
}

cd.validity<- function(d,keys,abs=TRUE)  {
#cd must be returned from Cohen.d
if(inherits(d, "cohen.d")) { cd <- d$cohen.d[,"effect",drop=FALSE]} else {cd <- as.matrix(d)}
#first find the correlation of the items with the criteria
#cd <- cd$cohen.d[,"effect",drop=FALSE]
n.keys <- length(keys)

key.n <- names(keys)
select <- lapply(keys,selectFromKeys)
valid <- matrix(NA,nrow=n.keys,ncol=1)
for(i in 1:n.keys) {
   list.i <-  unlist(keys[[i]])
        pos <- rep(1,length(keys[[i]]))
        
		neg <- grep("-",list.i)
		if(!is.null(neg)) pos [neg] <- -1
		
 if(all(select[[i]] %in% rownames(cd))) {valid[i,] <- mean(cd[select[[i]],"effect"]*pos,na.rm=TRUE) } else {valid[i] <- NA}  
 }
rownames(valid)<- key.n
#colnames(valid) <- colnames(criteria)
return(valid)
}