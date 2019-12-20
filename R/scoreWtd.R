"scoreWtd" <- function(weights,items,std=TRUE,sums=FALSE,impute="none"){
 vars <-rownames(weights)
 n.scales <- NCOL(weights)
 vnames <- colnames(weights)
 if(any(c("(Intercept)","Intercept") %in% vars)) {items <- data.frame(Intercept=1,items)
   colnames(items)[1] <- "(Intercept)"
    vars[1] <- "(Intercept)"}
 selected <-items[vars]   #just use those items that have weights
  switch (impute,
  mean ={   miss <- which(is.na(items),arr.ind=TRUE)
       		item.means <- colMeans(items,na.rm=TRUE)   #replace missing values with means
       		items[miss]<- item.means[miss[,2]]},
   median={ miss <- which(is.na(items),arr.ind=TRUE) 	
       		item.med   <- apply(items,2,median,na.rm=TRUE) #replace missing with medians
        	items[miss]<- item.med[miss[,2]]}
  )
 if(std) {z.scores <-scale(selected)} else z.scores <- selected
# wtd.scores <-z.scores %*% (weights)  #this is the most basic version, but doesn't handle any missing

wtd.scores <- matrix(rep(NA,n.scales * NROW(items)),ncol=n.scales)  #this is just a dummy array

if(sums) {
weights <- t(weights)
for(j in 1:n.scales)  {wtd.scores[,j] <- colSums(weights[j,] *t(z.scores),na.rm=TRUE)}

} else {
if(n.scales ==1) { wtd.scores[,1] <- colMeans(weights[,1] *t(z.scores),na.rm=TRUE)} else {
weights <- t(weights)
for(j in 1:n.scales) {
 wtd.scores[,j] <-  colMeans(weights[j,] *t(z.scores),na.rm=TRUE)
}}
}
colnames(wtd.scores) <- vnames 
return(wtd.scores)
}
 
#developed September, 2019 to get more precise weights (still not beta weights) to take advantage of large sample stability 
