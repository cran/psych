"score.items"  <- 
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE,impute="median",  min=NULL,max=NULL,digits=2,short=FALSE) {
    keys <- as.matrix(keys)   #just in case they were not matrices to start with
    items <- as.matrix(items)
    item.means <- colMeans(items,na.rm=TRUE)
    if (is.null(min)) {min <- min(items,na.rm=TRUE)}
    if (is.null(max)) {max <- max(items,na.rm=TRUE)}
     miss.rep <- rowSums(is.na(items))
    if(missing) {
        miss <- which(is.na(items),arr.ind=TRUE)
        if(impute=="mean") {
       item.means <- colMeans(items,na.rm=TRUE)   #replace missing values with means
       items[miss]<- item.means[miss[,2]]} else {
       item.med   <- apply(items,2,median,na.rm=TRUE) #replace missing with medians
        items[miss]<- item.med[miss[,2]]}
        }
           
                  
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
    n.subjects <- dim(items)[1]
     scores<- items %*%  keys  #this actually does all the work but doesn't handle missing values
    
    slabels <- colnames(keys)
    if (is.null(slabels)) {
    	if (totals) {slabels<- paste("S",1:n.keys,sep="")} else {
    	             slabels <- paste("A",1:n.keys,sep="")} }
   
    
    abskeys <- abs(keys)
    item.var <- diag(var(items,use="pairwise"))  #find the item variances
    cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
    var.scales <- diag(cov.scales)
    cor.scales <- cor(scores,use="pairwise")     #could do this as matrix operation, but why bother
    sum.item.var <- item.var %*% abskeys  
    num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
   alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
    colnames(alpha.scale) <- slabels
  alpha.scale[is.nan(alpha.scale)] <- 1
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   
   item.cor <- cor(items,scores,use="pairwise")    #this does not correct for item overlap
   if(is.null(ilabels)) {ilabels <-  paste("I",1:n.items,sep="")}
   rownames(item.cor) <- ilabels
   correction <- (colSums(abs(keys)-(keys))/2)*(max+min) #correct for flipping
   scores <- scores  + matrix(rep(correction,n.subjects),byrow=TRUE,nrow=n.subjects)
  
   
   if (!totals) {
    if(n.keys > 1) {scores <- scores %*% diag(1/num.item)   #find averages
                  }  else
                 {scores <- scores/num.item } }
    scale.cor <- correct.cor(cor.scales,t(alpha.scale))
     colnames(scores) <- slabels
   if (short) {results <- list(alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits),cor = round(cor.scales,digits) ,corrected = round(scale.cor,digits))} else {
   if (sum(miss.rep) >0) {results <-list(scores=scores,missing = miss.rep,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits),cor = round(cor.scales,digits) ,corrected = round(scale.cor,digits))} else{  
    results <- list(scores=scores,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits), cor = round(cor.scales,digits),corrected = round(scale.cor,digits))} }
    
    return(results)
 }
 #modified June 1 to add row names to items 
 #modified June 22 to add median imputation
 #modified August 8 to add colnames to scores
 #modified Sept 23, 2007 to allow for short output
