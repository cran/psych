"score.items"  <- 
 function (keys,items,totals=TRUE,ilabels=NULL, missing=TRUE, min=NULL,max=NULL,digits=2) {
    keys <- as.matrix(keys)   #just in case they were not matrices to start with
    items <- as.matrix(items)
    item.means <- colMeans(items,na.rm=TRUE)
    if (is.null(min)) {min <- min(items,na.rm=TRUE)}
    if (is.null(max)) {max <- max(items,na.rm=TRUE)}
    if(missing) { miss.rep <- rowSums(is.na(items))
       item.means <- colMeans(items,na.rm=TRUE)
            miss <- which(is.na(items),arr.ind=TRUE)
           items[miss]<- item.means[miss[,2]]}
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
    n.subjects <- dim(items)[1]
     scores<- items %*%  keys  #this actually does all the work but doesn't handle missing
    
    scores<- items %*%  keys  #this actually does all the work
    slabels <- colnames(keys)
    if (is.null(slabels)) {slabels<- paste("S",1:n.keys,sep="")}
    colnames(scores) <- slabels  
    
    abskeys <- abs(keys)
    item.var <- diag(var(items,use="pairwise"))  #find the item variances
    cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
    var.scales <- diag(cov.scales)
    cor.scales <- cor(scores,use="pairwise")     #could do this as matrix operation, but why bother
    sum.item.var <- item.var %*% abskeys  
    num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
    alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
    colnames(alpha.scale) <- slabels
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   item.cor <- cor(items,scores,use="pairwise")
   if(is.null(ilabels)) {ilabels <-  paste("I",1:n.items,sep="")}
   rownames(item.cor) <- ilabels
   correction <- (colSums(abs(keys)-(keys))/2)*(max+min)
   scores <- scores  + matrix(rep(correction,n.subjects),byrow=TRUE,nrow=n.subjects)
   
   if (!totals) {scores <- scores %*% diag(1/num.item)   #find averages
                 colnames(scores) <- paste("A",1:n.keys,sep="") }
    scale.cor <- correct.cor(cor.scales,t(alpha.scale))
   if (missing){results <-list(scores=scores,missing = miss.rep,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits),cor = round(cor.scales,digits) ,corrected = round(scale.cor,digits))} else{  
    results <- list(scores=scores,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits), cor = round(cor.scales,digits),corrected = round(scale.cor,digits))} 
 }
 
