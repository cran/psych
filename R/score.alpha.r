score.alpha <- function (keys,items,labels=NULL,totals=TRUE, digits=2) {
.Deprecated("score.alpha", msg = "score.alpha is deprecated.  Please use the scoreItems function")

    keys <- as.matrix(keys)   #just in case they were not matrices to start with
    items <- as.matrix(items)
    scores <- items %*%  keys  #this actually does all the work
    if (length(labels)>0) {colnames(scores) <- labels} #add labels
    abskeys <- abs(keys)
    item.var <- diag(var(items,use="pairwise"))  #find the item variances
    cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
    var.scales <- diag(cov.scales)
    cor.scales <- cor(scores,use="pairwise")     #could do this as matrix operation, but why bother
    sum.item.var <- item.var %*% abskeys  
    num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
    alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
   if (length(labels)>0) {colnames(alpha.scale) <- labels}
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   item.cor <- cor(items,scores,use="pairwise")
   if (!totals) scores <- scores/num.item   #find averages
    results <- list(scores=scores,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item, cor = round(cor.scales,digits), item.cor = round(item.cor,digits)) 
    class(results) <- "psych"
    return(results)
}
