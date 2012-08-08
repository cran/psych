 "score.multiple.choice" <-  
  function(key,data,score=TRUE,totals=FALSE,ilabels=NULL, missing=TRUE,impute="median", digits=2,short=TRUE) {
  #convert a data matrix or data with multiple choice responses to correct/incorrect
     cl <- match.call()
  if(!is.matrix(data)) {if(!is.data.frame(data)) {stop("data must be either a data frame or matrix!")} else data <- as.matrix(data)}
  nvar <- dim(data)[2]
  
  response.freq <- response.frequencies(data)
  alternatives <- dim(response.freq)[2]
  
  if(length(key)==nvar) {
      items <- t(t(data)==key[]) #scores as t/f
      items <- items + 0    #converts t/f to 1/0 } 
     }   else {stop("key must have as many elements as columns of 'data' ")}
 
if (score) {
  item.stats <- describe(items,ranges=FALSE)[,2:7]
  miss.rep <- rowSums(is.na(items))
    if(missing) {
        miss <- which(is.na(items),arr.ind=TRUE)
        if(impute=="mean") {
       item.means <- colMeans(items,na.rm=TRUE)     #replace missing with means
       items[miss]<- item.means[miss[,2]]} else {
       item.med   <- apply(items,2,median,na.rm=TRUE)  #or medians
       items[miss]<- item.med[miss[,2]]}
        } 

keys <- rep(1,nvar)      #now, score the items as the sum of correct
                      
 scores <- rowSums(items,na.rm=TRUE)
 slabels <- colnames(keys)
 if (is.null(slabels)) {
    	if (totals) {slabels<- paste("Totals") } else {
    	             slabels <- paste("Averages")} }
   names(scores) <- slabels

  r.items <- cov(items,use="pairwise")
  sum.item.var <- tr(r.items)
  var.scales <- sum(r.items)
  alpha.scale <- (var.scales - sum.item.var)*nvar/((nvar-1)*var.scales)
  av.r <- alpha.scale/(nvar - alpha.scale*(nvar-1))  #alpha 1 = average r
   item.cor <- cor(items,scores,use="pairwise")    #this does not correct for item overlap
   if(is.null(ilabels)) {ilabels <-  paste("I",1:nvar,sep="")}
   rownames(item.cor) <- ilabels
   
   if (!totals) {scores <- scores/nvar }
    item.stats <- cbind(key,response.freq,item.cor,item.stats)
    colnames(item.stats)[alternatives+2] <- "r"
   
   if(short) {results <- list(item.stats=round(item.stats,digits),alpha=round(alpha.scale,digits), av.r=round(av.r,digits),Call=cl)} else 
   if (sum(miss.rep) >0) {results <-list(scores=scores,missing = miss.rep,item.stats=round(item.stats,digits),alpha=round(alpha.scale,digits), av.r=round(av.r,digits))} else{  
    results <- list(scores=scores,item.stats=item.stats,alpha=round(alpha.scale,digits), av.r=round(av.r,digits),Call=cl)}  

 class(results) <- c("psych","mchoice")
 return(results) } else {return (items)}  
 }
 
 
  #introduce a function to get cell frequencies and compensate for possible different number of response alternatives
 

response.frequencies <-
    function (items, max = 10,uniqueitems=NULL)
{
    min.item <- min(items, na.rm = TRUE)
    max.item <- max(items, na.rm = TRUE)
   if(is.null(uniqueitems))  uniqueitems <- unique(as.vector(unlist(items)))
    if ((max.item - min.item > max) ||
        (nlevels(factor(items[, 1])) > max) ||
        length(uniqueitems) > max) {
        frequency <- NULL
    }
    else {
        n.var <- dim(items)[2]
        n.cases <- dim(items)[1]
        dummy <- matrix(rep(uniqueitems, n.var), ncol = n.var)
        colnames(dummy) <- names(items)
        xdum <- rbind(items, dummy)
        frequency <- apply(xdum, 2, table)
        frequency <- t(frequency - 1)
        responses <- rowSums(frequency)
        frequency <- frequency/responses
        miss <- 1 - responses/n.cases
        frequency <- cbind(frequency, miss)
    }
    return(frequency)
}
	
 #version of Sept 19, 2007
 #revised Sept 3, 2010 to count missing responses
 #revised Sept 7, 2010 to correctly handle missing data in terms of finding alpha and correlation
 #revised April 3, 2011 to incorporate very nice suggestion by Joshua Wiley to handle unique categories
 #revised August 4, 2012 to allow the specification of unique items-- useful for irt.responses