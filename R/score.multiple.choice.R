 "score.multiple.choice" <-  
  function(key,data,score=TRUE,totals=FALSE,ilabels=NULL, missing=TRUE,impute="median", digits=2,short=TRUE) {
  #convert a data matrix or data with multiple choice responses to correct/incorrect
  if(!is.matrix(data)) {if(!is.data.frame(data)) {stop("data must be either a data frame or matrix!")} else data <- as.matrix(data)}
  nvar <- dim(data)[2]
  
  #introduce a function to get cell frequencies and compensate for possible different number of response alternatives
  
  cell.frequencies <- function(x) {
	min.item <- min(x,na.rm=TRUE)
	max.item <- max(x,na.rm=TRUE)
	n.var <- dim(x)[2]
	n.cases <- dim(x)[1]
	dummy <- matrix(rep(min.item:max.item,n.var),ncol=n.var)
	colnames(dummy) <- names(x)
	xdum <- rbind(x,dummy)
	frequency <- apply(xdum,2,table)
	frequency <- (frequency -1)/n.cases
	return(t(frequency))
	}
  
  response.freq <- cell.frequencies(data)
  
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
                       #first, define a local function
 

  scores<- items %*%  keys  #this actually does all the work but doesn't handle missing values
    
    slabels <- colnames(keys)
    if (is.null(slabels)) {
    	if (totals) {slabels<- paste("Totals") } else {
    	             slabels <- paste("Averages")} }
   colnames(scores) <- slabels
  
   
    item.var <- diag(var(items,use="pairwise"))  #find the item variances
    var.scales <- var(scores)
    sum.item.var <- sum(item.var)
     alpha.scale <- (var.scales - sum.item.var)*nvar/((nvar-1)*var.scales)
     av.r <- alpha.scale/(nvar - alpha.scale*(nvar-1))  #alpha 1 = average r
       item.cor <- cor(items,scores,use="pairwise")    #this does not correct for item overlap
   if(is.null(ilabels)) {ilabels <-  paste("I",1:nvar,sep="")}
   rownames(item.cor) <- ilabels
   
   
    if (!totals) {scores <- scores/nvar }
    item.stats <- cbind(key,response.freq,item.cor,item.stats)
    colnames(item.stats)[alternatives+2] <- "r"
   
   if(short) {results <- list(item.stats=round(item.stats,digits),alpha=round(alpha.scale,digits), av.r=round(av.r,digits))} else 
   if (sum(miss.rep) >0) {results <-list(scores=scores,missing = miss.rep,item.stats=round(item.stats,digits),alpha=round(alpha.scale,digits), av.r=round(av.r,digits))} else{  
    results <- list(scores=scores,item.stats=item.stats,alpha=round(alpha.scale,digits), av.r=round(av.r,digits))}  
    


 return(results) } else {return (items)}  
 }
 #version of Sept 19, 2007
