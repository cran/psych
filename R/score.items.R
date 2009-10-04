"score.items"  <- 
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE,impute="median",  min=NULL,max=NULL,digits=2,short=FALSE) {
   cl <- match.call()
   raw.data <- TRUE
   keys <- as.matrix(keys)   #just in case they were not matrices to start with
    if ((dim(items)[1] == dim(items)[2])) {
     raw.data <- FALSE
     C <- as.matrix(items)
     cov.scales <- t(keys) %*% C %*% keys
     
    } else {
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
        items[miss]<- item.med[miss[,2]]}   #this only works if items is a matrix
        }
      
       scores<- items %*%  keys  #this actually does all the work but doesn't handle missing values
       C <- cov(items,use="pairwise")
       cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
       
       }
                  
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
    n.subjects <- dim(items)[1]
   
    
    slabels <- colnames(keys)
    if (is.null(slabels)) {
    	if (totals) {slabels<- paste("S",1:n.keys,sep="")} else {
    	             slabels <- paste("A",1:n.keys,sep="")} }
   
    
    abskeys <- abs(keys)

    item.var <- diag(C)  #find the item variances
    
    var.scales <- diag(cov.scales)
    cor.scales <- cov2cor(cov.scales)    
    sum.item.var <- item.var %*% abskeys  
   
    
    num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
   alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
    colnames(alpha.scale) <- slabels
  alpha.scale[is.nan(alpha.scale)] <- 1
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   
   #now find the Guttman 6 * reliability estimate as well as the corrected item-whole correlations
   
   if(raw.data) { item.cor <- cor(items,scores)} else {if (n.keys >1) {
         item.cor <- C %*% keys %*% diag(1/sqrt(var.scales))/sqrt(item.var)} else {item.cor <- C %*% keys /sqrt(var.scales * item.var)}}
    c.smc <- smc(C,TRUE)
    diag(C) <- c.smc
    sum.smc <- c.smc %*% abskeys
    G6 <- (var.scales - sum.item.var + sum.smc)/var.scales
    corrected.var <- diag(t(keys) %*%  C %*% keys)
    if(n.keys>1) {
    item.rc <- (C %*% keys) %*% sqrt(diag(1/corrected.var))/sqrt(item.var)} else {
      item.rc <- C %*% keys /sqrt(corrected.var*item.var) }
    colnames(item.rc) <- slabels
   
   
  
  if(is.null(ilabels)) {ilabels <- colnames(items) }
  if(is.null(ilabels)) {ilabels <-  paste("I",1:n.items,sep="")}
   rownames(item.rc) <- ilabels
  if(raw.data) {
  correction <- (colSums(abs(keys)-(keys))/2)*(max+min) #correct for flipping
   scores <- scores  + matrix(rep(correction,n.subjects),byrow=TRUE,nrow=n.subjects)
  
   
   if (!totals) {
    if(n.keys > 1) {scores <- scores %*% diag(1/num.item)   #find averages
                  }  else
                 {scores <- scores/num.item } }
            colnames(scores) <- slabels          
         } else {scores <- NULL}
    scale.cor <- correct.cor(cor.scales,t(alpha.scale))
  rownames(alpha.scale) <- "alpha"
  rownames(av.r) <- "average.r"
  rownames(G6) <- "Lambda.6"
   if (short | !raw.data) {results <- list(alpha=alpha.scale, av.r=av.r, n.items = num.item,  item.cor = item.cor,cor = cor.scales ,corrected = scale.cor,G6=G6,item.corrected = item.rc ,Call=cl)} else {
   if(raw.data) {if (sum(miss.rep) >0) {results <-list(scores=scores,missing = miss.rep,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits),cor = round(cor.scales,digits) ,corrected = round(scale.cor,digits),G6=G6,item.corrected = item.rc ,Call=cl)} else{  
    results <- list(scores=scores,alpha=round(alpha.scale,digits), av.r=round(av.r,digits), n.items = num.item,  item.cor = round(item.cor,digits), cor = round(cor.scales,digits),corrected = round(scale.cor,digits),G6=G6,item.corrected = item.rc ,Call=cl)} }
   }
   class(results) <- c("psych", "score.items")
    return(results)
 }
 #modified June 1 to add row names to items 
 #modified June 22 to add median imputation
 #modified August 8 to add colnames to scores
 #modified Sept 23, 2007 to allow for short output
 #modified December 10, 2007 to default to ilabels as colnames(items)
 #modified March, 2009 to better use the print.psych function
 #modified March 22, 2009 to add G6 and corrected for item overlap correlation
 #also allow for correlation/covariance matrix input
