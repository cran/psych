
#A number of estimates of unidimensinality
#Developed March 9. 2017

"unidim" <- function(x,keys.list =NULL,flip=FALSE) {
   cl <- match.call()
   
   n.keys <- 1
    all.x <- x
   results <- list()
   if(!is.null(keys.list)) {n.keys <- length(keys.list)
    } else {keys.list <- NULL }
   
   for(keys in 1:n.keys) { if(!is.null(keys.list)) {
     
   select <- keys.list[[keys]]
        flipper <- rep(1,length(select))
         flipper[grep("-",select)] <- -1
         if(is.numeric(select)) {select <- abs(select) } else {
         select <- sub("-","",unlist(select)) }
   
   x <- all.x[,select]}  else {flipper <- rep(1,ncol(x))} #this allows us to handle multiple scales 
   
 if(!isCorrelation(x) ) x <- cor(x,use="pairwise")
  f1 <- fa(x)
  g <- sum(f1$model)  # sum(f1$loadings %*% t(f1$loadings))
  n <- nrow(x)
  Vt <- sum(x)
   om.g <- g/Vt                          #model/ r
  om.t <- (Vt - sum(f1$uniqueness))/Vt   #total reliability 

 uni.orig <- g/ (Vt - sum(f1$uniqueness))  #raw unidimensionality
 
   
 #now, to find traditional alpha, we need to flip negative items
 if(flip | n.keys == 1) { flipper <- rep(1,n)
 flipper[sign(f1$loadings ) < 0] <- -1 }
  x <- diag(flipper) %*% x %*% diag(flipper)
  Vt <- sum(x)
 
 alpha.std <-  (1- n/Vt)*(n/(n-1))
 av.r <- (Vt-n)/(n*(n-1))
  omega.flip <- sum(diag(flipper) %*% f1$model %*% diag(flipper))/Vt
  omega.total.flip <-  (Vt - sum(f1$uniqueness))/Vt
  flipperped.loadings <- flipper * f1$loadings
  g.flipperped <- sum(flipperped.loadings%*% t(flipperped.loadings))
  uni.flipper <- g.flipperped/(Vt - sum(f1$uniqueness))

  stats <- list(uni=uni.orig,uni.flipper = uni.flipper,alpha=alpha.std,av.r = av.r,om.g=om.g, omega.pos = omega.flip,om.t=om.t,om.total.flip= omega.total.flip)
  if(!is.null(keys.list)) {results[[names(keys.list[keys])]]<- stats } else {results <- stats}
  }
  temp <- matrix(unlist(results),ncol=8,byrow=TRUE)
  colnames(temp) <- c("Raw Unidim","Adjusted","alpha","av.r","original model","adjusted model", "raw.total", "adjusted total")
  rownames(temp) <- names(keys.list)
  results <- list(uni=temp)
  results$Call <- cl
  class(results) <- c("psych","unidim")
  
  return(results)
  }
  
  print.psych.unidim <- function(x,digits=2) {
  cat("\nA measure of unidimensionality \n Call: ")
  print(x$Call)
  
  cat("\nUnidimensionality index = \n" )
  print(round(x$uni,digits=digits))
  
 cat("\nunidim adjusted index reverses negatively scored items.")
  cat("\nalpha ","  Based upon reverse scoring some items.")
  cat ("\naverage correlations are based upon reversed scored items") 
     }
  
  
  
  