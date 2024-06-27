#developed 6/4/24 
#correlation or covariance equivalent of FactoMineR::coeffRV or MatrixCorrelation::RV
#works with raw data or covariance or correlation matrices
#
RV <- function(x,y,xy=NULL, data= NULL, cor="cor",correct=0) {
    cl <- match.call()
    if(is.numeric(y )& !is.null(data)) y <- (data)[,y]
    if(is.numeric(x ) & !is.null(data)) x <- (data)[,x]
   if(is.null(xy)){ 
   switch(cor,
   cor ={xy <- cor(x,y,use="pairwise")
   x <- cor(x,use="pairwise")
   y <- cor(y,use="pairwise")
   },
    spearman ={xy <- cor(x,y,method="spearman",use="pairwise")
   x <- cor(x,method="spearman",use="pairwise")
   y <- cor(y,method="spearman",use="pairwise")
   },
   tet ={
   xy <- tetrachoric(x,y,correct=correct)$rho
   x <- tetrachoric(x,correct=correct)$rho
   y <- tetrachoric(y,correct=correct)$rho
   },
   poly = {xy <- polychoric(x,y,correct=correct)$rho
   x <- polychoric(x,correct=correct)$rho
   y <- polychoric(y,correct=correct)$rho
   },
   cov = {xy <- cov(x,y,use="pairwise")
   x <- cov(x,use="pairwise")
   y <- cov(y,use="pairwise")}
   )
  
   }
   
  m1 <- cbind(x,xy)
  m2 <- cbind(t(xy),y)
  m <- rbind(m1,m2) 
  
  Rset <- 1- det(m)/(det(x)* det(y))   
  Ru <- sum(xy)/sqrt(sum(x)*sum(y))          
 RV <- tr(xy %*% t(xy))/sqrt(tr(x %*% t(x)) * tr(y %*% t(y)))
 result <- list(RV=RV, Rset = Rset, Ru = Ru, Rx =x,Ry=y,Rxy=xy, Call=cl)
 class(result) <- c("psych","RV")
   return(result)
   }