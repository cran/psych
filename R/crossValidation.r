"crossValidation" <- function(model,data, options=NULL,select=NULL) {
 cl <- match.call()
obnames <- cs(setCor,bestScales )
     value <- inherits(model, obnames, which=TRUE)
			   if (any(value > 1)) {value <- obnames[which(value >0)]} else {value <- "none"}


switch(value, 
setCor = {wt <- model$coefficients
    if(rownames(wt)[1] =="(Intercept)") wt <- wt[-1,,drop=FALSE] },
    
bestScales ={#the default
if(is.null(options)){options <- "best.keys"}
   switch(options,
 best.keys = {  keys <- model$best.keys
  	select <- c(names(keys),selectFromKeys(keys))
  	data <- data[,select]
  	nvar <- ncol(data)
	 wt <- make.keys(data,keys)
    keys <- model$best.keys },
 weights =  {wt <- model$weights},     
 optimal.keys ={ 
      keys <- model$optimal.keys
  		select <- c(names(keys),selectFromKeys(keys))
  	data <- data[,select]
  	nvar <- ncol(data)
	 wt <- make.keys(data,keys)
     },
 optimal.weights ={ wt <- model$optimal.weights}
       )
       

   }, 
    
 none = {wt<- model} #from setCor
    )
    
   nvar <- nrow(wt)
   predictors <- rownames(wt)
   criteria <- colnames(wt)
   vars <- c(criteria,predictors)
    npred <-  length(predictors)
    ncrit <-length(criteria)
  if(isCorrelation(data)) {raw <- FALSE
    R <-data[vars,vars]
  }  else {raw<- TRUE
     data <- data[vars ]  #select just the variables in the model
  
    R <- cor(data,use="pairwise")
    }  #this is a matrix of the (criteria + predictors)
    cv <- rep(NA,ncrit)
    Mij <- matrix(NA,ncrit,ncrit)
    for(i in 1:ncrit)  {y <- wt[,i,drop=FALSE]
     wtm <- diag(wt[,i])
     Cxx <- wtm %*% R[predictors,predictors] %*% wtm
     Vxx <- sum(Cxx)
     
     
     
     Cxyi  <- wt[,i,drop=FALSE] * R[i,predictors]
     Cxxy <- wt[,i] *t(R[1:ncrit,predictors,drop=FALSE])
     cv[i]<-  sum(Cxyi)/sqrt(Vxx)

     #Cxyi <- (Cxyi)/sqrt(Vxx)
     Mij[,i] <-  colSums( Cxxy)/sqrt(Vxx)
   
     } 
    Mij <- t(Mij)
   

     names(cv) <- criteria

     #now do the raw correlation approach
    if(raw) {x.data <- scale(data[predictors])
   
     pred <- matrixMult.na(x.data,wt,scale = FALSE) 
     item.R <- cor(pred,data[criteria],use="pairwise")
      colnames(Mij) <- rownames(Mij)<-colnames(item.R)
     cv.items <-diag(item.R)
     crossValid <- data.frame(items = cv.items, mat = cv)} else {CrossValid <- cv}
     nvars <- length(vars)
     result <- list(crossV = crossValid,nvars=npred,item.R = item.R,mat.R = Mij,Call =cl)
     class(result) <- c("psych","crossV")
     return(result)
    
         }

#small function to call matplot with labels         
matPlot <- function(x, type = "b", minlength=6, xlas=0,legend=NULL,lab=NULL,pch=16,col=1:6,lty=NULL,...) {
     rownames(x)<- abbreviate(rownames(x),minlength=minlength)
      nvar <- NCOL(x)
       if(missing(pch)) pch <- seq(15,(15+nvar))
    matplot(x,type=type,xaxt="n",pch=pch,col=col,...)
    #now add the labels 
     at1 <- (1:nrow(x))
     labx <- (rownames(x))
     axis(1,at=at1,labels=labx,las=xlas)
     if(!is.null(legend)) {
    
     lab <- colnames(x)
    
    if(missing(lty)) lty <- 1:8
      legend.location <- c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center","none")
      legend(legend.location[legend],legend=lab,lty=lty,pch = pch,col=col)}
     }
     
 #Fixed 9/9/20 to handle case of single criterion (by using drop=FALSE in appropriate places)
