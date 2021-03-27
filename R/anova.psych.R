#A function to report the difference between two factor models
#adapted from John Fox's sem anova 
#modified November 29, 2019 to include anovas for setCor and mediate models 

anova.psych <- function(object,...) {
#if(length(class(object)) > 1)  { value <- class(object)[2] } else {value <- NA}

 if(length(class(object)) > 1)  {
    names <- cs(omega,fa, setCor,mediate)
    value <- inherits(object,names,which=TRUE)   # value <- class(x)[2]
    if(any(value > 1) ) { value <- names[which(value > 0)]} else {value <- "other"}
    
     } else {value <- "other"}

#this does the work for setCor and mediate or any model that returns SSR and dfs
small.function <- function(models,dfs,SSR,nvar=1) {
  #this next section is adapted  from anova.lm and anova.lmlist

if(nvar==1) { n.mod <- length(models)
 mods <- unlist(models)
 for(i in 1:n.mod) {
   temp <- unlist(mods[[i]])
   cat("Model",i, "= ")
   print(temp,rownames=FALSE)
    }
    }

 table <- data.frame(df=unlist(dfs),SSR=unlist(SSR))
 MSR <- table$SSR/table$df
  df <- table$df
  diffSS <- -diff(table$SSR)
  diffdf <- -diff(table$df) 
  
  #find the model with the most df
    biggest.df <- order(table$df)[1]
    scale <- table[biggest.df,"SSR"]/table[biggest.df,"df"]
   F <- (diffSS) /diffdf/scale
   
   prob <- pf(F,abs(diffdf),df[-1],lower.tail=FALSE)
   table <- data.frame(table,diff.df=c(NA,diffdf),diff.SS= c(NA,diffSS),F= c(NA,F),c(NA,prob))
   names(table) <- c("Res Df","Res SS", "Diff df","Diff SS","F","Pr(F > )")
   return(table)}
   
#and this does the work for fa and omega
another.function <- function(models,dfs,echis,chi,BICS) {
   mods <- unlist(models)
    n.mod <- length(models)
 mods <- unlist(models)
 for(i in 1:n.mod) {
   temp <- unlist(mods[[i]])
   cat("Model",i, "= ")
   print(temp,rownames=FALSE)
    }
  delta.df <- -diff(unlist(dfs))
  delta.chi <- -diff(unlist(chi))
 if(!is.null(echis) ) {delta.echi <- -diff(unlist(echis))} else {delta.echi <- NA}
  delta.bic <- diff(unlist(BICS))
  test.chi <- delta.chi/delta.df
test.echi <- delta.echi/delta.df
 p.delta <- pchisq(delta.chi, delta.df, lower.tail=FALSE)
 if(!is.null(echis) ){
  table <- data.frame(df=unlist(dfs),d.df=c(NA,delta.df),chiSq=unlist(chi), d.chiSq=c(NA,delta.chi),
  PR=c(NA,p.delta),test=c(NA,test.chi), empirical = unlist(echis),d.empirical=c(NA,delta.echi),test.echi=c(NA,test.echi),BIC=unlist(BICS),d.BIC = c(NA,delta.bic))} else {
  table <- data.frame(df=unlist(dfs),d.df=c(NA,delta.df),chiSq=unlist(chi), d.chiSq=c(NA,delta.chi),
  PR=c(NA,p.delta),test=c(NA,test.chi),BIC=unlist(BICS),d.BIC = c(NA,delta.bic))}

 table <- round(table,2)
return(table)


}
 

switch(value,

mediate ={
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)	
 	 dfs <- lapply(objects, function(x) x$cprime.reg$df)
  	SSR <- lapply(objects, function(x) x$cprime.reg$SE.resid^2 * x$cprime.reg$df)
 	 models <- lapply(objects, function(x) x$Call)
 	 table <- small.function(models=models,dfs=dfs,SSR=SSR)
	 }
},

setCor ={
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)	
 	 dfs <- lapply(objects, function(x) x$df[2])
 	 
  	 SSR <- lapply(objects, function(x) x$SE.resid^2 * x$df[2])
 	 models <- lapply(objects, function(x) x$Call)
     if (length(SSR) ==1) {table <-  small.function(models=models,dfs=dfs,SSR=SSR)} else  {
 	 table <- list()
 	 nvar <- length(SSR[[1]])
 	 ssrm <- matrix(unlist(SSR),ncol=nvar,byrow=TRUE)
 	 for (i in (1:nvar)) {
 	 table[[names(SSR[[1]][i])]] <-  small.function(models=models,dfs=dfs,SSR=ssrm[,i],nvar=i)
 	 } 
 	 
	 }
	 }
	 return(table)
   },


fa = {
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)
 	 n.models <- length(objects)
 	 
 	 echis <- lapply(objects,function(x) x$chi)	
 	 BICS <-  lapply(objects,function(x) x$BIC)	
     dofs <-   lapply(objects,function(x) x$dof)
     chi <-  lapply(objects,function(x) x$STATISTIC)
     models <- lapply(objects, function(x) x$Call)
     nechi <- length (echis)
     nBics <- length(BICS)
     nchi <- length(chi)
     
     if(nechi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
      if(nchi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
     table <- another.function(models,dfs=dofs,echis=echis,chi = chi,BICS = BICS)	
}
},

omega = {   #should change this to include more than 2 models (see above )
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)
 	 n.models <- length(objects)

 	# echis <- lapply(objects,function(x) x$schmid$chi)	
 	 BICS <-  lapply(objects,function(x) x$schmid$BIC)	
     dofs <-   lapply(objects,function(x) x$schmid$dof)
     chi <-  lapply(objects,function(x) x$schmid$STATISTIC)
     models <- lapply(objects, function(x) x$Call)
    # nechi <- length (echis)
     nBics <- length(BICS)
     nchi <- length(chi)
   table <- another.function(models,dfs=dofs,echis=NULL,chi = chi,BICS = BICS)	
   
   }
}
)
return(table)
			}
