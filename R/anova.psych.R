#A function to report the difference between two factor models
#adapted from John Fox's sem anova 
#modified November 29, 2019 to include anovas for setCor and mediate models 
#modified June,16, 2024 to include iclust models
#modified December 22, 2025 to handle CFA models
anova.psych <- function(object,...) {

Call <- match.call()
n.models <- length(Call)-1

#figure out which kind of anova to do
 if(length(class(object)) > 1)  {
    names <- cs(omega,fa, lmCor,mediate, iclust,CFA, CFA.bifactor)
    value <- inherits(object,names,which=TRUE)   # value <- class(x)[2]
    if(any(value > 1) ) { value <- names[which(value > 0)]} else {value <- "other"}
    
     } else {value <- "other"}
if(value=="CFA.bifactor") {value <- "fa"}
if(value=="CFA") {value <- "fa"}
if(value=="omega") {value <- "fa"}
        
#this does the work for lmCors and mediate or any model that returns SSR and dfs
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

 table.df<- data.frame(df=unlist(dfs),SSR=unlist(SSR))
 MSR <- table.df$SSR/table.df$df
  df <- table.df$df
  diffSS <- -diff(table.df$SSR)
  diffdf <- -diff(table.df$df) 
  
  #find the model with the most df
    biggest.df <- order(table.df$df)[1]
    scale <- table.df[biggest.df,"SSR"]/table.df[biggest.df,"df"]
   F <- (diffSS) /diffdf/scale
   
   prob <- pf(F,abs(diffdf),df[-1],lower.tail=FALSE)
   table.df<- data.frame(table.df,diff.df=c(NA,diffdf),diff.SS= c(NA,diffSS),F= c(NA,F),c(NA,prob))
   names(table.df) <- c("Res Df","Res SS", "Diff df","Diff SS","F","Pr(F > )")
   return(table.df)}
   
#and this does the work for fa and omega
another.function <- function(models,dfs,echis,chi,BICS, RMSEA, rms, null.chi,null.dof) {
   mod0 <- "Null Model"
   models <-c(mod0,models)
   mods <- unlist(models)
    n.mod <- length(models) 
    

   dfs <- c(null.dof[[1]],dfs)
    chi <- c( null.chi[[1]],chi)
 mods <- unlist(models)
  
 for(i in 1:(n.mod)) {
   temp <- unlist(mods[[i]])
   cat("Model",i, "= ")
   print(temp,rownames=FALSE)
    }
   
  delta.df <- c(NA,-diff(unlist(dfs)))
  if(!is.null(chi[[1]])){delta.chi<- c(NA,diff(unlist(chi)))} else {delta.chi <- NA
              chi[1:length(chi) ] <- NA
}

 if(!is.null(echis) ) {delta.echi <- c(NA,NA,-diff(unlist(echis)))} else {delta.echi <- NA}
  if(!is.null(BICS[[1]])){delta.bic <- c(NA,NA,diff(unlist(BICS)))} else {delta.bic <- NA
              BICS[1:length(BICS) ] <- NA
}
  #delta.bic <- diff(unlist(BICS))
  test.chi <- c(delta.chi)/delta.df
test.echi <- c(delta.echi)/delta.df
 p.delta <- pchisq(-delta.chi, delta.df, lower.tail=FALSE)
 
 if(is.null(RMSEA[[1]])) RMSEA[1:length(RMSEA)] <- NA

  if(!is.null(echis) ){

  table.df <- data.frame(df=unlist(dfs),chiSq=unlist(chi),Prchi=1-pchisq(unlist(chi), unlist(dfs)),d.df=c(delta.df), d.chiSq=c(delta.chi),
  PR=c(p.delta),
  #test=c(test.chi) , empirical = c(NA,unlist(echis)), 
  #d.empirical=c(delta.echi)
  #test.echi=c(test.echi),
  BIC=c(NA,unlist(BICS))
 # ,  d.BIC = c(delta.bic)
  ,RMSEA=c(NA,unlist(RMSEA)), rms =c(NA,unlist(rms)))
  
  
   } else {
  table.df<- data.frame(df=unlist(dfs),chiSq=unlist(chi),Prchi=1-pchisq(unlist(chi), unlist(dfs)),d.df=c(NA,delta.df), d.chiSq=c(NA,delta.chi),
  PR=c(NA,p.delta),test=c(NA,test.chi),BIC=unlist(BICS),d.BIC = c(NA,delta.bic),
  RMSEA=unlist(RMSEA), rms =unlist(rms))}

 table.df <- round(table.df,2)
 rownames(table.df)[1] <- noquote(as.character("NULL"))
 for (i in  2:n.mod ) {rownames(table.df)[i] <- noquote(as.character(Call[[i]]))}
# table.df <-  order(table.df)
return(table.df)
}
 

switch(value,

mediate ={
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)	
 	 dfs <- lapply(objects, function(x) x$cprime.reg$df)
  	SSR <- lapply(objects, function(x) x$cprime.reg$SE.resid^2 * x$cprime.reg$df)
 	 models <- lapply(objects, function(x) x$Call)
 	 table.df<- small.function(models=models,dfs=dfs,SSR=SSR)
	 }
},

lmCor ={
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)	
 	 dfs <- lapply(objects, function(x) x$df[2])
 	 
  	 SSR <- lapply(objects, function(x) x$SE.resid^2 * x$df[2])
 	 models <- lapply(objects, function(x) x$Call)
     if (length(SSR) ==1) {table.df<-  small.function(models=models,dfs=dfs,SSR=SSR)} else  {
 	 table.df<- list()
 	 nvar <- length(SSR[[1]])
 	 ssrm <- matrix(unlist(SSR),ncol=nvar,byrow=TRUE)
 	 for (i in (1:nvar)) {
 	 table.df[[names(SSR[[1]][i])]] <-  small.function(models=models,dfs=dfs,SSR=ssrm[,i],nvar=i)
 	 } 
 	 
	 }
	 }
	 return(table.df)
   },


fa = {
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)
 	 n.models <- length(objects)
 	 rms <- lapply(objects,function(x) x$stats$rms)
 	 RMSEA <- lapply(objects,function(x) x$stats$RMSEA[1])
 	 echis <- lapply(objects,function(x) x$stats$chi)	
 	 BICS <-  lapply(objects,function(x) x$stats$BIC)	
     dofs <-   lapply(objects,function(x) x$stats$dof)
     chi <-  lapply(objects,function(x) x$stats$STATISTIC)
     null.chi <- lapply(objects,function(x) x$stats$null.chisq)
     null.dof <- lapply(objects,function(x) x$stats$null.dof)
     models <- lapply(objects, function(x) x$Call)
     nechi <- length (echis)
     nBics <- length(BICS)
     nchi <- length(chi)
     
     if(nechi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
      if(nchi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
     table.df<- another.function(models,dfs=dofs,echis=echis,chi = chi,BICS = BICS,RMSEA=RMSEA,rms=rms, null.chi=null.chi,null.dof=null.dof)	
}
},

CFA = {  #not used because we treat it like fa
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)
 	 n.models <- length(objects)
 	 rms <- lapply(objects,function(x) x$stats$rms)
 	 RMSEA <- lapply(objects,function(x) x$stats$RMSEA[1])
 	echis <- lapply(objects,function(x)  x$stats$chi)	
 	 BICS <-  lapply(objects,function(x)  x$stats$BIC)	
     dofs <-   lapply(objects,function(x)  x$stats$dof)
     chi <-  lapply(objects,function(x)  x$stats$STATISTIC)
     models <- lapply(objects, function(x) x$Call)
     nechi <- length (echis)
     nBics <- length(BICS)
     nchi <- length(chi)
     
     if(nechi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
      if(nchi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
     table.df<- another.function(models,dfs=dofs,echis=echis,chi = chi,BICS = BICS,RMSEA=RMSEA,rms=rms)	
}
},


iclust = {
 if (length(list(object, ...)) > 1L)  {
 	 objects <- list(object,...)
 	 n.models <- length(objects)
 	  rms <- lapply(objects,function(x) x$stats$rms)
 	 RMSEA <- lapply(objects,function(x) x$stats$RMSEA[1])
 	 echis <- lapply(objects,function(x) x$stats$chi)	
 	 BICS <-  lapply(objects,function(x) x$stats$BIC)	
     dofs <-   lapply(objects,function(x) x$stats$dof)
     chi <-  lapply(objects,function(x) x$stats$STATISTIC)
     models <- lapply(objects, function(x) x$call)
     nechi <- length (echis)
     nBics <- length(BICS)
     nchi <- length(chi)
     
     if(nechi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
      if(nchi != n.models)  {stop("You do not seem to have chi square values for one of the models ")}
     table.df<- another.function(models,dfs=dofs,echis=echis,chi = chi,BICS = BICS, RMSEA=RMSEA,rms=rms)	
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
   table.df <- another.function(models,dfs=dofs,echis=NULL,chi = chi,BICS = BICS)	
   
   }
}
)     #end of switch

return(table.df)
			}
