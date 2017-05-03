#Developed February, 2017 
#closely follows chapters by Pat Shrout and Sean Lane in terms of the statistics
"multilevel.reliability" <- "mlr" <- 
 function(x,grp="id",Time="time",items=c(3:5),alpha=TRUE,icc=FALSE,aov=TRUE,lmer=FALSE,
  lme = TRUE,long=FALSE,values=NA,na.action="na.omit",plot=FALSE,
  main="Lattice Plot by subjects over time")  {
 cl <- match.call() 
 
s.lmer <- lmer.MS <- MS_id <- s.aov <-NULL 
MS.df <- data.frame( matrix(NA,nrow=8,ncol=2))  #a filler in case we don't have it
colnames(MS.df) <- c("Variance","Percent")

 if(!long) {#the normal case is wide data which we analyze and then convert to long
 
#first check if we should reverse any items and convert location numbers (if specified) to location names
 n.items <- length(items)

  if(is.character(items)) {
    temp <- rep(1,n.items)
  
   temp [strtrim(items,1)=="-"] <- -1
   if(any(temp < 0) )  {items <- sub("-","",items) }
   } else {temp <- sign(items)
      items <- colnames(x)[abs(items)] 
    }
  if(any(temp < 0)) {
   
   min.item <- min(x[items],na.rm=TRUE)
   max.item <- max(x[items],na.rm=TRUE)
   x[items[temp <0]] <- max.item- x[items[temp <0]] + min.item 
   }
   
wide <- x[items]
rwname <- unique(x[grp])
n.obs <- nrow(rwname)
n.time <- length(table(x[Time]))
n.items <- ncol(wide)

if(alpha) {
alpha.by.person <- by(x,x[grp],function(x) alphaBy(x[items]))
rnames <- paste0("ID",names(alpha.by.person))
alpha.by.person <- matrix(unlist(alpha.by.person),ncol=4,byrow=TRUE)
colnames(alpha.by.person) <- c("Raw alpha","Std. alpha","av.r","signal/noise")
rownames(alpha.by.person) <- rnames
} else {alpha.by.person <- NA}

if(icc) {   #this takes a long time, and can be skipped
icc.by.person <- by(x,x[grp],function(x) ICC(x[items]))
icc.by.time <- by(x,x[Time],function(x) ICC(x[items]))

icc.person.sum <- matrix(NA,nrow=n.obs,ncol=6)  #add more columns when we think about what to put there
icc.time.sum <- matrix(NA,nrow=n.time,ncol=6)

 for(person in 1:n.obs) {icc.person.sum[person,] <- c(icc.by.person[[person]][["results"]][["ICC"]][c(3,6)],
 icc.by.person[[person]][["results"]][["lower bound"]][3],icc.by.person[[person]][["results"]][["upper bound"]][3],
 icc.by.person[[person]][["results"]][["lower bound"]][6],icc.by.person[[person]][["results"]][["upper bound"]][6])
   }
 rownames(icc.person.sum) <- paste0("ID",names(icc.by.person))
 colnames(icc.person.sum) <- c("ICC single" ,"ICC summed","Lower Bound ICC13","Upper Bound ICC13", "Lower Bound ICC23","Upper Bound ICC23")
 for(time in 1:n.time) {icc.time.sum[time,] <- c(icc.by.time[[time]][["results"]][["ICC"]][c(3,6)],
           icc.by.time[[time]][["results"]][["lower bound"]][3],icc.by.time[[time]][["results"]][["upper bound"]][3],
 icc.by.time[[time]][["results"]][["lower bound"]][6],icc.by.time[[time]][["results"]][["upper bound"]][6])         
  }  
  rownames(icc.time.sum) <- paste0("time",1:n.time)
   colnames(icc.time.sum) <- c("ICC single" ,"ICC summed","Lower Bound ICC13","Upper Bound ICC13", "Lower Bound ICC23","Upper Bound ICC23")


} else {icc.by.person <- icc.by.time <- icc.time.sum <- icc.person.sum <-  NA}
#this next part treats the possibility of missing times 

 long <- NULL
  long.list <- by(x,x[grp],function(xx) {xx 
     y.df <- data.frame(id = as.factor(xx[,grp]), time=as.factor(xx[,Time]),stack(xx[items]))
          }) 
 for (i in 1:n.obs) {long <- rbind(long,long.list[[i]]) }
colnames(long)[4] <- "items"
} else {#we have long data already, but need to add a few values to make it work
   long <- x
   n.items <- length(table(long[items]))  #this does not 
   n.obs <- length(table(long[grp]))
   n.time <- length(table(long[Time]))
   long <- long[c(grp,Time,items,values)]
   colnames(long) <- c("id","time","items","values")
   alpha.by.person <- icc.person.sum <- icc.time.sum <-  icc.by.person <- icc.by.time <- NULL  }
if(lmer) { 
if (!requireNamespace('lme4')) {stop("I am sorry, to do a NREML  requires the lme4 package to be installed")} 
  mod.lmer <- lme4::lmer(values ~ 1 + (1 | id) + (1 | time) + (1 | items) + (1 | id:time)+ (1 | id:items)+ (1 | items :time),
  data=long,na.action=na.action) #might want to add control option, but probably not needed
  vc <- lme4::VarCorr(mod.lmer)
  MS_id <- vc$id[1,1]
  MS_time <- vc$time[1,1]
  MS_items <- vc$items[1,1]
  MS_pxt <- vc[[1]][[1]]
  MS_pxitem <- vc[[2]][[1]]
  MS_txitem <- vc[[3]][[1]]
  error <- MS_resid <- (attributes(vc)$sc)^2
  s.lmer <- s.aov <- summary(mod.lmer)
 MS.df <- data.frame(variance= c(MS_id, MS_time ,MS_items, MS_pxt, MS_pxitem, MS_txitem, MS_resid,NA))

 rownames(MS.df) <- c("ID","Time","Items","ID x time", "ID x items", "time x items", "Residual","Total") 

 MS.df["Total",]  <- sum(MS.df[1:7,1],na.rm=TRUE)
MS.df["Percent"] <- MS.df/MS.df["Total",1]
lmer.MS <- MS.df  #save these
 } 
if(aov) {
aov.x <- aov(values ~  id + time + items + time * id + time * items + items * id , data = long)
s.aov <- summary(aov.x)
stats <- matrix(unlist(s.aov),ncol=5, byrow=FALSE)
colnames(stats) <- c("df","SS", "MS", "F", "p")
rownames(stats) <- c("id","time","items", "id x time","time x items","id x items","residuals")

MS_id <- (stats["id","MS"] - stats["id x time","MS"] - stats["id x items","MS"] + stats["residuals","MS"]) /( n.time * n.items)
MS_time <- (stats["time","MS"] -  stats["id x time","MS"] - stats["time x items","MS"] + stats["residuals","MS"])/(n.obs*n.items)
MS_items <- (stats["items","MS"] -  stats["id x items","MS"] - stats["time x items","MS"] + stats["residuals","MS"])/(n.obs*n.time)
MS_pxt <- (stats["id x time", "MS"] - stats["residuals","MS"])/( n.items)
MS_pxitem <- (stats["id x items", "MS"] - stats["residuals","MS"])/( n.time)
MS_txitem <- (stats["time x items", "MS"] - stats["residuals","MS"])/(n.obs)

MS.df <- data.frame(variance= c(MS_id, MS_time ,MS_items, MS_pxt, MS_pxitem, MS_txitem, stats["residuals","MS"],NA))
error <- stats["residuals","MS"] 

if(any(MS.df[1:7,1] < 0) & is.null(lmer.MS)){ warning("Some of the variance estimates from ANOVA are negative.  This is probably due to missing values and an unbalanced design.  You should consider using the lme option")
   } else {if(!is.null(lmer.MS)) {
   MS_id <- lmer.MS[1,1]
   MS_time <- lmer.MS[2,1]
   MS_items <- lmer.MS[3,1]
   MS_pxt <- lmer.MS[4,1]
  MS_pxitem <- lmer.MS[5,1]
   MS_txitem <- lmer.MS[6,1]
   error <- MS_resid <- lmer.MS[7,1]
   MS.df[1] <- lmer.MS[1]
 
    }}
    


rownames(MS.df) <- c("ID","Time","Items","ID x time", "ID x items", "time x items", "Residual","Total")
MS.df["Total",]  <- sum(MS.df[1:7,1],na.rm=TRUE)
MS.df["Percent"] <- MS.df/MS.df["Total",1]
}  

if(!is.null(MS_id)) {
#now find the reliabilities  -- note the typo in equation 7 in Lane and Shrout
Rkf <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + error/(n.time * n.items)))
R1r <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + MS_time + MS_pxt + error/( n.items)))  #per Sean Lane
Rkr <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + MS_time/n.time + MS_pxt/n.time + error/( n.time * n.items)))
Rc <- (MS_pxt)/(MS_pxt + error/n.items)
} else { Rkf <- R1r <- Rkr <-Rc <- NULL}
if(lme | lmer) {

#we find these using a nested structure from lme or lmer (if available)

if(!lme & lmer) {#use lmer to do the nested test
   mod.lmer <- lme4::lmer(values ~ 1 + (1 | id/time),data=long,na.action=na.action)
   s.lme <- summary(mod.lmer)
   vc <-lme4::VarCorr(mod.lmer)
   vid <- vc$id[1,1]
   vtime_id <- vc$time[1,1]
   vres <- (attributes(vc)$sc)^2
    } else {#use lme to do the nested test
mod.lme <- nlme::lme(values ~ 1 , random = list(id =~ 1 ,time =~ 1 | id:items), data=long,na.action=na.action)
s.lme <- summary(mod.lme)
vc <- suppressWarnings(matrix(as.numeric(nlme::VarCorr(mod.lme)),ncol=2))
vid <- vc[2,1]
vtime_id <- vc[4,1]
vres <- vc[5,1]}
Rkrn <- vid/(vid + vtime_id/(n.time) +  vres/(n.time * n.items))  #this are the nested terms
Rcn <-  vtime_id/ (vtime_id + vres/n.items) 
MS.df ["id",1] <- vid 
MS.df ["id(time)",1] <- vtime_id
MS.df["residual",1] <- vres
MS.df["total",1] <- vid + vtime_id + vres
MS.df ["id",2] <- vid/ MS.df["total",1]
MS.df ["id(time)",2] <- vtime_id/MS.df["total",1]
MS.df["residual",2] <- vres/MS.df["total",1]
MS.df["total",2] <- MS.df["total",1]/MS.df["total",1]

} else {MS.df ["id",1] <- 
MS.df ["id(time)",1] <- 
MS.df["residual",1] <-
MS.df["total",1] <- 
MS.df ["id",2] <- 
MS.df ["id(time)",2] <- 
MS.df["residual",2] <- NA
MS.df["total",2] <- MS.df["total",1]/MS.df["total",1]

  s.lme <- Rkrn <- Rcn <- NULL}
# 
if(aov || lmer ||lme) {result <- list(n.obs = n.obs, n.time = n.time, n.items=n.items, components = MS.df,RkF =Rkf,R1R = R1r,RkR = Rkr,Rc=Rc,RkRn=Rkrn,Rcn = Rcn, ANOVA=s.aov,s.lmer =s.lmer,s.lme= s.lme,alpha=alpha.by.person, summary.by.person = icc.person.sum,summary.by.time=icc.time.sum, ICC.by.person = icc.by.person,ICC.by.time=icc.by.time,lmer=lmer,long = long,Call=cl)
}  else  {result <- list(n.obs = n.obs, n.time = n.time, n.items=n.items,alpha=alpha.by.person,lmer=lmer,long=long,Call=cl) } 
  
if(plot) { plot1<- xyplot(values ~ time | id, group=items, data=long, type = "b",as.table=TRUE,strip=strip.custom(strip.names=TRUE,strip.levels=TRUE),col=c("blue","red","black","grey"))
print(plot1)}
class(result) <- c("psych","multilevel")
return(result)
}

"print.psych.multilevel" <- function(x,digits=2,all=FALSE,short=TRUE) {
cat("\nMultilevel Generalizability analysis ",x$title," \n")
	cat("Call: ")
	print(x$Call)
	
	cat("\nThe data had ",x$n.obs, " observations taken over ", x$n.time ," time intervals for ", x$n.items, "items.\n")
	
	mat <- list(n.obs = x$n.obs,n.time = x$n.time,n.items = x$n.items)    #save these 
	cat("\n Alternative estimates of reliabilty based upon Generalizability theory\n")
	
	if(!is.null(x$RkF)){ cat("\nRkF  = ",round(x$RkF,digits) , "Reliability of average of all ratings across all items and  times (Fixed time effects)")
	 mat["RkF"] <- x$RkF}
	if(!is.null(x$R1R)) {cat("\nR1R  = ",round(x$R1R,digits),"Generalizability of a single time point across all items (Random time effects)")
		 mat["R1R"] <- x$R1R}
	
	if(!is.null(x$RkR)) {cat("\nRkR  = ",round(x$RkR,digits),"Generalizability of average time points across all items (Random time effects)")
	  mat["RkR"] <- x$RkR}
	if(!is.null(x$Rc))  {cat("\nRc   = ",round(x$Rc,digits),"Generalizability of change (fixed time points, fixed items) ")
	   mat["Rc"] <- x$Rc}
	if(!is.null(x$RkRn) ) {cat("\nRkRn = ",round(x$RkRn,digits),"Generalizability of between person differences averaged over time (time nested within people)")
	    mat["RkRn"] <- x$RkRn}
	if(!is.null(x$Rcn)) {cat("\nRcn  = ",round(x$Rcn,digits),"Generalizability of within person variations averaged over items  (time nested within people)")
	   mat["Rcn"] <- x$Rcn}
	
	if(!x$lmer && !is.null(x$RkF) ) {cat("\n\n These reliabilities are derived from the components of variance estimated by ANOVA \n") 
		if(!is.null(x$components)) {
	     if(!any(is.na(x$components[1:8,1])) & any(x$components[1:8,1] < 0))  { warning("The ANOVA based estimates are suspect, probably due to missing data, try using lmer")}
	        }} else {
	 if(x$lmer ) { cat("\n\n These reliabilities are derived from the components of variance estimated by lmer \n")}}
	if(!is.null(x$components) && !is.na(x$components[1,1]))  {      print(round(x$components[1:8,],digits=digits))
	mat["components"] <- list( x$components[1:8,1])}
		if(!is.null(x$components) && !is.na(x$components[9,1] )) {
		if(!x$lmer) {cat("\n The nested components of variance estimated from lme are:\n")} else {cat("\n The nested components of variance estimated from lmer are:\n")}
	   print(x$components[9:12,],digits=digits)
	   mat["lmer"] = list(x$components[9:12,1])} else {cat("\nNested components were not found because lme was not used\n")}
	
	
	if(!short) {cat("\n\n Three way ANOVA or lmer analysis \n")
	print(x$ANOVA,digits=digits)
	 cat("\nvariance components from lme(r)\n")
	 print(x$s.lme,digits=digits)
	 cat("\n Alpha reliability by subjects)\n")
	 print(x$alpha,digits) }
	
	 if(all) {
	          cat("\n Intraclass Correlations by subjects (over time and items) \n")
	          print(x$summary.by.person,digits)
	          cat("\n Intraclass Correlations by time (over subjects and items) \n")
	          print(x$summary.by.time,digits) }
	    
	         if(short) { cat("\nTo see the ANOVA and alpha by subject, use the short = FALSE option.")}
	         if(!all) {cat("\n To see the summaries of the ICCs by subject and time, use all=TRUE")}
	         cat("\n To see specific objects select from the following list:\n",names(x)[-c(1:10)])
	invisible(mat)
	         
	}
	
	
	"alphaBy" <- function(x) {
    n <- dim(x)[2]
    C <- cov(x,use="pairwise")
    R <- cov2cor(C)
    alpha.raw <- (1- tr(C)/sum(C))*(n/(n-1))
    sumR <- sum(R)
    alpha.std <-  (1- n/sum(R))*(n/(n-1))
   # smc.R <- smc(R)
   # G6 <- (1- (n-sum(smc.R))/sumR)
    av.r <- (sumR-n)/(n*(n-1))
   # mod1 <- matrix(av.r,n,n)
   # Res1 <- R - mod1
   # GF1 =  1- sum(Res1^2)/sum(R^2)
   # Rd <- R - diag(R)
   # diag(Res1) <- 0
   # GF1.off <- 1 - sum(Res1^2)/sum(Rd^2)  
    sn <- n*av.r/(1-av.r)
   # Q = (2 * n^2/((n-1)^2*(sum(C)^3))) * (sum(C) * (tr(C^2) + (tr(C))^2) - 2*(tr(C) * sum(C^2))) #corrected 1/15/16 
   # Q = (2 * n^2/((n - 1)^2 * (sum(C)^3))) * (sum(C) * (tr(C%*%C) +  (tr(C))^2) - 2 * (tr(C) * sum(C%*%C)))   #correction from Tamaki Hattori
    result <- list(raw=alpha.raw,std=alpha.std,av.r=av.r,sn =sn)
    return(result)
    }

"mlArrange" <- 
function(x,grp="id",Time="time",items=c(3:5),extra=NULL)  {
 n.items <- length(items)

  if(is.character(items)) {
    temp <- rep(1,n.items)
  
   temp [strtrim(items,1)=="-"] <- -1
   if(any(temp < 0) )  {items <- sub("-","",items) }
   } else {temp <- sign(items)
      items <- colnames(x)[abs(items)] 
    }
  if(any(temp < 0)) {
   
   min.item <- min(x[items],na.rm=TRUE)
   max.item <- max(x[items],na.rm=TRUE)
   x[items[temp <0]] <- max.item- x[items[temp <0]] + min.item 
   }
   
wide <- x[items]
rwname <- unique(x[grp])
n.obs <- nrow(rwname)
n.time <- nrow(unique(x[Time]))
n.items <- ncol(wide)

 long <- NULL
 if(is.null(extra)) {
  long.list <- by(x,x[grp],function(xx) {xx 
     y.df <- data.frame(id = (xx[,grp]), time=(xx[,Time]),stack(xx[items]))
          }) 
 for (i in 1:n.obs) {long <- rbind(long,long.list[[i]]) }
colnames(long)[4] <- "items" } else {
 long.list <- by(x,x[grp],function(xx) {xx 
     y.df <- cbind((xx[,grp]), (xx[,Time]),stack(xx[items]), (xx[,extra]),row.names=NULL)
          }) 
 for (i in 1:n.obs) {long <- rbind(long,long.list[[i]]) }
 colnames(long)[1:4] <- c("id","time","values", "items")
 colnames(long)[5:(4+length(extra))] <- colnames(x)[extra]
    }
results <- long 
return(long)
}

"mlPlot" <- function(x,grp="id",Time="time",items=c(3:5),extra=NULL,col=c("blue","red","black","grey"), main="Lattice Plot by subjects over time",...) {
long <- mlArrange(x =x,grp=grp,Time=Time,items=items,extra=extra)
 plot1<- xyplot(values ~ time | id, group=items, data=long, type = "b",as.table=TRUE,strip=strip.custom(strip.names=TRUE,strip.levels=TRUE),col=col,main=main,...)
print(plot1)
invisible(long)
}

"print.psych.multilevel.mat" <- function(x,digits=2,all=FALSE,short=TRUE) {
cat("\nMultilevel Generalizability analysis ",x$title," \n")
if(length(x) < 21 ) items <- length(x)

temp <- matrix(unlist(x),ncol=items,byrow=FALSE)

rownames(temp) <- c("n.obs","n.time","n,items","RkF","R1R","RkR","Rc","RkRn","Rcn","ID","Time","Items","ID x Time","ID x items","Time x items","residual","Total","Id","ID (time)","residual","total")

	
	
	cat("\nThe data had ",x$n.obs, " observations taken over ", x$n.time ," time intervals for ", x$n.items, "items.\n")
	
	mat <- list(n.obs = x$n.obs,n.time = x$n.time,n.items = x$n.items)    #save these 
	cat("\n Alternative estimates of reliabilty based upon Generalizability theory\n")
	
	if(!is.null(x$RkF)){ cat("\nRkF  = ",round(x$RkF,digits) , "Reliability of average of all ratings across all items and  times (Fixed time effects)")
	 mat["RkF"] <- x$RkF}
	if(!is.null(x$R1R)) {cat("\nR1R  = ",round(x$R1R,digits),"Generalizability of a single time point across all items (Random time effects)")
		 mat["R1R"] <- x$R1R}
	
	if(!is.null(x$RkR)) {cat("\nRkR  = ",round(x$RkR,digits),"Generalizability of average time points across all items (Random time effects)")
	  mat["RkR"] <- x$RkR}
	if(!is.null(x$Rc))  {cat("\nRc   = ",round(x$Rc,digits),"Generalizability of change (fixed time points, fixed items) ")
	   mat["Rc"] <- x$Rc}
	if(!is.null(x$RkRn) ) {cat("\nRkRn = ",round(x$RkRn,digits),"Generalizability of between person differences averaged over time (time nested within people)")
	    mat["RkRn"] <- x$RkRn}
	if(!is.null(x$Rcn)) {cat("\nRcn  = ",round(x$Rcn,digits),"Generalizability of within person variations averaged over items  (time nested within people)")
	   mat["Rcn"] <- x$Rcn}
	
	if(!x$lmer && !is.null(x$RkF) ) {cat("\n\n These reliabilities are derived from the components of variance estimated by ANOVA \n") 
		if(!is.null(x$components)) {
	     if(!any(is.na(x$components[1:8,1])) & any(x$components[1:8,1] < 0))  { warning("The ANOVA based estimates are suspect, probably due to missing data, try using lmer")}
	        }} else {
	 if(x$lmer ) { cat("\n\n These reliabilities are derived from the components of variance estimated by lmer \n")}}
	if(!is.null(x$components) && !is.na(x$components[1,1]))  {      print(round(x$components[1:8,],digits=digits))
	mat["components"] <- list( x$components[1:8,])}
		if(!is.null(x$components) && !is.na(x$components[9,1] )) {
		if(!x$lmer) {cat("\n The nested components of variance estimated from lme are:\n")} else {cat("\n The nested components of variance estimated from lmer are:\n")}
	   print(x$components[9:12,],digits=digits)
	   mat["lmer"] = list(x$components[9:12,])} else {cat("\nNested components were not found because lme was not used\n")}
	
	
	if(!short) {cat("\n\n Three way ANOVA or lmer analysis \n")
	print(x$ANOVA,digits=digits)
	 cat("\nvariance components from lme(r)\n")
	 print(x$s.lme,digits=digits)
	 cat("\n Alpha reliability by subjects)\n")
	 print(x$alpha,digits) }
	
	 if(all) {
	          cat("\n Intraclass Correlations by subjects (over time and items) \n")
	          print(x$summary.by.person,digits)
	          cat("\n Intraclass Correlations by time (over subjects and items) \n")
	          print(x$summary.by.time,digits) }
	    
	         if(short) { cat("\nTo see the ANOVA and alpha by subject, use the short = FALSE option.")}
	         if(!all) {cat("\n To see the summaries of the ICCs by subject and time, use all=TRUE")}
	         cat("\n To see specific objects select from the following list:\n",names(x)[-c(1:10)])
	invisible(mat)
	         
	}
	
