#developed July 4, 2012
#modified July 9, 2014 to allow polychorics within groups
#modified June 2, 2015 to include covariance, pearson, spearman, poly ,etc. in correlations
#Fixed March 3, 2017 to not weight empty cells in finding ICCs
#some ideas taken from Bliese multilevel package (specifically, the WABA results)
#modifed November 3, 2018 to allow a single DV (in case cohen.d is just comparing two groups on one DV)
#corrected January 1, 2019 to allow for grouping variables that are characters
"statsBy" <-
   function (data,
            group,
            cors=FALSE,
            cor="cor",
             method="pearson",
             use="pairwise", 
             poly=FALSE,
             na.rm=TRUE,
             alpha=.05,
             minlength=5) { #  
 cl <- match.call()
    
valid <- function(x) { #count the number of valid cases 
        sum(!is.na(x))
    }
    
#define a function to count pairwise observations
 pairwise <- function(x) {n <- t(!is.na(x)) %*% (!is.na(x))
              n}
#get the grouping information

if(length(group) < NROW(data) ){   #added 01/01/19 to handle the case of non-numeric grouping data
if(is.character(group) ) {
gr <- which(colnames(data) %in% group) } else {gr <- group} 
}   else {data <- cbind(data,group)
    group <- "group"
    gr <- which(colnames(data) %in% group)}   
z1 <- data[,group]
       z <- z1
       cnames <- colnames(data)
       for (i in 1:ncol(data)) {if(is.factor(data[,i]) || is.logical(data[,i])) {
             data[,i] <- as.numeric(data[,i])
            # colnames(data)[i] <- paste(cnames[i],"*",sep="")
             }}
       xvals <- list()
       #find the statistics by group
               temp <- by(data,z,colMeans,na.rm=na.rm)
               rowname <- dimnames(temp)[[1]]
               if(length(dimnames(temp))> 1){            #if we have multiple criteria, we need to name them
               for (i in 2:length(dimnames(temp))) {
               rowname <- paste0(rep(rowname,each=length(dimnames(temp)[[i]])),"-",dimnames(temp)[[i]])} 
               }

               rownn <- lapply(temp,is.null)     #drop rows without values
               if(sum(as.integer(rownn)) > 0) {
              	 rowname <-  rowname[-which(rownn==TRUE)] }  # look for missing criteria           
               xvals$mean <- t(matrix(unlist(temp),nrow=ncol(data)))              
               xvals$sd <-t(matrix(unlist(by(data,z,function(x) sapply(x,sd,na.rm=na.rm))),nrow=ncol(data)))
               xvals$n <- t(matrix(unlist(by(data,z,function(x) sapply(x,valid))),nrow=ncol(data)))
              
               
               colnames(xvals$mean) <- colnames(xvals$sd) <- colnames(xvals$n) <-  colnames(data)
               rownames(xvals$mean) <-  rownames(xvals$sd) <- rownames(xvals$n) <- rowname
                             # nH <- harmonic.mean(xvals$n)     #this will be 0 if any is zero -- but is not used anyway
               nG <- colSums(!is.na(xvals$mean))         #fixed this so it is just for the cells that are not NA
               GM <- colSums(xvals$mean*xvals$n,na.rm=na.rm)/colSums(xvals$n,na.rm=na.rm) 
               MSb <- colSums(xvals$n*t((t(xvals$mean) - GM)^2),na.rm=na.rm)/(nG-1) #weight means by n
               MSw <- colSums(xvals$sd^2*(xvals$n-1*(xvals$n>0)),na.rm=na.rm)/(colSums(xvals$n-1*(xvals$n>0)))#find the pooled sd   #fix this for 0 cell size
               
               xvals$F <- MSb/MSw
               N <- colSums(xvals$n)  #overall N for each variable
             
              npr <- (colSums(xvals$n-1*(xvals$n > 0))+colSums(xvals$n >0))/(colSums(xvals$n >0))
               xvals$ICC1 <- (MSb-MSw)/(MSb + MSw*(npr-1))
               xvals$ICC2 <- (MSb-MSw)/(MSb)
               
               #now, figure out the cis for the ICCs  
               #taken from the ICC function
                F11 <- MSb/MSw
  #df11n <- n.obs-1
  #df11d <- n.obs*(nj-1)
  df11n <- nG - 1
  df11d <- nG* (npr-1) 
 
 # p11 <- 1-pf(F11,df11n,df11d)
    p11 <-  -expm1(pf(F11,df11n,df11d,log.p=TRUE))  
  #F21 <- MSB/MSE
  df21n <- N - 1
  df21d <-  N * (npr-1) 
 # p21 <- 1-pf(F21,df21n,df21d)
 # F31 <- F21
 F1L <- F11 / qf(1-alpha/2,df11n,df11d)  
 F1U <- F11 * qf(1-alpha/2,df11d,df11n)
 L1 <- (F1L-1)/(F1L+(npr-1))
 U1 <- (F1U -1)/(F1U+(npr-1))
 L2 <- 1-1/F1L
 U2 <- 1 -1/F1U
 xvals$ci1 <-as.matrix(data.frame(L1 = L1,U1 = U1) ) 
 xvals$ci2 <- as.matrix(data.frame(L2 = L2,U2 = U2) )    
 
    #if we want within group correlations, then find them  
     # if(cors) {if(!poly) { r <- by(data,z,function(x) cor(x[-gr],use="pairwise",method=method)) } else { r <- by(data,z,function(x) polychoric(x[-gr])$rho)}
               
      #added 02/06/15
       if(cors) {if (poly) {cor <- "poly"}
       switch(cor, 
       cor = {r <- by(data,z,function(x) cor(x[-gr],use=use,method=method))},
       cov = {r <- by(data,z,function(x) cov(x[-gr],use=use))
              covar <- TRUE},
       tet = {r <- by(data,z,function(x) tetrachoric(x[-gr])$rho)},
       poly = {r <- by(data,z,function(x) polychoric(x[-gr])$rho)},
       mixed = {r <- by(data,z,function(x) mixedCor(x[-gr])$rho)}
       )         
            
             nWg <- by(data,z,function(x) pairwise(x[-gr, -gr]))
              nvars <-  ncol(r[[1]])
              xvals$r <- r   #store them as square matrices
             
              length.r <- length(r)
              # xvals$r.ci <- r
             
              xvals$r.ci <- mapply(cor.Ci, r=r, n=nWg)
              attributes(xvals$r.ci) <- attributes(r)
              lower <- lapply(r,function(x) if(!is.null(x)){ x[lower.tri(x)]})
              xvals$within <- t(matrix(unlist(lower),nrow=nvars*(nvars-1)/2))  #string them out as well
             
              cnR <- abbreviate(cnames[-gr],minlength=minlength) 
           
      k <- 1
       colnames(xvals$within) <- paste("V",1:ncol(xvals$within))
      for(i in 1:(nvars-1)) {for (j in (i+1):nvars) {
      	colnames(xvals$within)[k] <- paste(cnR[i],cnR[j],sep="-")
     	 k<- k +1 }}
     	  
     	# rownames(xvals$within) <- paste0("z",names(xvals$r))
          rownames(xvals$within) <- rowname
             
             wt <- by(data,z,function(x) pairwiseCount(x[-gr]))
             lower.wt <- t(matrix(unlist(lapply(wt,function(x) if(!is.null(x)) { x[lower.tri(x)]})    )  ,nrow=nvars*(nvars-1)/2))
             lower.wt <- t(t(lower.wt)/colSums(lower.wt,na.rm=TRUE))
             pool  <- colSums( lower.wt * xvals$within,na.rm=TRUE)
             pool.sd <- apply(xvals$within, 2,FUN=sd, na.rm=TRUE)
             xvals$pooled <- matrix(0,nvars,nvars)
             xvals$pooled[lower.tri(xvals$pooled)] <- pool  
             xvals$pooled <- xvals$pooled + t(xvals$pooled)  #changed, May 12 to properly reflect values
             diag(xvals$pooled) <- 1
             xvals$sd.r <-  matrix(NaN,nvars,nvars)
             xvals$sd.r[lower.tri(xvals$sd.r)] <- pool.sd
             xvals$sd.r[upper.tri(xvals$sd.r)] <- pool.sd
             colnames(xvals$pooled) <- rownames (xvals$pooled) <- cnames[-gr]
              }

              nvar <- ncol(data)-length(group) #we have dropped the grouping variable
            #   if(!poly) {xvals$raw <- cor(data,use="pairwise",method=method)} else {xvals$raw <- polychoric(data)$rho}
               
        ##added 02/06/15
         if (poly) cor <- "poly"
            switch(cor, 
       			cor = {xvals$raw  <- cor(data,use=use,method=method)},
       			cov = {xvals$raw  <- cov(data,use=use) 
              covar <- TRUE},
             poly=  {xvals$raw <- polychoric(data)$rho},
             tet = {xvals$raw <- tetrachoric(data)$rho},
             mixed = {xvals$raw <- mixed.cor(data)$rho}   
       )
       
             new.data <- as.matrix( merge(xvals$mean,data,by=group,suffixes =c(".bg",""))) #drop the grouping variable(s) 
             new.data <- new.data[,(length(group)+1):ncol(new.data)]
             
             diffs <- new.data[,(nvar+1):ncol(new.data),drop=FALSE] - new.data[,1:nvar]  #Difference of raw data within group - within group mean
             #why don't we just use scale for the data within each group?

             colnames(diffs) <- paste(colnames(new.data)[(nvar + 1):ncol(new.data)], ".wg", sep = "")
             if(nvar > 1) {xvals$rbg <- cor(new.data[,1:nvar],use="pairwise",method=method)} else {xvals$rbg <- NA}  #the between group (means)
             nBg <- pairwise(xvals$mean)[-gr,-gr]
            if(any(nBg > 2)) {xvals$ci.bg <-cor.Ci(xvals$rbg,nBg,alpha=alpha,minlength=minlength)} else {xvals.ci.bg=NA}  #added the test for all nBg  <= 2 so we don't call cor.Ci and throw an error  11/2/19
             #t <- rep(NA,(length(nG)-length(gr)))
             
             if(all(nG > 2)) {
             t <- (xvals$rbg*sqrt(nG[-gr]-2))/sqrt(1-xvals$rbg^2)
            # if(any(nG < 3) ) {#warning("Number of groups must be at least 2")
            
              xvals$pbg <- 2*(1 - pt(abs(t),(nG-2)))
              } else { t <- xvals$pbg <- NA }
                      #     xvals$rwg <- cor(diffs,use="pairwise",method=method)  #the within group (differences)
            if(cor %in% c("tet","poly","mixed","mixed.cor") ) cor <- "cor"
            switch(cor, 
      			 cor = {xvals$rwg  <- cor(diffs,use=use,method=method)},
      			 cov = {xvals$rwg  <- cov(diffs,use=use) 
              covar <- TRUE}
       )
       
               xvals$nw <- pairwise(diffs)
               rwg <- cov2cor(xvals$rwg)
               t <- (rwg*sqrt(xvals$nw -2))/sqrt(1-rwg^2)
                if(any(NCOL(rwg) > 2)) {xvals$ci.wg <- cor.Ci(rwg,xvals$nw,alpha=alpha,minlength=minlength)} else {xvals.ci.wg=NA}  #added the test for all nBg  <= 2 so we don't call cor.Ci and throw an error  11/2/19
              # xvals$ci.wg <- cor.Ci(rwg,xvals$nw,alpha=alpha,minlength=minlength)   #this is a local function
               
               
            if(all(nG > 2)) {
            xvals$pwg <- 2*(1 - pt(abs(t),(N[-gr] - nG[-gr] -2)))} else {xvals$pwg <- NA}
            # colnames(xvals$rwg) <- rownames(xvals$rwg) <- paste(colnames(xvals$rwg),".wg",sep="")
             xvals$etabg <- diag(cor(new.data[,1:(nvar)],new.data[,(nvar+1):ncol(new.data)],use="pairwise",method=method) )#the means with the data
             xvals$etawg <- diag(cor(new.data[,(nvar+1):ncol(new.data)],diffs,use="pairwise",method=method)) #the deviations and the data
            names(xvals$etabg)  <- colnames(xvals$rbg)
         
            xvals$nwg <- N - nG
            xvals$nG <- nG
            
            xvals$Call <- cl
    statsBy <- xvals
    class(statsBy) <- c("psych","statsBy")
    return(statsBy)
    }
    
cor.Ci <- function(r,n,alpha=.05,minlength=10) {
   cl <- match.call()

z.r <- fisherz(r)
sd <-1/sqrt(n-3)
nvar <- NCOL(r)
ci.lower <- fisherz2r(z.r  +( qnorm(alpha/2) * sd))
ci.upper <- fisherz2r(z.r + (qnorm(1-alpha/2) * sd))
cnR <- abbreviate(colnames(r),minlength=minlength)
ci.lower <- ci.lower[lower.tri(ci.lower)]
ci.upper <- ci.upper[lower.tri(ci.upper)]
r.ci <- data.frame(lower=ci.lower,r=r[lower.tri(r)],upper=ci.upper)
k <- 1
 for(i in 1:(nvar-1)) {for (j in (i+1):nvar) {
      rownames(r.ci)[k] <- paste(cnR[i],cnR[j],sep="-")
 k <- k + 1     
}}
result <-list(r.ci=r.ci)
class(result) <- cs(psych,corCi )
return(result)
}
