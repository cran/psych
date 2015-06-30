#developed July 4, 2012
#modified July 9, 2014 to allow polychorics within groups
#modified June 2, 2015 to include covariance, pearson, spearman, poly ,etc. in correlations
#some ideas taken from Bliese multilevel package (specifically, the WABA results)
"statsBy" <-
   function (data,group,cors=FALSE, cor="cor", method="pearson",use="pairwise", poly=FALSE,na.rm=TRUE) { #  
    cl <- match.call()
  valid <- function(x) { #count the number of valid cases 
        sum(!is.na(x))
    }
#defing a function to count pairwise observations
  pairwise <- function(x) {n <- t(!is.na(x)) %*% (!is.na(x))
              n}
#get the grouping information           
gr <- which(colnames(data) %in% group)      
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
               
               rownn <- lapply(temp,is.null)
               if(sum(as.integer(rownn)) > 0) {
              	 rown <-  names(temp)[-which(rownn==TRUE)] } else {rown <- names(temp) }              
               xvals$mean <- t(matrix(unlist(temp),nrow=ncol(data)))              
               xvals$sd <-t(matrix(unlist(by(data,z,function(x) sapply(x,sd,na.rm=na.rm))),nrow=ncol(data)))
               xvals$n <- t(matrix(unlist(by(data,z,function(x) sapply(x,valid))),nrow=ncol(data)))
              
               
               colnames(xvals$mean) <- colnames(xvals$sd) <- colnames(xvals$n) <-  colnames(data)
               rownames(xvals$mean) <-  rownames(xvals$sd) <- rownames(xvals$n) <- rown
                              nH <- harmonic.mean(xvals$n)
               nG <- nrow(xvals$mean)         #we need to fix this so it is just for the cells that are not NA
               GM <- colSums(xvals$mean*xvals$n,na.rm=na.rm)/colSums(xvals$n,na.rm=na.rm) 
               MSb <- colSums(xvals$n*t((t(xvals$mean) - GM)^2),na.rm=na.rm)/(nG-1) #weight means by n
               MSw <- colSums(xvals$sd^2*(xvals$n-1),na.rm=na.rm)/(colSums(xvals$n-1))#find the pooled sd
               
               xvals$F <- MSb/MSw
               N <- colSums(xvals$n)
             
              npr <- (colSums(xvals$n-1)+nrow(xvals$n))/(nrow(xvals$n))
               xvals$ICC1 <- (MSb-MSw)/(MSb + MSw*(npr-1))
               xvals$ICC2 <- (MSb-MSw)/(MSb)
               

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
       mixed = {r <- by(data,z,function(x) mixed.cor(x[-gr])$rho)}
       )         
                 
              nvars <-  ncol(r[[1]])
              xvals$r <- r   #store them as square matrices
              lower <- lapply(r,function(x) x[lower.tri(x)])
              xvals$within <- t(matrix(unlist(lower),nrow=nvars*(nvars-1)/2))  #string them out as well
             
              cnR <- abbreviate(cnames[-gr],minlength=5) 
             
      k <- 1
       colnames(xvals$within) <- paste("V",1:ncol(xvals$within))
      for(i in 1:(nvars-1)) {for (j in (i+1):nvars) {
      	colnames(xvals$within)[k] <- paste(cnR[i],cnR[j],sep="-")
     	 k<- k +1 }} 
             
             wt <- by(data,z,function(x) count.pairwise(x[-gr]))
             lower.wt <- t(matrix(unlist(lapply(wt,function(x) x[lower.tri(x)])    )  ,nrow=nvars*(nvars-1)/2))
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
             
             diffs <- new.data[,(nvar+1):ncol(new.data)] - new.data[,1:nvar]
             colnames(diffs) <- paste(colnames(new.data)[(nvar + 1):ncol(new.data)], ".wg", sep = "")
             xvals$rbg <- cor(new.data[,1:nvar],use="pairwise",method=method)  #the between group (means)
             t <- (xvals$rbg*sqrt(nG-2))/sqrt(1-xvals$rbg^2)
             
             if(nG > 2) {xvals$pbg <- 2*(1 - pt(abs(t),(nG-2)))} else {xvals$pbg <- NA}
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
             
            xvals$pwg <- 2*(1 - pt(abs(t),(N - nG -2)))
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
