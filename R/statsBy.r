 #developed July 4, 2012
 #some ideas taken from Bliese multilevel package (specifically, the WABA results)
   "statsBy" <-
   function (data,group,cors=FALSE,method="pearson") { #  
    cl <- match.call()
  valid <- function(x) { #count the number of valid cases 
        sum(!is.na(x))
    }
  pairwise <- function(x) {n <- t(!is.na(x)) %*% (!is.na(x))
              n}
       gr <- which(colnames(data) == group)
       
       z1 <- data[,group]
       z <- z1
       cnames <- colnames(data)
       for (i in 1:ncol(data)) {if(is.factor(data[,i]) || is.logical(data[,i])) {
             data[,i] <- as.numeric(data[,i])
            # colnames(data)[i] <- paste(cnames[i],"*",sep="")
             }}
       xvals <- list()
               temp <- by(data,z,colMeans,na.rm=TRUE)
               
               rownn <- lapply(temp,is.null)
               if(sum(as.integer(rownn)) > 0) {
               rown <-  names(temp)[-which(rownn==TRUE)] } else {rown <- names(temp) }              
               xvals$mean <- t(matrix(unlist(temp),nrow=ncol(data)))              
               xvals$sd <-t(matrix(unlist(by(data,z,function(x) sapply(x,sd,na.rm=TRUE))),nrow=ncol(data)))
               xvals$n <- t(matrix(unlist(by(data,z,function(x) sapply(x,valid))),nrow=ncol(data)))
              
               
               colnames(xvals$mean) <- colnames(xvals$sd) <- colnames(xvals$n) <-  colnames(data)
               rownames(xvals$mean) <-  rownames(xvals$sd) <- rownames(xvals$n) <- rown
                              nH <- harmonic.mean(xvals$n)
               nG <- nrow(xvals$mean)
               GM <- colSums(xvals$mean*xvals$n)/colSums(xvals$n) 
               MSb <- colSums(xvals$n*t((t(xvals$mean) - GM)^2))/(nG-1) #weight means by n
               MSw <- colSums(xvals$sd^2*(xvals$n-1))/(colSums(xvals$n-1)) #find the pooled sd
               xvals$F <- MSb/MSw
               N <- colSums(xvals$n)
              # npr <-(N^2 - colSums(xvals$n))/(N *(nrow(xvals$n) -1))
              # npr <- harmonic.mean(xvals$n-1)
              npr <- (colSums(xvals$n-1)+nrow(xvals$n))/(nrow(xvals$n))
               xvals$ICC1 <- (MSb-MSw)/(MSb + MSw*(npr-1))
               xvals$ICC2 <- (MSb-MSw)/(MSb)
                             if(cors) { r <- by(data,z,function(x) cor(x[-1],use="pairwise",method=method))
              nvars <-  ncol(r[[1]])
              xvals$r <- r
              lower <- lapply(r,function(x) x[lower.tri(x)])
             xvals$within <- t(matrix(unlist(lower),nrow=nvars*(nvars-1)/2))
             wt <- by(data,z,function(x) count.pairwise(x[-1]))
             lower.wt <- t(matrix(unlist(lapply(wt,function(x) x[lower.tri(x)])    )  ,nrow=nvars*(nvars-1)/2))
             lower.wt <- t(t(lower.wt)/colSums(lower.wt,na.rm=TRUE))
             pool  <- colSums( lower.wt * xvals$within,na.rm=TRUE)
             pool.sd <- apply(xvals$within, 2,FUN=sd, na.rm=TRUE)
             xvals$pooled <- matrix(NA,nvars,nvars)
             xvals$pooled[lower.tri(xvals$pooled)] <- pool
             xvals$pooled[upper.tri(xvals$pooled)]  <- pool
             diag(xvals$pooled) <- 1
             xvals$sd.r <-  matrix(NA,nvars,nvars)
             xvals$sd.r[lower.tri(xvals$sd.r)] <- pool.sd
             xvals$sd.r[upper.tri(xvals$sd.r)] <- pool.sd
             colnames(xvals$pooled) <- rownames (xvals$pooled) <- cnames[-1]
              }

              nvar <- ncol(data)-length(group) #we have dropped the grouping variable
               xvals$raw <- cor(data,use="pairwise")
             new.data <- as.matrix( merge(xvals$mean,data,by=group,suffixes =c(".bg",""))) #drop the grouping variable(s) 
             new.data <- new.data[,(length(group)+1):ncol(new.data)]
             
             diffs <- new.data[,(nvar+1):ncol(new.data)] - new.data[,1:nvar]
             colnames(diffs) <- paste(colnames(new.data)[(nvar + 1):ncol(new.data)], ".wg", sep = "")
             xvals$rbg <- cor(new.data[,1:nvar],use="pairwise",method=method)  #the between group (means)
             t <- (xvals$rbg*sqrt(nG-2))/sqrt(1-xvals$rbg^2)
            xvals$pbg <- 2*(1 - pt(abs(t),(nG-2)))
             xvals$rwg <- cor(diffs,use="pairwise",method=method)  #the within group (differences)
             xvals$nw <- pairwise(diffs)

               t <- (xvals$rwg*sqrt(xvals$nw -2))/sqrt(1-xvals$rwg^2)
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
