cluster.fit <- function(original,load,clusters,diagonal=FALSE) {
  Pattern <- TRUE 
  df <- nrow(original) * (ncol(original)-1)/2   

 sqoriginal <- original*original    #squared correlations
 totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
        load <- as.matrix(load) 
        clusters <- as.matrix(clusters)
        model <- load %*% t(load)               #reproduce the correlation matrix by the factor law R=  FF' except these are not orthogonal 
        residual <- original-model              #find the residual  R* = R - FF'
        sqresid <- residual*residual            #square the residuals
        totalresid <- sum(sqresid)- diagonal * sum(diag(sqresid) )      #sum squared residuals - the main diagonal
        fit <- 1-totalresid/totaloriginal       #fit is 1-sumsquared residuals/sumsquared original     (of off diagonal elements)
        
        
         covar <- t(clusters) %*% original %*% clusters    #matrix algebra is our friend
       
        phi <- cov2cor(covar)  
        phi.inv <- try(solve(phi),TRUE)
        if(class(phi.inv) == as.character("try-error")) {Pattern <- FALSE
                    message("Could not invert cluster intercorrelation matrix, pattern matrix not found")
                    } #can not invert matrix
    
        if(Pattern) {
      pattern <- load %*%  phi.inv
       model2 <- pattern  %*% t(load) 
       residual <- original - model2
       sqresid <- residual*residual
       totalresid <- sum(sqresid) -(1-diagonal) * sum(diag(sqresid))   #changed Sept 2, 2012 to make more sense (i.e. don't count diagonal if diagonal is false)
       patternrmse <- sqrt(totalresid/(2*df))
       fit2 <- 1-totalresid/totaloriginal } else {fit2 <- NULL 
           patternrmse <- 0} 
        
       clusters <- abs(clusters)               #why do I do this?
       model.1 <- (load * clusters) %*% phi %*% t(load*clusters) #because the items are already signed
       residual <- original - model.1
       sqresid <- residual*residual            #square the residuals
        totalresid <- sum(sqresid)- diagonal * sum(diag(sqresid) )      #sum squared residuals - the main diagonal
        fit.1 <- 1-totalresid/totaloriginal       #fit is 1-sumsquared residuals/sumsquared original     (of off diagonal elements 
        clusterrmse <- sqrt(totalresid/(2*df))
     
        
 cluster.fit <- list(clusterfit=fit.1,structurefit=fit,patternfit=fit2,clusterrmse=clusterrmse,patternrmse=patternrmse)
 }
 