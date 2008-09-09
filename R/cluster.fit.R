cluster.fit <- function(original,load,clusters,diagonal=FALSE,digits=2) {
      

 sqoriginal <- original*original    #squared correlations
 totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
        load <- as.matrix(load) 
        clusters <- as.matrix(clusters)
        model <- load %*% t(load)               #reproduce the correlation matrix by the factor law R=  FF'
        residual <- original-model              #find the residual  R* = R - FF'
        sqresid <- residual*residual            #square the residuals
        totalresid <- sum(sqresid)- diagonal * sum(diag(sqresid) )      #sum squared residuals - the main diagonal
        fit <- 1-totalresid/totaloriginal       #fit is 1-sumsquared residuals/sumsquared original     (of off diagonal elements
        
        
         covar <- t(clusters) %*% original %*% clusters    #matrix algebra is our friend
        var <- diag(covar)
        sd.inv <- 1/sqrt(var)
       ident.sd <- diag(sd.inv,ncol = length(sd.inv))
    phi <- ident.sd %*% covar  %*% ident.sd
       phi.inv <- solve(phi)
       pattern <- load %*%  phi.inv
       model2 <- pattern  %*% t(load)
       residual <- original - model2
       sqresid <- residual*residual
       totalresid <- sum(sqresid) -diagonal * sum(diag(sqresid))
       fit2 <- 1-totalresid/totaloriginal 
        
       clusters <- abs(clusters)               #why do I do this?
       model.1 <- (load * clusters) %*%  t(load*clusters)
       residual <- original - model.1
       sqresid <- residual*residual            #square the residuals
        totalresid <- sum(sqresid)- diagonal * sum(diag(sqresid) )      #sum squared residuals - the main diagonal
        fit.1 <- 1-totalresid/totaloriginal       #fit is 1-sumsquared residuals/sumsquared original     (of off diagonal elements 
        
        
 cluster.fit <- list(clusterfit=round(fit.1,digits),structurefit=round(fit,digits),patternfit =round(fit2,digits))
 }
 