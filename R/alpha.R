"alpha" <- 
    function(x,keys=NULL,title=NULL, na.rm=TRUE) {  #find coefficient alpha given a data frame or a matrix
    
    alpha.1 <- function(C,R) {
    n <- dim(C)[2]
    alpha.raw <- (1- tr(C)/sum(C))*(n/(n-1))
    alpha.std <-  (1- n/sum(R))*(n/(n-1))
    smc.R <- smc(R)
    G6 <- (1- (n-sum(smc.R))/sum(R))
    av.r <- (sum(R)-n)/(n*(n-1))
    result <- list(raw=alpha.raw,std=alpha.std,G6,av.r)
    return(result)
    }
    cl <- match.call()
    if(!is.matrix(x) && !is.data.frame(x)) stop('Data must either be a data frame or a matrix')
    nvar <- dim(x)[2]
    nsub <- dim(x)[1]
    if (nsub !=nvar)  {
         C <- cov(x,use="pairwise") 
         if(!is.null(keys)) {
         			keys<- as.vector(keys)
        			 key.d <- diag(keys)
                     C <- key.d %*% C %*% key.d}
         total <- rowSums(x,na.rm=na.rm)
         mean.t <- mean(total,na.rm=na.rm)
         sdev <- sd(total,na.rm=na.rm)
         t.valid <- colSums(!is.na(x))} else { C <- as.matrix(x) 
                                               if(!is.null(keys)) {key.d <- diag(keys)
                            C <- key.d %*% C %*% key.d}}
         R <- cov2cor(C)
         drop.item <-list()
         alpha.total <- alpha.1(C,R)
         if(nvar>2) {
         for (i in 1:nvar) {
         drop.item[[i]] <- alpha.1(C[-i,-i],R[-i,-i])
         } 
         }
        by.item <- data.frame(matrix(unlist(drop.item),ncol=4,byrow=TRUE))
        colnames(by.item) <- c("raw_alpha","std.alpha","G6(smc)","average_r")
        rownames(by.item) <- colnames(x)
       
        Vt <- sum(R)
        item.r <- colSums(R)/sqrt(Vt)
     #correct for item overlap 
        RC <-R
        diag(RC) <- smc(R)
        Vtc <- sum(RC)
        item.rc <-colSums(RC)/sqrt(Vtc)
     #  
        item.means <- colMeans(x, na.rm=na.rm )
        item.sd <-  SD(x,na.rm=na.rm)
        if(nsub > nvar) {
        	alpha.total <- data.frame(alpha.total,mean=mean.t,sd=sdev)
        	colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)","average_r","mean","sd")
        	rownames(alpha.total) <- ""
        	stats <- data.frame(n=t.valid,r =item.r,r.cor = item.rc,mean=item.means,sd=item.sd)} else {
        	alpha.total <- data.frame(alpha.total)
        	        colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)" ,"average_r")
        	        rownames(alpha.total) <- ""

                  stats <- data.frame(r =item.r,r.cor = item.rc)
        	}
       	rownames(stats) <- colnames(x)
        result <- list(total=alpha.total,alpha.drop=by.item,item.stats=stats,call=cl,title=title)
        class(result) <- c("psych","alpha")
        return(result)
     
    }
    
