#ICLUST  - a function to form homogeneous item composites
# originally based upon Revelle, W. (1979). Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74.
# but much changed over the years
# pseudo code
#	find similarity matrix
#		original is either covariance or correlation
#		corrected is disattenuated
#find most similar pair
#if size of pair > min size, apply beta criterion
#	if beta new > min(beta 1, beta 2) combine pair
#update similarity matrix
#repeat until finished
#then report various summary statistics


#ICLUST is the main function and calls other routines
"ICLUST" <- 
 function (r.mat,nclusters=0,alpha=3,beta=1,beta.size=4,alpha.size=3,correct=TRUE, 
 correct.cluster=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2,labels=NULL,cut=0,n.iterations=0,title="ICLUST",plot=TRUE,weighted=TRUE,cor.gen =TRUE,SMC=TRUE,purify=TRUE,diagonal=FALSE ) {
iclust(r.mat,nclusters,alpha,beta,beta.size,alpha.size,correct,correct.cluster,reverse,beta.min,output,digits,labels,cut,n.iterations,title,plot,weighted,cor.gen,SMC,purify,diagonal)}


"iclust" <- 
 function (r.mat,nclusters=0,alpha=3,beta=1,beta.size=4,alpha.size=3,correct=TRUE,correct.cluster=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2,labels=NULL,cut=0,n.iterations=0,title="ICLUST",plot=TRUE,weighted=TRUE,cor.gen =TRUE,SMC=TRUE,purify=TRUE,diagonal=FALSE ) {#should allow for raw data, correlation or covariances

 #ICLUST.options <- list(n.clus=1,alpha=3,beta=1,beta.size=4,alpha.size=3,correct=TRUE,correct.cluster=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2,cor.gen=TRUE) 
 cl <- match.call()
 if(is.null(labels)) {labels <- colnames(r.mat)} else {if((length(labels)==1) && (!labels)) labels <- NULL}  #modified October 11, 2011
 
 ICLUST.debug <- FALSE
 
	ICLUST.options <- list(n.clus=nclusters,alpha=alpha,beta=beta,beta.size=beta.size,alpha.size=alpha.size,correct=correct,correct.cluster=correct.cluster,reverse=reverse,beta.min=beta.min,output=output,digits=digits,weighted=weighted,cor.gen=cor.gen,SMC=SMC) 
	
	if(dim(r.mat)[1]!=dim(r.mat)[2]) {r.mat <- cor(r.mat,use="pairwise") } else {r.mat <- cov2cor(r.mat)}    #cluster correlation matrices, find correlations if not square matrix -- added the conversion from covariances to correlations, March, 2012
	if(!is.matrix(r.mat)) {r.mat <- as.matrix(r.mat)}    # for the case where we read in a correlation matrix as a data.frame
	nvar <- dim(r.mat)[2] 
	if(nvar < 3 ) {message("Cluster analysis of items is only meaningful for more than 2 variables. Otherwise, you will find one cluster that is just the composite of the two.  Beta = Alpha = 2*r/(1+r).  Have you made a mistake? \n Try calling the alpha function to give some trivial statistics.")
	  stop() }
	if(is.null(colnames(r.mat))) {colnames(r.mat) <- paste("V",1:nvar,sep="")} 
	
	if(is.null(rownames(r.mat))) {rownames(r.mat) <- paste("V",1:nvar,sep="")} 
	
	#added 24/4/15 to check for bad data
	if(any(is.na(r.mat))) {
       bad <- TRUE
      tempr <-r.mat
      wcl <-NULL
     while(bad) {
     	wc <- table(which(is.na(tempr), arr.ind=TRUE))  #find the correlations that are NA
    	wcl <- c(wcl,as.numeric(names(which(wc==max(wc)))))
    	tempr <- r.mat[-wcl,-wcl]
    	if(any(is.na(tempr))) {bad <- TRUE} else {bad <- FALSE}
         }

     	cat('\nLikely variables with missing values are ',colnames(r.mat)[wcl],' \n')
      	 stop("I am sorry: missing values (NAs) in the correlation matrix do not allow me to continue.\nPlease drop those variables and try again." )
       }
       
	#added this 30/12/13 to improve speed
	if(ICLUST.options$SMC) {smc.items <- smc(r.mat)} else {smc.items <- rep(1,nvar)}  
	
	iclust.results <- ICLUST.cluster(r.mat,ICLUST.options,smc.items) #ICLUST.cluster does all the work - the answers are in iclust.results
	
	loads <- cluster.loadings(iclust.results$clusters,r.mat,SMC=SMC) #summarize the results by using cluster.loadings  -- these are the original values
	
	if(is.matrix(iclust.results$clusters) ) {
		eigenvalue <- diag(t(loads$pattern) %*% loads$loading)
		sorted.cluster.keys.ord <- order(eigenvalue,decreasing=TRUE)
		sorted.cluster.keys <- iclust.results$clusters[,sorted.cluster.keys.ord]
	loads <- cluster.loadings(sorted.cluster.keys,r.mat,SMC=SMC)  #these are the original cluster loadings with clusters sorted by eigenvalues   
	iclust.results$clusters <- sorted.cluster.keys
	cluster.beta <- iclust.results$results[colnames(sorted.cluster.keys),"beta"]
		names(cluster.beta) <- colnames(sorted.cluster.keys)} else {sorted.cluster.keys <- iclust.results$clusters} #these are fine
	
	fits <- cluster.fit(r.mat,as.matrix(loads$loadings),iclust.results$clusters,diagonal)  #check this 
	sorted <- ICLUST.sort(ic.load=loads,labels=labels,cut=cut) #sort the loadings (again?) This is done for sorted output if desired
	
	if(is.matrix(sorted.cluster.keys) ) {cluster.beta <- iclust.results$results[colnames(sorted.cluster.keys),"beta"]
		names(cluster.beta) <- colnames(sorted.cluster.keys) 
		} else {
			number.of.clusters <- dim(iclust.results$results)[1]
			cluster.beta <- iclust.results$results[number.of.clusters,"beta"]}
	
	#now, iterate the cluster solution to clean it up (if desired)
	
		clusters <- as.matrix(iclust.results$clusters)     #just in case there is only one cluster  -- these are now sorted by eigen value
		
		if (dim(clusters)[2]==0 ) {warning('no items meet the specification time1')}
		old.clusters <- clusters
		old.fit <- fits$clusterfit
		if (ICLUST.debug) {print(paste('clusters ',clusters))}
		if(purify) {clusters <- factor2cluster(loads,cut=cut)    #this will assign items to the clusters with the highest loadings --  might be different from original solution
			clusters <- as.matrix(clusters) }  #in case only one cluster 
		if (dim(clusters)[2]==0 ) {warning('no items meet the specification stage 2',immediate.=TRUE)}
		if (ICLUST.debug) {print(paste('clusters ',clusters))
		 print(paste('loads ',loads))}
		loads <- cluster.loadings(clusters,r.mat,SMC=SMC) 
		
		if (n.iterations > 0) {  #it is possible to iterate the solution to perhaps improve it 
		for (steps in 1:n.iterations) {   #
			loads <- cluster.loadings(clusters,r.mat,SMC=SMC)
			
			clusters <- factor2cluster(loads,cut=cut)
			if(dim(clusters)[2]!=dim(old.clusters)[2]) {change <- 999 
			loads <- cluster.loadings(clusters,r.mat,SMC=SMC) 
			 } else {
			change <- sum(abs(clusters)-abs(old.clusters)) }  #how many items are changing?
			fit <- cluster.fit(r.mat,as.matrix(loads$loadings),clusters,diagonal)
		old.clusters <- clusters
		print(paste("iterations ",steps," change in clusters ", change, "current fit " , fit$clusterfit))
		if ((abs(change) < 1) | (fit$clusterfit <= old.fit)) {break}    #stop iterating if it gets worse or there is no change in cluster definitions
		old.fit <- fit$cluster.fit
					}
		}
   
	p.fit <- cluster.fit(r.mat,as.matrix(loads$loadings),clusters,diagonal)
	p.sorted <- ICLUST.sort(ic.load=loads,labels=labels,cut=cut,keys=TRUE)   #at this point, the clusters have been cleaned up, but are not in a sorted order.  Sort them
	
	purified <- cluster.cor(p.sorted$clusters,r.mat,SMC=SMC,item.smc=smc.items)
	class(loads$loadings) <- "loading"
	result <- list(title=title,clusters=iclust.results$clusters,corrected=loads$corrected,loadings=loads$loadings,pattern=loads$pattern,G6 = loads$G6,fit=fits,results=iclust.results$results,cor=loads$cor,Phi=loads$cor,alpha=loads$alpha,beta=cluster.beta,av.r = loads$av.r,size=loads$size,
	sorted=sorted,
	p.fit = p.fit,p.sorted = p.sorted,purified=purified,purify=purify,call=cl)
	#if(plot && requireNamespace('Rgraphviz')) {ICLUST.rgraph(result,labels=labels,title=title,digits=digits)}
	if(plot) iclust.diagram(result,labels=labels,main=title,digits=digits)
	class(result) <- c("psych","iclust")

	return(result)
}   


 

 
