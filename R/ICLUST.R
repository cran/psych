#ICLUST  - a function to form homogeneous item composites
# based upon Revelle, W. (1979). Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74.
#
# psudo code
#	find similarity matrix
#		original is either covariance or correlation
#		corrected is disattenuated
#find most similar pair
#if size of pair > min size, apply beta criterion
#	if beta new > min(beta 1, beta 2) combine pair
#update similarity matrix
#repeat until finished
#then report various summary statistics
#example code
#r.mat<- Harman74.cor$cov
# print(ICLUST(r.mat),digits=2)

#ICLUST is the main function and calls other routines

"ICLUST" <- 
 function (r.mat,nclusters=0,alpha=3,beta=1,beta.size=4,alpha.size=3,correct=TRUE,correct.cluster=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2,labels=NULL,cut=0,n.iterations=0,title="ICLUST",plot=TRUE) {#should allow for raw data, correlation or covariances

 #ICLUST.options <- list(n.clus=1,alpha=3,beta=1,beta.size=4,alpha.size=3,correct=TRUE,correct.cluster=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2) 
 cl <- match.call()
 if(is.null(labels)) {labels <- colnames(r.mat)} else {if(!labels) labels<- NULL}
 ICLUST.debug <- FALSE
	ICLUST.options <- list(n.clus=nclusters,alpha=alpha,beta=beta,beta.size=beta.size,alpha.size=alpha.size,correct=correct,correct.cluster=correct.cluster,reverse=reverse,beta.min=beta.min,output=output,digits=digits) 
	if(dim(r.mat)[1]!=dim(r.mat)[2]) {r.mat <- cor(r.mat,use="pairwise") }    #cluster correlation matrices, find correlations if not square matrix
	if(!is.matrix(r.mat)) {r.mat <- as.matrix(r.mat)}    # for the case where we read in a correlation matrix as a data.frame
	iclust.results <- ICLUST.cluster(r.mat,ICLUST.options) #this does all the work - the answers are in iclust.results
	loads <- cluster.loadings(iclust.results$clusters,r.mat) #summarize the results by using cluster.loadings
	
	if(is.matrix(iclust.results$cluster) ) {
		eigenvalue <- diag(t(loads$loading) %*% loads$loading)
		sorted.cluster.keys.ord <- order(eigenvalue,decreasing=TRUE)
		sorted.cluster.keys <- iclust.results$clusters[,sorted.cluster.keys.ord]
	loads <- cluster.loadings(sorted.cluster.keys,r.mat) } else {sorted.cluster.keys <- iclust.results$clusters}
	
	fits <- cluster.fit(r.mat,as.matrix(loads$loadings),iclust.results$clusters)
	sorted <- ICLUST.sort(ic.load=loads,labels=labels,cut=cut) #sort the loadings (again?  I think this might not be necessary anymore
	
	if(is.matrix(sorted.cluster.keys) ) {cluster.beta <- iclust.results$results[colnames(sorted.cluster.keys),"beta"]
	names(cluster.beta) <- colnames(sorted.cluster.keys) } else {
	number.of.clusters <- dim(iclust.results$results)[1]
	cluster.beta <- iclust.results$results[number.of.clusters,"beta"]}
	
	#now, iterate the cluster solution to clean it up (if desired)
	
		clusters <- as.matrix(iclust.results$clusters)     #just in case there is only one cluster
		if (dim(clusters)[2]==0 ) {warning('no items meet the specification time1')}
		old.clusters <- clusters
		old.fit <- fits$clusterfit
		if (ICLUST.debug) {print(paste('clusters ',clusters))}
		clusters <- factor2cluster(loads,cut=cut)
		clusters <- as.matrix(clusters)   #in case only one cluster 
		if (dim(clusters)[2]==0 ) {warning('no items meet the specification stage 2',immediate.=TRUE)}
		if (ICLUST.debug) {print(paste('clusters ',clusters))
		 print(paste('loads ',loads))}
		loads <- cluster.loadings(clusters,r.mat)
		
		if (n.iterations > 0) {  #it is possible to iterate the solution to perhaps improve it 
		for (steps in 1:n.iterations) {   #
			loads <- cluster.loadings(clusters,r.mat)
			
			clusters <- factor2cluster(loads,cut=cut)
			if(dim(clusters)[2]!=dim(old.clusters)[2]) {change <- 999 
			loads <- cluster.loadings(clusters,r.mat) 
			 } else {
			change <- sum(abs(clusters)-abs(old.clusters)) }  #how many items are changing?
			fit <- cluster.fit(r.mat,as.matrix(loads$loadings),clusters)
		old.clusters <- clusters
		print(paste("iterations ",steps," change in clusters ", change, "current fit " , fit$clusterfit))
		if ((abs(change) < 1) | (fit$clusterfit <= old.fit)) {break}    #stop iterating if it gets worse or there is no change in cluster definitions
		old.fit <- fit$cluster.fit
					}
		}

	p.fit <- cluster.fit(r.mat,as.matrix(loads$loadings),clusters)
	p.sorted <- ICLUST.sort(ic.load=loads,labels=labels,cut=cut,keys=TRUE)
	
	purified <- cluster.cor(p.sorted$clusters,r.mat)
	class(loads$loadings) <- "loading"
	result <- list(title=title,clusters=iclust.results$clusters,corrected=loads$corrected,loadings=loads$loadings,pattern=loads$pattern,G6 = loads$G6,fit=fits,results=iclust.results$results,cor=loads$cor,alpha=loads$alpha,beta=cluster.beta,av.r = loads$av.r,size=loads$size,sorted=sorted,p.fit = p.fit,p.sorted = p.sorted,purified=purified,call=cl)
	if(plot && require(Rgraphviz)) {ICLUST.rgraph(result,labels=labels,title=title)}
	class(result) <- c("psych","iclust")
	return(result)
}   


 

 
