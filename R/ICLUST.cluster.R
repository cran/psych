ICLUST.cluster <- function (r.mat,ICLUST.options) {#should allow for raw data, correlation or covariances

#options:  alpha =1  (minimum alpha)  2 (average alpha)    3 (maximum alpha)
#          beta =1   (minimum beta)   2 (average beta)     3 (maximum beta)
#          correct  for reliability
#          reverse score items if negative correlations
#          stop clustering if beta for new clusters < beta.min
#          output =1 (short)    2  (show steps)     3 show rejects as we go

#initialize  various arrays and get ready for the first pass
output <- ICLUST.options$output
num.var <- nrow(r.mat)
keep.clustering <- TRUE          #used to determine when we are finished clustering
results <- data.frame(matrix(rep(0,18*(num.var-1)),ncol=18))
names(results) <- c("Item/Cluster", "Item/Cluster","similarity","correlation","alpha1","alpha2",
"beta1","beta2","size1","size2","rbar1","rbar2","r1","r2","alpha","beta","rbar","size")
rownames(results) <- paste("C",1:(num.var-1),sep="")

clusters <- diag(1,nrow =nrow(r.mat))    #original cluster structure is 1 item clusters
rownames(clusters) <- rownames(r.mat)
colnames(clusters) <- paste("V",1:num.var,sep="")
count=1

#master loop 

while (keep.clustering) {   #loop until we figure out we should stop
#find similiarities
 #we will do most of the work on a copy of the r.mat
cluster.stats <- cluster.cor(clusters,r.mat,FALSE)  
sim.mat <- cluster.stats$cor    #the correlation matrix
diag(sim.mat) <- 0   #we don't want 1's on the diagonal to mess up the maximum 

#two ways to estimate reliability -- for 1 item clusters, max correlation, for >1, alpha
#this use of initial max should be an option
if (ICLUST.options$correct) {   #find the largest and smallest similarities for each variable      
	row.range <- apply(sim.mat,1,range,na.rm=TRUE)     
	row.max <- pmax(abs(row.range[1,]),abs(row.range[2,]))  #find the largest absolute similarity
   } else {row.max <- rep(1, nrow(sim)) }         #don't correct for largest similarity
   
   item.rel <-  cluster.stats$alpha
for (i in 1: length(item.rel)) { if (cluster.stats$size[i]<2) {
	item.rel[i] <- row.max[i]
	#figure out item betas here?
	}}
sq.max <- diag(1/sqrt(item.rel))     #used to correct for reliabilities
 #this is the corrected for maximum r similarities
if (ICLUST.options$correct) {sim <- sq.max %*% sim.mat %*% sq.max  
     } else {sim <- sim.mat}
diag(sim) <- NA                  #we need to not consider the diagonal when looking for maxima
#find the most similar pair    and apply tests if we should combine
test.alpha <- FALSE
test.beta <- FALSE

while(!(test.alpha&test.beta)){

max.cell <- which.max(sim)             #global maximum
if (length(max.cell) < 1) {
   keep.clustering <- FALSE
   break}   #there are no non-NA values left
sign.max <- 1
if ( ICLUST.options$reverse ) {            #normal case is to reflect if necessary
	min.cell <- which.min(sim)             #location of global minimum
	if (sim[max.cell] <  abs(sim[min.cell] )) { 
 	  	sign.max <- -1
		max.cell <- min.cell }
	if (sim[max.cell] < 0.0) {sign.max <- -1 }}  
	               #this is a weird case where all the similarities are negative 
	
max.col <- trunc(max.cell/nrow(sim))+1    #is in which row and column?
max.row <- max.cell - (max.col-1)*nrow(sim) #need to fix the case of first column
if (max.row < 1) {max.row <- nrow(sim)
 max.col <- max.col-1 }

#combine these two rows if the various criterion are passed
beta.combined <-  2* sign.max*sim.mat[max.cell]/(1+sign.max* sim.mat[max.cell])   #unweighted beta

size1 <- cluster.stats$size[max.row]
if(size1 < 2) {V1 <- 1
	beta1 <-  item.rel[max.row]
	alpha1 <-  item.rel[max.row]
	rbar1 <- item.rel[max.row] } else {
	rbar1 <- results[cluster.names[max.row],"rbar"]
	beta1 <- results[cluster.names[max.row],"beta"]
	alpha1 <- results[cluster.names[max.row],"alpha"]}
	
 V1 <- size1 + size1*(size1-1) * rbar1
 
size2 <- cluster.stats$size[max.col]
if(size2 < 2) {V2 <- 1 
	beta2 <-  item.rel[max.col] 
	alpha2 <-  item.rel[max.col] 
	rbar2 <- item.rel[max.col] } else {
	rbar2 <- results[cluster.names[max.col],"rbar"]
	beta2 <- results[cluster.names[max.col],"beta"]
	alpha2 <- results[cluster.names[max.col],"alpha"]}

 V2 <- size2 + size2 * (size2-1) * rbar2
	
Cov12 <- sign.max* sim.mat[max.cell] * sqrt(V1*V2)
V12 <- V1 + V2 + 2 * Cov12
size12 <- size1 + size2
alpha <- (V12 - size12)*(size12/(size12-1))/V12
rbar <- alpha/(size12-alpha*(size12-1))

#what is the correlation of this new cluster with the two subclusters?
#this considers item overlap problems
c1 <-  sign.max*rbar1*size1*size1 + sign.max* Cov12   #corrects for item overlap
c2 <-  rbar2*size2*size2 + Cov12     #only flip one of the two correlations with the combined cluster
if(size1 > size2) {r1 <- c1/sqrt(V1*V12)
r2 <- sign.max* c2/sqrt(V2*V12) } else {r1 <-sign.max* c1/sqrt(V1*V12)    
                                     #flip the smaller of the two clusters
r2 <-  c2/sqrt(V2*V12) }

#test if we should combine these two clusters  
#first, does alpha increase?
test.alpha <- TRUE

if (ICLUST.options$alpha>0) { #should we apply any tests?
if (ICLUST.options$alpha.size < min(size1,size2)) {
  switch(ICLUST.options$alpha, {if (alpha < min(alpha1,alpha2)) {if (output>2) {print(
  paste ('do not combine ', cluster.names[max.row],"with", cluster.names[max.col],
  'new alpha =', round (alpha,2),'old alpha1 =',round( alpha1,2),"old alpha2 =",round(alpha2,2)))}
												test.alpha <- FALSE }},
  {if (alpha < mean(alpha1,alpha2)) {if (output>2) {print(paste ('do not combine ',
  cluster.names[max.row],"with", cluster.names[max.col],'new alpha =', round (alpha,2),
  'old alpha1 =',round( alpha1,2),"old alpha2 =",round(alpha2,2)))}
												test.alpha <- FALSE  }},
  {if (alpha < max(alpha1,alpha2)) {if (output>2) {print(paste ('do not combine ',
  cluster.names[max.row],"with", cluster.names[max.col],'new alpha =', round (alpha,2),
  'old alpha1 =',round( alpha1,2),"old alpha2 =",round(alpha2,2)))}
												test.alpha <- FALSE  }}) #end switch  
  }   #end if options$alpha.size
  }

#second, does beta increase ?
test.beta <- TRUE          
if (ICLUST.options$beta>0) { #should we apply any tests?
if (ICLUST.options$beta.size < min(size1,size2)) {
  switch(ICLUST.options$beta, {if (beta.combined < min(beta1,beta2)) {if (output>2) {print(
  paste ('do not combine ', cluster.names[max.row],"with", cluster.names[max.col],'new beta =',
  round (beta.combined,2),'old beta1 =',round( beta1,2),"old beta2 =",round(beta2,2)))}
												test.beta <- FALSE }},
  {if (beta.combined < mean(beta1,beta2)) {if (output>2) {print(paste ('do not combine ', 
  cluster.names[max.row],"with", cluster.names[max.col],'new beta =', round (beta.combined,2),
  'old beta1 =',round( beta1,2),"old beta2 =",round(beta2,2)))}
												test.beta <- FALSE  }},
  {if (beta.combined < max(beta1,beta2)) {if (output>2) {print(paste ('do not combine ',
  cluster.names[max.row],"with", cluster.names[max.col],'new beta =', round (beta.combined,2),
  'old beta1 =',round( beta1,2),"old beta2 =",round(beta2,2)))}
												test.beta <- FALSE  }}) #end switch  
												
  }   #end if options$beta.size

}




if(test.beta&test.alpha)   {break } else  {
if (beta.combined < ICLUST.options$beta.min) {
			keep.clustering <- FALSE       #the most similiar pair is not very similar, we should quit
			break}  else {sim[max.row,max.col] <- NA 
		sim[max.col,max.row] <- NA }
    }   #end of test.beta&test.alpha
    
}    #end of while test.alpha&test.beta.loop

#combine and summarize
if (keep.clustering) 
   {          # we have based the alpha and beta tests, now combine these two variables
clusters[,max.row] <- clusters[,max.row] + sign.max *  clusters[,max.col]  
cluster.names <- colnames(clusters)

#summarize the results
results[count,1] <- cluster.names[max.row]
results[count,2] <- cluster.names[max.col]
results[count,"similarity"] <- sim[max.cell]
results[count,"correlation"] <- sim.mat[max.cell]
results[count,"alpha1"] <- item.rel[max.row]
results[count,"alpha2"] <- item.rel[max.col]
size1 <- cluster.stats$size[max.row]
size2 <- cluster.stats$size[max.col]
results[count,"size1"] <- size1
results[count,"size2"] <- size2
results[count,"beta1"] <-  beta1
results[count,"beta2"] <-  beta2
results[count,"rbar1"] <-  rbar1
results[count,"rbar2"] <-  rbar2
results[count,"r1"] <- r1
results[count,"r2"] <- r2


results[count,"beta"] <- beta.combined

results[count,'alpha'] <- alpha
results[count,'rbar'] <- rbar
results[count,"size"] <- size12
results[count,3:18] <- round(results[count,3:18],ICLUST.options$digits)
#update
cluster.names[max.row] <- paste("C",count,sep="")
colnames(clusters) <- cluster.names
clusters <- clusters[,-max.col]
cluster.names<- colnames(clusters)

#row.max <- row.max[-max.col]    



	}  #end of combine section

if(output > 1) print(results[count,],digits=2)
count=count+1 
if ((num.var - count) < ICLUST.options$n.clus) {keep.clustering <- FALSE}
if(num.var - count < 1) {keep.clustering <- FALSE}   #only one cluster left
	
}  #end of keep clustering loop

ICLUST.cluster <- list(results=results,clusters=clusters,number <- num.var - count)
}   # end ICLUST.cluster