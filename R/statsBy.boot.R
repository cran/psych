statsBy.boot <-    function (data,group,ntrials=10,cors=FALSE,replace=TRUE,method="pearson") { #  
    cl <- match.call()
    result <- vector("list",ntrials) #supposedly allocates more memory
    for (i in 1:ntrials) {
    progressBar(i,ntrials,"statsBy")
    data[,group] <- sample(data[,group],size=nrow(data),replace=replace)
    result[[i]] <- statsBy(data,group,cors=cors,method=method)
    }
return(result)
}
   
statsBy.boot.summary <- function(res.list,var="ICC2") {
nreps <- length(res.list)
nvar <- length(res.list[[1]][[var]])
cnames <- names(res.list[[1]][[var]])
temp <- matrix(NaN,ncol=nvar,nrow=nreps)
colnames(temp) <- cnames
for(i in 1:nreps){
temp[i,] <- res.list[[i]][[var]]
}
return(temp)
}
 