grouplag <- function(x,gr,varnum=NULL) {
result <- by(x,gr,lag,varnum)
return(result)
}

lag <- function(x,varnum) {if(!is.null(varnum))  x <- x[-varnum]
cname <- colnames(x)
results <- cbind(x, x[row(x[,1,drop=FALSE])+1,])
colnames(results) <- c(cname,paste(cname,2,sep=""))
return(results)
}

