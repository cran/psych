"describe.by" <-
function (x,group,mat=FALSE,...) {               #data are x, grouping variable is group
answer <- by(x,group,describe,...)
if (mat) { ncol <- dim(answer[[1]])[2]
	n.var <- dim(answer[[1]])[1]
	if (n.var ==1) {
	names  <- colnames(answer[[1]])
	row.names <-names(answer)
	answer <- matrix(unlist(answer),ncol=ncol,byrow=TRUE) 
	colnames(answer) <- names
	rownames(answer) <- row.names } else {stop("matrix output not yet implemented for more than one variable ") 
	                  mat.ans <- matrix(NA,ncol=ncol,nrow=length(answer))
}  } 
return(answer)}