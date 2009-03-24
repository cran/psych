#modified March 4, 2009 for matrix output
"describe.by" <-
function (x,group=NULL,mat=FALSE,...) {               #data are x, grouping variable is group
if(is.null(group)) {answer <- describe(x)
   warning("no grouping variable requested")} else {
answer <- by(x,group,describe,...)
       class(answer) <- c("psych","describe")}

if (mat) { ncol <- length(answer[[1]])  #the more complicated case. How to reorder a list of data.frames
	n.var <- length(answer[[1]][[2]])
	n.groups <- dim(answer)
	
	names  <- names(answer[[1]])
	row.names <-attr(answer[[1]],"row.names")
  
	 mat.ans <- matrix(NA,ncol=(ncol+1),nrow=n.var*n.groups)
	 colnames(mat.ans) <- c("Group",names)
	 rn <- 1:n.var*n.groups
	 k <- 1
	
	 for (var  in 1:n.var) {
	  for (group in 1:n.groups) {
	  rn[k] <- paste(row.names[var],group,sep="")
	    mat.ans[k,1] <- group
	    for (stat in 1:ncol) {
	    mat.ans[k,stat+1] <- answer[[group]][[stat]][var] } 
	     k <- k+ 1}
	   }
	  answer <- data.frame( mat.ans) 
	  rownames(answer) <- rn
	   class(answer) <- c("psych","describe","list")}

return(answer)}