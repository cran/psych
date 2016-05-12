"block.random" <- 
function(n,ncond=NULL) {
if(is.null(ncond)) {ncond <- 2
                   IVs <- 1
                   conditions <- c(ncond)
                   } else {
                        if (length(ncond) < 2) {
                         IVs<- 1
                         conditions <- c(ncond)} else {
                   		  IVs <- length(ncond)
                  	 	  conditions <- ncond
                  	 	  ncond <- prod(ncond)}
                  	      }
if(floor(n/ncond) * ncond != n) {stop("number of subjects much be a multiple of number of conditions")}
blocks <- matrix(rep(NA,n*(1+IVs)),ncol=1+IVs)
rcond <- rep(NA,n)
if(is.null(names(conditions))) {colnames(blocks) <- c("blocks",paste("IV",1:IVs,sep=""))}  else {colnames(blocks) <- c("blocks",names(conditions))}
rownames(blocks) <- paste("S",1:n,sep="")
for (block in 1:(n/ncond)) {
	blocks[((block-1)*ncond+1):(block*ncond),1] <- block
	rcond [((block-1)*ncond+1):(block*ncond)] <- sample(ncond,replace=FALSE)
}
for (i in 1:IVs) {if(i<2) { blocks[,i+1]<-  ceiling( (rcond %% conditions[i] + 1))} else { blocks[,i+1]<-  ceiling( (rcond %% prod(conditions[1:i])  +1 ) /prod(conditions[1:(i-1)]))
}}
return(blocks)}