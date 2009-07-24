#polar   a function to convert a factor loadings matrix to polar coordinates
#Version of Sept 9, 2007
#Slightly revised July 25, 2009 
"polar" <- 
function(f,sort=TRUE)  #return all pairwise polar coordinates
	{ if (!is.matrix(f) && !is.data.frame(f) ) {fload <-f$loadings} else {fload <- f}
	
	nf <- dim(fload)[2]
	n.var <- dim(fload)[1]
	polar  <- matrix(0,nrow=n.var,ncol = nf * (nf-1)+1)
	if (!is.null(rownames(fload))) {rownames(polar) <- rownames(fload)} else {rownames(polar) <- paste("v",1:n.var,sep="") }
	colnames(polar) <- rep(NA,nf*(nf-1)+1)  #just give it something to play with
	
	polar[,1] <- seq(1:n.var)
	colnames(polar)[1] <- "Var"
	k <- 2
	kk <- nf*(nf-1)/2
	for (i in 2:nf) {
	   for (j in 1:(i-1)) {
	  	 	vector.length <- fload[,i]^2 + fload[,j]^2
			theta=sign(fload[,i])*180*acos(fload[,j]/sqrt(vector.length))/pi #vector angle (-180: 180) 
			polar[,k] <- theta %% 360
			polar[,k+kk] <- vector.length
			colnames(polar)[k] <- paste("theta",i,j,sep="")
			colnames(polar)[k+kk] <- paste("vecl",i ,j,sep="")
			k <- k + 1
			}
		}
	
	if (sort) {polar <- polar[order(polar[,2]),] }
	return(polar) 
	} 
	