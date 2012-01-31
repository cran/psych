#A number of useful helper functions
#added January, 2012

#a dummy function to allow the help to find misc
"psych.misc" <-
function() {}

"lower.mat" <- 
function(R,digits=2) {
   lowleft <- lower.tri(R,diag=TRUE)
   nvar <- ncol(R)
	nc <- digits+3
	width <- getOption("width")
	k1 <- width/(nc+2)
   if(is.null(colnames(R))) {colnames(R) <- paste("C",1:nvar,sep="")}
   if(is.null(rownames(R))) {rownames(R) <- paste("R",1:nvar,sep="")}
	colnames(R) <- abbreviate(colnames(R),minlength=digits+3)
	
	nvar <- ncol(R)
	nc <- digits+3
	#k1 <- width/(nc+2)
	
	if(k1 * nvar < width) {k1 <- nvar}  
	k1 <- floor(k1)
	fx <- format(round(R,digits=digits))
	 if(nrow(R) == ncol(R) ) {fx[!lowleft] <- ""}
	 for(k in seq(0,nvar,k1)) { if(k<nvar) {
	print(fx[(k+1):nvar,(k+1):min((k1+k),nvar)],quote=FALSE)}
	}}
	
	
	
#adapted from utils::txtProgressBar 
"progressBar" <- 
function(value,max,label=NULL) {
width <- 80
char="."
nw <- nchar(char, "w")
nb <- round(width * value/max )
        pc <- round(100 * value/max)
        cat(paste(c("\r  ",label," |", rep.int(" ", nw * width + 6)), collapse = ""))
        cat(paste(c("\r  ",label," |", rep.int(char, nb), rep.int(" ", 
            nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""))
}

