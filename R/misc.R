#A number of useful helper functions
#added January, 2012

#a dummy function to allow the help to find misc
"psych.misc" <-
function() {}

"lowerMat" <- 
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
	
"lowerCor" <- 
function(x,digits=2,use="pairwise",method="pearson") {
   R <- cor(x,use=use,method=method)
   lowerMat(R,digits)
   invisible(R)
   }
	
	
#adapted from utils::txtProgressBar
#modified August 10, 2012 to print just 100 times. 
"progressBar" <- 
function(value,max,label=NULL) {
if(class(stdout())[1]=="terminal")  #only print to the screen, not to a file
 { pc <- round(100 * value/max)
if(ceiling(100 * value/max)==floor(100 * value/max)) {
width <- 80
char="."
nw <- nchar(char, "w")
nb <- round(width * value/max )
       
        cat(paste(c("\r  ",label," |", rep.int(" ", nw * width + 6)), collapse = ""))
        cat(paste(c("\r  ",label," |", rep.int(char, nb), rep.int(" ", 
            nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""))}
            }
}

