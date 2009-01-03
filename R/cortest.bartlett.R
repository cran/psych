"cortest.bartlett" <-
function(R,n=NULL) {
#message("Bartlett's test that correlation matrix is an identity matrix")
if (dim(R)[1] != dim(R)[2]) {n <- dim(R)[1] 
                             message("R was not square, finding R from data")
                             R <- cor(R,use="pairwise")}
                             
p <- dim(R)[2]
if(!is.matrix(R) ) R <- as.matrix(R)  #converts data.frames to matrices if needed
if(is.null(n)) {n <- 100 
                warning("n not specified, 100 used") }
detR <- det(R)
statistic  <- -log(detR) *(n -1 - (2*p + 5)/6)
df <- p * (p-1)/2
pval <- pchisq(statistic,df,lower.tail=FALSE)

bartlett <- list(chisq = statistic, p.value =pval, df =df)
return(bartlett) }

