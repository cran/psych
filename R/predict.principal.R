#written November 22, 2010
#revised April 7, 2012 to treat single case problem
#revised October 30, 2019 to include bestScales options
"predict.psych" <-
function(object,data,old.data,options=NULL,missing=FALSE,impute="none",...) {
  obnames <- cs(fa,bestScales,setCor,pca, principal )
     value <- inherits(object, obnames, which=TRUE)
			   if (any(value > 1)) {value <- obnames[which(value >0)]} else {value <- "none"}

if(value %in% cs(factor,pca,principal,omega)) value <- "fa"

switch(value, 

fa = {
	data <- as.matrix(data)
	if(ncol(data) ==1) data <- t(data)
	if(missing(old.data)) {data <- scale(data)} else {
	stats <- describe(old.data)
	data <- scale(data,center=stats$mean,scale=stats$sd)}
	wt <- object$weights
	if(impute !="none") data <- impute.na(data,impute)
	if(missing) {pred <- matrixMult.na(data,wt)} else {
	pred <- data %*% wt}
	},


bestScales = {
if(!is.null(options)) {keys<- options} else {keys <- "best.keys"}
 if(impute != "none") {
 #for speed we want to just impute those items that will be scored
 select <- selectFromKeys(keys)
 data <- impute.na(data,impute)}
switch(keys,
 best.keys = {keys <- object$best.keys
              scores <- scoreVeryFast(keys,data)},
 weights =  {keys <- object$weights
         scores <- scoreWtd(keys,data)},
 optimal.keys ={ keys <- object$optimal.keys
        scores <- scoreVeryFast(keys,data)},
 optimal.weights ={ keys <- object$optimal.weights
       scores <- scoreWtd(keys,data)}
       )

criteria <- data[names(keys)]
bwt <- object$final.stats$r *  object$final.stats$crit.sd/ object$final.stats$sd
xmean <- object$final.stats$mean
ymean <- object$final.stats$crit.mean
pred <-  t(bwt *(t(scores) - xmean) + ymean )	
},

#added January 5, 2020
setCor = {
   data <- as.matrix(data)
   if(ncol(data) ==1) data <- t(data)
   vars <- rownames(object$coefficients)
   vars <- vars[ vars %in% colnames(data)]
   data <- data[,vars,drop=FALSE]  
	if(missing(old.data)) {data <- scale(data)} else {
	stats <- describe(old.data)
	data <- scale(data,center=stats$mean,scale=stats$sd)}
	wt <- object$coefficients[vars,]  #don't use the intercept
	if(impute !="none") data <- impute.na(data,impute)
	if(missing) {pred <- matrixMult.na(data,wt)} else {
	pred <- data %*% wt}
   }
)

return(pred)}



#these next two do not standardize the prediction
"predict.principal" <-
function(object,data) {
wt <- object$weights
data <- as.matrix(data)
pred <- data %*% wt
return(pred)
}

"predict.fa" <-
function(object,data) {
wt <- object$weights
data <- as.matrix(data)
pred <- data %*% wt
return(pred)
}

predict.setCor.best <- function(object,data,p=.01){
 data <- as.matrix(data)
 vars <- rownames(object$coefficients)
 vars <- vars[ vars %in% colnames(data)]
 data <- data[,vars,drop=FALSE]  
wt <- object$coefficients
prob <- object$Probability
wt <- wt* (prob < p)
wt <- wt[vars,]
pred <- data %*% wt
}

#do matrix multiplication with missing data
matrixMult.na <- function(x,y,scale=TRUE) {
nvar <- ncol(x)
if(nvar != nrow(y) ) stop("matrices are not compatible")#matrices are not compatible

if(scale) x <- scale(x) #zero center and standaridize
tx <- t(x) #we want to do  it on the transposed matrix
ny <- ncol(y)
result <- matrix(NA,nrow = nrow(x),ncol= ncol(y))
result <- apply(y,2,function(x ) colMeans(x * tx,na.rm=TRUE))
return((result))
}

"impute.na" <- function(x,impute="mean") {
  miss <- which(is.na(x),arr.ind=TRUE)
   if(impute=="mean") {
       		item.means <- colMeans(x,na.rm=TRUE)   #replace missing values with means
       		x[miss]<- item.means[miss[,2]]} else { 
       		item.med   <- apply(x,2,median,na.rm=TRUE) #replace missing with medians
        	x[miss]<- item.med[miss[,2]]} 
    return(x)
 }
 
 
 cor.na <- function(x,y=NULL,scale=TRUE) {
nvar <- ncol(x)
if(scale) {sx <- scale(x)} else {sx <- x}
if(is.null(y)) y <- sx
if(nvar != ncol(y) ) stop("matrices are not compatible")#matrices are not compatible

#if(scale) x <- scale(x) #zero center and standaridize
 #we want to do  it on the transposed matrix
ny <- ncol(y)
result <- matrix(NA,nrow = ncol(x),ncol= ncol(y))
result <- apply(y,2,function(x ) colMeans(x * sx,na.rm=TRUE))
return((result))
}