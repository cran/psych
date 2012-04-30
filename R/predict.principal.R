#written November 22, 2010
#revised April 7, 2012 to treat single case problem

"predict.psych" <-
function(object,data,old.data,...) {
data <- as.matrix(data)
if(ncol(data) ==1) data <- t(data)
if(missing(old.data)) {data <- scale(data)} else {
	stats <- describe(old.data)
	data <- scale(data,center=stats$mean,scale=stats$sd)}
wt <- object$weights
pred <- data %*% wt
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