#written November 22, 2010

"predict.psych" <-
function(object,data,old.data,...) {
if(missing(old.data)) {data <- scale(data)} else {
	stats <- describe(old.data)
	data <- scale(data,center=stats$mean,scale=stats$sd)}
wt <- object$weights
pred <- data %*% wt
return(pred)}




"predict.principal" <-
function(object,data) {
pc <- object$weights
data <- as.matrix(data)
pred <- data %*% pc
return(pred)
}

"predict.fa" <-
function(object,data) {
wt <- object$weights
data <- as.matrix(data)
pred <- data %*% wt
return(pred)
}