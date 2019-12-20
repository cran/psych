#written November 22, 2010
#revised April 7, 2012 to treat single case problem
#revised October 30, 2019 to include bestScales options
"predict.psych" <-
function(object,data,old.data,options=NULL,...) {
  obnames <- cs(fa,bestScales )
     value <- inherits(object, obnames, which=TRUE)
			   if (any(value > 1)) {value <- obnames[which(value >0)]} else {value <- "none"}

if(value !="bestScales") value <- "fa"

switch(value, 

fa = {
	data <- as.matrix(data)
	if(ncol(data) ==1) data <- t(data)
	if(missing(old.data)) {data <- scale(data)} else {
	stats <- describe(old.data)
	data <- scale(data,center=stats$mean,scale=stats$sd)}
	wt <- object$weights
	pred <- data %*% wt},


bestScales = {
if(!is.null(options)) {keys<- options} else {keys <- "best.keys"}
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
})

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