"residuals.psych" <- 
function(object,...) {
result <- NULL
if(length(class(object)) > 1)  { value <- class(object)[2] } else {stop("No appropriate residual found")}
switch(value,

fa = {residual <- object$residual},
principal = {residual <- object$residual},
omega = {residual <- object$stats$residual},
irt.fa ={residual <- object$fa$residual},
extension = {residual <- object$resid})

class(residual) <- c("psych","residuals")
return(residual)
}
#added January 30, 2012
"resid.psych" <- 
function(object,...) {
    residuals(object)
}




#added Feb 4, 2012
#modified April 15, 2016 to add chisquare and histograms as well as to identify variables
"plot.residuals" <-
function(x,main,type=c("qq","chi","hist","cor"),std,bad=4,numbers=TRUE,upper=FALSE,diag=FALSE,...) {
if(missing(type)) type <- "qq"
nr <- nrow(x)
nc <- ncol(x)
if(!is.null(rownames(x))) {rname <- rownames(x)} else {rname <- paste0("V",1:nr)}
diag(x) <- NA

switch(type,
hist = {
if(missing(std)) std <- FALSE
x <- x[lower.tri(x,diag=TRUE)]
std.x <- x/sd(x,na.rm=TRUE)
if(std) {if(missing(main)) main <- "Histogram of standardized residuals"
   hist(std.x,main=main,...)} else {
if(missing(main)) main <- "Histogram of residuals"
hist(x,main=main,...)}},

qq= {  if(missing(std)) std <- TRUE
		x <- x[lower.tri(x,diag=TRUE)]
		if(std) {
     	if(missing(main)) main <- "Plot of standardized residuals"
		std.x <- x/sd(x,na.rm=TRUE)
		xy <- qqnorm(std.x,main=main)
		qqline(std.x)
		 worst <- order(abs(std.x), decreasing = TRUE)
		} else {
		  if(missing(main)) main <- "Plot of raw residuals"
		xy <- qqnorm(x,main=main,...)
		qqline(x)
		worst <- order(abs(x), decreasing = TRUE)}

    worstItems <- arrayInd(worst[1:bad],c(nr,nc))
    text(xy$x[worst[1:bad]],xy$y[worst[1:bad]],paste(rname[worstItems[,2]],rname[worstItems[,1]]),pos=4,...)
},

chi = {#note that xy reported for qqplot is already sorted
 if(missing(std)) std <- TRUE
   x <- x[lower.tri(x,diag=TRUE)] 
  if(std) {x <- x/sd(x,na.rm=TRUE)
           if(missing(main)) main <- "Plot of squared standardized residuals"} else {
           if(missing(main)) main <- "Plot of squared residuals"}
  nx <- length(x) - nr
xy <-   qqplot(qchisq(ppoints(nx),df=1),y=x^2,main=main,ylab="Quantiles of Squared residuals",xlab="Expected value for quantile")
  qqline(x^2,distribution=function(p) qchisq(p,df=1))
   worst <- order(abs(x^2), decreasing = TRUE)
 worstItems <- arrayInd(worst[1:5],c(nr,nc))
 text(xy$x[nx:(nx-4)],xy$y[nx:(nx-4)],paste(rname[worstItems[,2]],rname[worstItems[,1]]),pos=2,...)
  },

cor= {if(missing(main)) main <- "Plot of residual correlations"
cor.plot(x,main=main,numbers=numbers,upper=upper,diag=diag)})

}

