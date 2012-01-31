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
"plot.residuals" <-
function(x,main="QQ plot of residuals",qq=TRUE,...) {
if(qq) {
x <- x[lower.tri(x)]
std.x <- x/sd(x)
qqnorm(std.x,main=main)
qqline(std.x)
} else {
cor.plot(x,main=main)}}