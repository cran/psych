"VSS.scree" <-
function(rx,main="scree plot") {
 nvar <- dim(rx)[2]
 if (nvar != dim(rx)[1]) {rx <- cor(rx,use="pairwise")}
 values <-  eigen(rx)$values
plot(values,type="b", main = main,ylab="Eigen values of components",xlab=" component number")
abline(h=1) }



"scree" <-
function(rx,factors=TRUE,pc=TRUE,main="Scree plot",hline=NULL,add=FALSE) {
cl <- match.call()
 nvar <- dim(rx)[2]
 if (nvar != dim(rx)[1]) {rx <- cor(rx,use="pairwise")}
if(pc) {values <-  eigen(rx)$values
        if(factors) {ylab="Eigen values of factors and components"
         xlab="factor or component number"} else {ylab="Eigen values of components"
         xlab=" component number"}
        } else {values <- fa(rx)$values
                ylab="Eigen values of factors"
                xlab=" factor number"
                                           factors <- FALSE }
max.value <- max(values)

if(!add) {plot(values,type="b", main = main ,pch=16,ylim=c(0,max.value),ylab=ylab,xlab=xlab)} else {
          points(values,type="b",  ,pch=16)  }
          
if(factors) {
   fv <- fa(rx)$values
   points(fv,type="b",pch=21,lty="dotted")
   } else {fv <- NULL}
   if(is.null(hline)) {abline(h=1)} else {abline(h=hline) }

if(factors & pc) { legend("topright", c("PC ","FA"),pch=c(16,21),
       text.col = "green4", lty = c("solid","dotted"),
       merge = TRUE, bg = 'gray90')}
if(pc) {
     results <- list(fv = fv, pcv=values,call=cl)} else {
     results <- list(fv=values, pcv=NULL,call=cl) }
class(results) <- c("psych","scree")
 invisible(results)
       }
