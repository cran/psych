"VSS.scree" <-
function(rx,main="scree plot") {
 nvar <- dim(rx)[2]
 if (nvar != dim(rx)[1]) {rx <- cor(rx,use="pairwise")}
 values <-  eigen(rx)$values
plot(values,type="b", main = main )
abline(h=1) }



"scree" <-
function(rx,factors=TRUE,main="Scree plot") {
 nvar <- dim(rx)[2]
 if (nvar != dim(rx)[1]) {rx <- cor(rx,use="pairwise")}
 values <-  eigen(rx)$values
plot(values,type="b", main = main ,pch=16,ylab="Eigen values of factors and components",xlab="factor or component number")
if(factors) {
   fv <- fa(rx)$values
   points(fv,type="b",pch=21,lty="dotted")}
abline(h=1)

 legend("topright", c("PC ","FA"),pch=c(16,21),
       text.col = "green4", lty = c("solid","dotted"),
       merge = TRUE, bg = 'gray90')
       }
