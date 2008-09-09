"VSS.scree" <-
function(rx,main="scree plot") {
 nvar <- dim(rx)[2]
 if (nvar != dim(rx)[1]) {rx <- cor(rx,use="pairwise")}
 values <-  eigen(rx)$values
plot(values,type="b", main = main )
abline(h=1) }
