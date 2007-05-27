"poly.mat" <- 
function(x,short=TRUE,std.err=FALSE,ML=FALSE) {
require(polycor)  #John Fox's Polycor package
xm <- as.matrix(x)   
xm.cat <- matrix(as.factor(xm),ncol=dim(xm)[2])
colnames(xm.cat) <- colnames(xm)
r.het <- hetcor(xm.cat,std.err=std.err,ML=ML)
rownames(r.het$correlations) <- colnames(r.het$correlations) <- colnames(xm)
if(short) {return(r.het$correlations)} else {return(r.het)}
}