"lavaan.diagram" <- 
function(fit,model="cfa",...) {
if(model=="cfa") {structure.diagram(fx=fit@Fit@GLIST[[1]]$lambda,Phi=fit@Fit@GLIST[[1]]$psi,Rx=fit@Fit@GLIST[[1]]$theta,...)}
else {structure.diagram(fx=fit@Fit@GLIST[[1]]$lambda,Phi=fit@Fit@GLIST[[1]]$beta,Rx=fit@Fit@GLIST[[1]]$theta,...)}
}
