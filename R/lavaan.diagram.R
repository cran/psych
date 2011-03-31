"lavaan.diagram" <- 
function(fit,title,...) {
if (is.null(fit@Model@GLIST$beta)) {model <- "cfa"} else {model <- "sem"}
if(missing(title)) {if(model =="cfa") { title="Confirmatory structure" } else {title = "Structural model"}
                       }
        fx=fit@Model@GLIST$lambda
        colnames(fx) <- fit@Model@dimNames$lambda[[2]]
         Phi <-  fit@Model@GLIST$psi
                  Rx <- fit@Model@GLIST$theta
                  v.labels <- fit@Model@dimNames$lambda[[1]]
  if(model=="cfa") {                
structure.diagram(fx=fx,Phi=Phi,Rx=Rx,labels=v.labels,main=title,...)}
else {structure.diagram(fx=fx,Phi=fit@Model@GLIST$beta,Rx=Rx,labels=v.labels,main=title,...) }
}