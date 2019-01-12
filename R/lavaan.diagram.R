"lavaan.diagram" <- 
function(fit,main,e.size=.1,...) {
if (is.null(fit@Model@GLIST$beta)) {model <- "cfa"} else {model <- "sem"}
if(missing(main)) {if(model =="cfa") { main="Confirmatory structure" } else {main = "Structural model"}
                       }
        mimic <- fit@Model@fixed.x
        
       if(!mimic) {#the normal case, either a cfa or a sem
        fx=fit@Model@GLIST$lambda
       # colnames(fx) <- fit@Model@dimNames$lambda[[2]]
       colnames(fx) <- fit@Model@dimNames[[1]][[2]]
       rownames(fx) <- fit@Model@dimNames[[1]][[1]]
        if(model=="sem") { fit@Model@GLIST$beta} else {  Phi <-  fit@Model@GLIST$psi}
                  Rx <- fit@Model@GLIST$theta
                  v.labels <-fit@Model@dimNames[[1]][[1]]
        } else {
         #mimic
        y.vars <- fit@Model@x.user.idx[[1]]

        nx  <- fit@Model@ov.x.dummy.lv.idx[[1]]
        fy <- as.matrix(fit@Model@GLIST$lambda[y.vars,-nx],drop=FALSE)
        fx <- as.matrix(fit@Model@GLIST$beta[-nx,nx],drop=FALSE)
        colnames(fy) <- fit@Model@dimNames[[1]][[2]][-nx]
       rownames(fy) <- fit@Model@dimNames[[1]][[1]][y.vars]
       rownames(fx) <- fit@Model@dimNames[[1]][[2]][-nx]
       colnames(fx) <- fit@Model@dimNames[[1]][[1]][-y.vars]
        # v.labels <-fit@Model@dimNames[[1]][[1]]
        v.labels <- c(rownames(fx),rownames(fy))
         Rx <- fit@Model@GLIST$theta
         Phi <- fit@Model@GLIST$beta
        }

  if(model=="cfa") {                
structure.diagram(fx=fx,Phi=Phi,Rx=Rx,labels=v.labels,main=main,e.size=e.size,...)}
else {
if(mimic) {
fx <- t(fx)
structure.diagram(fx=fx,fy=fy,Rx=Rx,main=main,e.size=e.size,...)} else {
#a cfa model
#Phi <-  t(fit@Model@GLIST$beta) 
Phi <-  (fit@Model@GLIST$beta) 

structure.diagram(fx=fx,Phi=Phi,Rx=Rx,labels=v.labels,main=main,e.size=e.size,...) }
}
}
#modified 11/6/14 to draw the regression paths
#modified 11/14/18 to properly do mimic models

#created August 17, 2017 to allow sem.diagrams and graphs from sem output 
"sem.diagram" <- function(fit,main="A SEM from the sem package",...) {
  nvar <- ncol(fit$S)
  var.names <- fit$var.names
  tot.var <- ncol(fit$A)
  num.factors <- length(var.names) - nvar
  fx <- fit$A[1:nvar,(nvar+1):ncol(fit$A)]
  Phi <- fit$P[(nvar+1):tot.var,(nvar+1):tot.var]
 structure.diagram(fx,Phi,main=main,...)
 }
 

"sem.graph" <- function(fit,out.file=NULL,main="A SEM from the sem package",...) {
  nvar <- ncol(fit$S)
  var.names <- fit$var.names
  tot.var <- ncol(fit$A)
  num.factors <- length(var.names) - nvar
  fx <- fit$A[1:nvar,(nvar+1):ncol(fit$A)]
  Phi <- fit$P[(nvar+1):tot.var,(nvar+1):tot.var]
 structure.graph(fx,Phi,out.file=out.file,title=main,...)
 }
  