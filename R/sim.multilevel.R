sim.multilevel <- function(nvar=9,ngroups=4,ncases=16,rwg,rbg,eta) {
 e.wg <- eigen(rwg)
 v.wg <- pmax(e.wg$values,0)
 etabg <- sqrt(1-eta^2)
 e.bg <- eigen(rbg)
 v.bg <- pmax(e.bg$values,0)
 wg<- matrix(rnorm(nvar*ncases),ncases)
 wg <- scale(wg)
 wg <- t(e.wg$vectors %*% sqrt(diag(v.wg)) %*% t(wg))
 bg <- matrix(rnorm(nvar*ngroups),ngroups)
 bg <- scale(bg)
 bg <- e.bg$vectors %*% sqrt(diag(v.bg)) %*% t(bg)
 bg <- matrix(rep(bg, (ncases/ngroups)),nrow=ncases,byrow=TRUE)
 gr <- rep((1:ngroups),(ncases/ngroups))
 XY <- wg %*% diag(eta^2)   +  bg %*% diag(etabg^2) 
 XY <- cbind(gr,XY)
 colnames(XY) <- c("Group",paste("V",1:nvar,sep=""))
 result <- list(wg=wg,bg=bg,xy=XY)
 }
 
 