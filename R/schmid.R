#corrected estimate of communality, May 21, 2007
#removed "x" from factanal call, June 14, 2007
#added ability to do 2 factors by treating them with equal loadings Jan 2008
#added use of simplimax rotation June 2008
"schmid" <-
function (model, nfactors = 3, pc = "pa",  digits=2,rotate="oblimin",...) 
{
 #model is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      if(!require(GPArotation)) {stop("I am sorry, you need to have the  GPArotation package installed")}
      nvar <-dim(model)[2]
      if(dim(model)[1] != dim(model)[2]) model <- cor(model,use="pairwise")
     if (pc=="pc") {
        fact <- principal(model, nfactors,...)
    } else {if (pc=="pa") {fact <- factor.pa(model, nfactors,...) } else {
        fact <- factanal(covmat = model, factors = nfactors)
    }}
    orth.load <- loadings(fact)
    colnames(orth.load)  <- paste("F",1:nfactors,sep="")
      if (rotate == "simplimax") {obminfact <- simplimax(orth.load)} else {obminfact <- oblimin(orth.load)}
    rownames(obminfact$loadings) <- attr(model,"dimnames")[[1]]
    fload <- obminfact$loadings
    #factr <- t(obminfact$Th) %*% (obminfact$Th)
    factr <- obminfact$Phi
    colnames(factr) <- rownames(factr) <- paste("F",1:nfactors,sep="")
   if (nfactors>2) {
       gfactor <- factanal( covmat = factr, factors = 1)
       gload <- loadings(gfactor) } else {gload<- c(NA,NA)
       gload[1] <- sqrt(abs(factr[1,2]))
       gload[2] <- sign(factr[1,2])*sqrt(abs(factr[1,2]))
       
              warning("Three factors are required for identification -- general factor loadings set to be equal. Proceed with caution.")}  

    gprimaryload <- fload %*% gload
    colnames(gprimaryload) <- "g"
    u2 <- 1 - diag(orth.load %*% t(orth.load)) 
    h2 <- 1 - u2                         
    uniq <- 1 - fload^2
    guniq <- 1 - gprimaryload^2
    Ig <- matrix(0, ncol = nfactors, nrow = nfactors)
    diag(Ig) <- gload
    primeload <- fload %*% Ig

    uniq2 <- 1 - uniq - primeload^2
    sm <- sqrt(uniq2)
    colnames(sm) <- paste("F",1:nfactors,"*",sep="")
    if (!is.null(digits)) {schmid <- list(sl = cbind(round(gprimaryload,digits), round(sm,digits),h2=round( h2,digits), u2=round(u2,digits)), orthog = round(orth.load,digits), oblique = round(fload,digits),
        phi = round(factr,digits), gloading = round(gload,digits) )} else {
    schmid <- list(sl = cbind(gprimaryload, sm, h2, u2), orthog = fact$loadings, oblique = fload, phi = factr, gloading = gload)}
}
