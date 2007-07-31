#corrected estimate of communality, May 21, 2007
#removed "x" from factanal call, June 14, 2007
"schmid" <-
function (model, nfactors = 3, pc = "pa",...) 
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
       obminfact <- oblimin(orth.load)
    rownames(obminfact$loadings) <- attr(model,"dimnames")[[1]]
    fload <- obminfact$loadings
    factr <- t(obminfact$Th) %*% (obminfact$Th)
    gfactor <- factanal( covmat = factr, factors = 1)
    gload <- loadings(gfactor)
    gprimaryload <- fload %*% gload
    colnames(gprimaryload) <- "g factor"
    u2 <- 1 - diag(orth.load %*% t(orth.load)) 
    h2 <- 1 - u2                         
    uniq <- 1 - fload^2
    guniq <- 1 - gprimaryload^2
    Ig <- matrix(0, ncol = nfactors, nrow = nfactors)
    diag(Ig) <- gload
    primeload <- fload %*% Ig
    uniq2 <- 1 - uniq - primeload^2
    sm <- sqrt(uniq2)
    schmid <- list(sl = cbind(gprimaryload, sm, h2, u2), orthog = fact$loadings, oblique = fload, 
        fcor = factr, gloading = gload)
}
