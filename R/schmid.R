"schmid" <-
function (model, nfactors = 3, pc = "pa",...) 
{
    if (pc=="pc") {
        fact <- principal(model, nfactors,...)
    } else {if (pc=="pa") {fact <- factor.pa(model, nfactors,...) } else {
        fact <- factanal(x, covmat = model, factors = nfactors,...)
    }}
       obminfact <- oblimin(loadings(fact))
    rownames(obminfact$loadings) <- attr(model,"dimnames")[[1]]
    fload <- obminfact$loadings
    factr <- t(obminfact$Th) %*% (obminfact$Th)
    gfactor <- factanal(x, covmat = factr, factors = 1)
    gload <- loadings(gfactor)
    gprimaryload <- fload %*% gload
    colnames(gprimaryload) <- "g factor"
    u2 <- 1 - diag(fload %*% t(fload))
    h2 <- 1 - u2
    uniq <- 1 - fload^2
    guniq <- 1 - gprimaryload^2
    Ig <- matrix(0, ncol = nfactors, nrow = nfactors)
    diag(Ig) <- gload
    primeload <- fload %*% Ig
    uniq2 <- 1 - uniq - primeload^2
    sm <- sqrt(uniq2)
    schmid <- list(sl = cbind(gprimaryload, sm, h2, u2), orthog = fload, 
        fcor = factr, gloading = gload)
}

