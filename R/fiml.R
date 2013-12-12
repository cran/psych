#December 12, 2013
#taken almost completely from the matching lavaan functions
#which unfortunately, are not public functions
#some of the functionality of lavaan has been dropped 
#The relevant functions from lavaan are
#getMissingPatterns
#getMissingPatternStats 
#estimate.moments.fiml 
#minimize.this.function
#first.derivative.param
# first.derivative.param.numerical
#estimator.FIML
#derivative.FIML
#vech

corFiml <- 
function (x, covar = FALSE,show=FALSE) 
{ if (!is.matrix(x)) 
            x <- as.matrix(x)
        Mp <- getMissingPatterns(x)
        if (length(Mp$empty.idx) > 0L) {
            x <- x[-Mp$empty.idx, , drop = FALSE]
        }
        mpat <- getMissingPatternStats(X = x, Mp = Mp)
        if(show) {return(Mp$pat) } else {
        moments <- estimate.moments.fiml(X = x, M = mpat)
        colnames(moments$sigma) <- rownames(moments$sigma) <- colnames(x)
        cor <- cov2cor(moments$sigma)
        if (covar) {
            return(list(mean = moments$mu, cor = cor, cov = moments$sigma, 
                  fx = moments$fx))
            } else {return(cor)}
}
}



getMissingPatterns <- 
function (X) 
{
    nobs <- nrow(X)
    nvar <- ncol(X)
    MISSING <- 1L * is.na(X)  #convert to number
    coverage <- crossprod(1 - MISSING)/nobs
    #this next step looks for the missing cases and removes someone with all missing
    id <- apply(MISSING, MARGIN = 1, function(x) {
        if (sum(x) == length(x)) {
            out <- "empty"
        }
        else {
            paste(x, collapse = "")
        }
    })
    empty.idx <- which(id == "empty")
    if (length(empty.idx) > 0) {
        MISSING <- MISSING[-empty.idx, ]
        X <- X[-empty.idx, ]
        id <- id[-empty.idx]
        nobs <- nobs - length(empty.idx)
    }
    
    TABLE <- sort(table(id), decreasing = TRUE)
    order <- names(TABLE)
    npatterns <- length(TABLE)
    pat <- 1L - MISSING[match(order, id), , drop = FALSE]
    storage.mode(pat) <- "logical"
    row.names(pat) <- as.character(TABLE)
    out <- list(nobs = nobs, nvar = nvar, coverage = coverage, 
        id = id, npatterns = npatterns, order = order, pat = pat, 
        empty.idx = empty.idx)
    out
}


getMissingPatternStats <- 
function (X = NULL, Mp = NULL) 
{ npatterns <- Mp$npatterns
    id <- Mp$id
    order <- Mp$order
    pat <- Mp$pat
    data <- vector("list", length = npatterns)
    for (p in 1:npatterns) {
        row.idx <- which(id == order[p])
        nobs <- length(row.idx)
        Xp <- X[row.idx, pat[p, ], drop = FALSE]
        if (nobs > 1) {
            M <- colMeans(Xp)
            S <- crossprod(Xp)/nobs - tcrossprod(M)
        }
        else {
            S <- 0
            M <- as.numeric(Xp)
        }
        data[[p]] <- list(X = Xp, SX = S, MX = M, nobs = nobs, 
            var.idx = pat[p, ])
    }
    data
}


estimate.moments.fiml <- 
function (X = NULL, M = NULL) 
{
    nvar <- ncol(X)
    pstar <- nvar * (nvar + 1)/2
    start.cov <- cov(X, use = "p")
    dimnames(start.cov) <- NULL
    start.mean <- apply(X, 2, mean, na.rm = TRUE)
    names(start.mean) <- NULL
    lower.idx <- which(lower.tri(start.cov, diag = TRUE))
    upper.idx <- which(upper.tri(t(start.cov), diag = TRUE))
    
x2param <- function(x) {
        mu <- x[1:nvar]
        sigma.el <- x[-(1:nvar)]
        sigma <- matrix(0, nvar, nvar)
        sigma[lower.idx] <- sigma.el
        sigma[upper.idx] <- t(sigma)[upper.idx]
        list(mu = mu, sigma = sigma)
    }
    
minimize.this.function <- function(x) {
        out <- x2param(x)
        ev <- eigen(out$sigma)$values
        if (any(ev < 0)) {
            return(Inf)
        }        
        fx <- estimator.FIML(Sigma.hat = out$sigma, Mu.hat = out$mu, 
            M = M)
        fx
    }
    
first.derivative.param <- function(x) {
        out <- x2param(x)
        dx.out <- derivative.FIML(Sigma.hat = out$sigma, Mu.hat = out$mu,  M = M)
        dx <- c(dx.out$dx.mu, vech(dx.out$dx.Sigma))
        dx
    }
    
start.x <- c(start.mean, vech(start.cov))
    iter.max <- 500
    optim.out <- nlminb(start = start.x, objective = minimize.this.function, 
        gradient = first.derivative.param, control = list(iter.max = iter.max, 
            eval.max = iter.max * 2, trace = 0))
   
    x <- optim.out$par
    fx <- optim.out$objective
    out <- x2param(x)
    sigma <- out$sigma
    mu <- out$mu
  
    list(sigma = sigma, mu = mu, fx = fx)
}


estimator.FIML <- 
function (Sigma.hat = NULL, Mu.hat = NULL, M = NULL, h1 = NULL) 
{
    npatterns <- length(M)
    fx.p <- numeric(npatterns)
    w.p <- numeric(npatterns)
    for (p in 1:npatterns) {
        SX <- M[[p]][["SX"]]
        MX <- M[[p]][["MX"]]
        w.p[p] <- nobs <- M[[p]][["nobs"]]
        var.idx <- M[[p]][["var.idx"]]
        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx], logdet = TRUE)
        Sigma.log.det <- attr(Sigma.inv, "logdet")
        Mu <- Mu.hat[var.idx]
        TT <- SX + tcrossprod(MX - Mu)
        trace <- sum(Sigma.inv * TT)
        fx.p[p] <- Sigma.log.det + trace
    }
    fx <- weighted.mean(fx.p, w = w.p)
    if (!is.null(h1)) {
        fx <- fx - h1
        if (fx < 0) 
            fx <- 0
    }
    fx
}

inv.chol <- 
function (S, logdet = FALSE) 
{
    cS <- chol(S)
    S.inv <- chol2inv(cS)
    if (logdet) {
        attr(S.inv, "logdet") <- sum(log(diag(cS)^2))
    }
    S.inv
}

derivative.FIML <- 
function (Sigma.hat, Mu.hat, M) 
{
    ntotal <- sum(sapply(M, "[[", "nobs"))
    nvar <- length(Mu.hat)
    npatterns <- length(M)
    dx.Sigma <- matrix(0, nvar, nvar)
    dx.Mu <- matrix(0, nvar, 1)
    for (p in 1:npatterns) {
        SX <- M[[p]][["SX"]]
        MX <- M[[p]][["MX"]]
        nobs <- M[[p]][["nobs"]]
        var.idx <- M[[p]][["var.idx"]]
        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx], logdet = FALSE)
        Mu <- Mu.hat[var.idx]
        TT <- SX + tcrossprod(MX - Mu)
        dx.Mu[var.idx, 1] <- (dx.Mu[var.idx, 1] + nobs/ntotal * 
            -2 * t(t(MX - Mu) %*% Sigma.inv))
        dx.Sigma[var.idx, var.idx] <- (dx.Sigma[var.idx, var.idx] - 
            nobs/ntotal * 2 * (Sigma.inv %*% (TT - Sigma.hat[var.idx, 
                var.idx]) %*% Sigma.inv))
    }
    diag(dx.Sigma) <- diag(dx.Sigma)/2
    out <- list(dx.mu = dx.Mu, dx.Sigma = dx.Sigma)
    out
}


vech <- 
function (S, diagonal = TRUE) 
{
    ROW <- row(S)
    COL <- col(S)
    if (diagonal) 
        S[ROW >= COL]
    else S[ROW > COL]
}



