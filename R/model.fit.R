model.fit <- function(R,model,N, digits=5, conf.level=.90, ...) {
#summary.sem adapted to factanal or factor.pa output
   # norm.res <- normalized.residuals(object)
  #  se <- sqrt(diag(object$cov))
  se <- sqrt(diag(R))
  z <- R$coeff/se       
  #  z <- object$coeff/se
  #  n.fix <- object$n.fix
  #still under development
  object <- list()
  object$raw <- FALSE
    n <- dim(R)[2]
    n.fix <- 0
    t <- 0
    S <- R
    C <- model %*% t(model)
    diag(C) <- 1
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    dfNull <- n*(n - 1)/2
    invC <- solve(C)
    CSC <- invC %*% (S - C)
    CSC <- CSC %*% CSC
    CS <- invC %*% S
    CS <- CS %*% CS
    chisq <- object$criterion * (N - (!object$raw))
    chisqNull <- object$chisqNull
    GFI <- if (!object$raw) 1 - sum(diag(CSC))/sum(diag(CS)) else NA
    if ((!object$raw) && df > 0){
        AGFI <- 1 - (n*(n + 1)/(2*df))*(1 - GFI)
        NFI <- (chisqNull - chisq)/chisqNull
        NNFI <- (chisqNull/dfNull - chisq/df)/(chisqNull/dfNull - 1)
        L1 <- max(chisq - df, 0)
        L0 <- max(L1, chisqNull - dfNull)
        CFI <- 1 - L1/L0
        RMSEA <- sqrt(max(object$criterion/df - 1/(N - (!object$raw)), 0))
        tail <- (1 - conf.level)/2 
        max <- 0
        while (max > 1){
            res <- optimize(function(lam) (tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
            if (sqrt(res$objective) < tail/100) break
            max <- max/2
            }
        lam.U <- if (max <= 1) NA else res$minimum
        max <- max(max, 1)
        while (max > 1){
            res <- optimize(function(lam) (1 - tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
            if (sqrt(res$objective) < tail/100) break
            max <- max/2
            }
        lam.L <- if (max <= 1) NA else res$minimum
        RMSEA.U <- sqrt(lam.U/((N - (!object$raw))*df))
        RMSEA.L <- sqrt(lam.L/((N - (!object$raw))*df))
        }
    else RMSEA.U <- RMSEA.L <- RMSEA <- NFI <- NNFI <- CFI <- AGFI <- NA
    RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
    if (!is.null(object$coeff)){
        var.names <- rownames(object$A)
        ram <- object$ram[object$par.posn, , drop=FALSE]
        par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
                        var.names[ram[,3]])
        coeff <- data.frame(object$coeff, se, z, 2*(1 - pnorm(abs(z))), par.code)
        names(coeff) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)", " ")
        row.names(coeff) <- names(object$coeff)
        }
    else coeff <- NULL
    BIC <- chisq - df * log(N)
    SRMR <- sqrt(sum(standardized.residuals(object)^2 * lower.tri(diag(n), diag=TRUE))/(n*(n + 1)/2))
 
   norm.res <- NA
    ans <- list(chisq=chisq, df=df, chisqNull=chisqNull, dfNull=dfNull,
        GFI=GFI, AGFI=AGFI, RMSEA=RMSEA, NFI=NFI, NNFI=NNFI, CFI=CFI, BIC=BIC, SRMR=SRMR, 
        norm.res=norm.res, coeff=coeff, digits=digits, 
        iterations=object$iterations, aliased=object$aliased, raw=object$raw)
    class(ans) <- "summary.sem"
    ans
    }
    