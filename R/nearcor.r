# Copyright (2007) Jens Oehlschlägel
# GPL licence, no warranty, use at your own risk

#! \name{nearcor}
#! \alias{nearcor}
#! \title{ function to find the nearest proper correlation matrix given an improper one }
#! \description{
#!   This function smooths a improper correlation matrix as it can result from \code{\link{cor}} with \code{use="pairwise.complete.obs"} or \code{\link[polycor]{hetcor}}.
#! }
#! \usage{
#! nearcor(R, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, maxits = 100, verbose = FALSE)
#! }
#! \arguments{
#!   \item{R}{ a square symmetric approximate correlation matrix }
#!   \item{eig.tol}{ defines relative positiveness of eigenvalues compared to largest, default=1.0e-6 }
#!   \item{conv.tol}{ convergence tolerance for algorithm, default=1.0e-7  }
#!   \item{posd.tol}{ tolerance for enforcing positive definiteness, default=1.0e-8 }
#!   \item{maxits}{ maximum number of iterations allowed }
#!   \item{verbose}{ set to TRUE to verbose convergence }
#! }
#! \details{
#!   This implements the algorithm of Higham (2002), then forces symmetry, then forces positive definiteness using code from \code{\link[sfsmisc]{posdefify}}.
#!   This implementation does not make use of direct LAPACK access for tuning purposes as in the MATLAB code of Lucas (2001).
#!   The algorithm of Knol DL and ten Berge (1989) (not implemented here) is more general in (1) that it allows contraints to fix some rows (and columns) of the matrix and (2) to force the smallest eigenvalue to have a certain value.
#! }
#! \value{
#!   A LIST, with components
#!   \item{cor}{resulting correlation matrix}
#!   \item{fnorm}{Froebenius norm of difference of input and output}
#!   \item{iterations}{number of iterations used}
#!   \item{converged}{logical}
#! }
#! \references{
#!        Knol, DL and ten Berge, JMF (1989). Least-squares approximation of an improper correlation matrix by a proper one.  Psychometrika, 54, 53-61.
#!   \cr  Higham (2002). Computing the nearest correlation matrix - a problem from finance, IMA Journal of Numerical Analysis, 22, 329-343.
#!   \cr  Lucas (2001). Computing nearest covariance and correlation matrices. A thesis submitted to the University of Manchester for the degree of Master of Science in the Faculty of Science and Engeneering.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[polycor]{hetcor}}, \code{\link{eigen}}, \code{\link[sfsmisc]{posdefify}} }
#! \examples{
#!   cat("pr is the example matrix used in Knol DL, ten Berge (1989)\n")
#!   pr <- structure(c(1, 0.477, 0.644, 0.478, 0.651, 0.826, 0.477, 1, 0.516,
#!   0.233, 0.682, 0.75, 0.644, 0.516, 1, 0.599, 0.581, 0.742, 0.478,
#!   0.233, 0.599, 1, 0.741, 0.8, 0.651, 0.682, 0.581, 0.741, 1, 0.798,
#!   0.826, 0.75, 0.742, 0.8, 0.798, 1), .Dim = c(6, 6))
#!
#!   nr <- nearcor(pr)$cor
#!   plot(pr[lower.tri(pr)],nr[lower.tri(nr)])
#!   round(cbind(eigen(pr)$values, eigen(nr)$values), 8)
#!
#!   cat("The following will fail:\n")
#!   try(factanal(cov=pr, factors=2))
#!   cat("and this should work\n")
#!   try(factanal(cov=nr, factors=2))
#!
#!   \dontrun{
#!     library(polycor)
#!
#!     n <- 400
#!     x <- rnorm(n)
#!     y <- rnorm(n)
#!
#!     x1 <- (x + rnorm(n))/2
#!     x2 <- (x + rnorm(n))/2
#!     x3 <- (x + rnorm(n))/2
#!     x4 <- (x + rnorm(n))/2
#!
#!     y1 <- (y + rnorm(n))/2
#!     y2 <- (y + rnorm(n))/2
#!     y3 <- (y + rnorm(n))/2
#!     y4 <- (y + rnorm(n))/2
#!
#!     dat <- data.frame(x1, x2, x3, x4, y1, y2, y3, y4)
#!
#!     x1 <- ordered(as.integer(x1 > 0))
#!     x2 <- ordered(as.integer(x2 > 0))
#!     x3 <- ordered(as.integer(x3 > 1))
#!     x4 <- ordered(as.integer(x4 > -1))
#!
#!     y1 <- ordered(as.integer(y1 > 0))
#!     y2 <- ordered(as.integer(y2 > 0))
#!     y3 <- ordered(as.integer(y3 > 1))
#!     y4 <- ordered(as.integer(y4 > -1))
#!
#!     odat <- data.frame(x1, x2, x3, x4, y1, y2, y3, y4)
#!
#!     xcor <- cor(dat)
#!     pcor <- cor(odat)
#!     hcor <- hetcor(odat, ML=TRUE, std.err=FALSE)$correlations
#!     ncor <- nearcor(hcor)$cor
#!
#!     try(factanal(covmat=xcor, factors=2, n.obs=n))
#!     try(factanal(covmat=pcor, factors=2, n.obs=n))
#!     try(factanal(covmat=hcor, factors=2, n.obs=n))
#!     try(factanal(covmat=ncor, factors=2, n.obs=n))
#!   }
#! }
#! \keyword{algebra}
#! \keyword{array}

nearcor <- function(  # Computes the nearest correlation matrix to an approximate correlation matrix, i.e. not positive semidefinite.
  R                   # n-by-n approx correlation matrix
, eig.tol   = 1.0e-6  # defines relative positiveness of eigenvalues compared to largest
, conv.tol  = 1.0e-7  # convergence tolerance for algorithm
, posd.tol  = 1.0e-8  # tolerance for enforcing positive definiteness
, maxits    = 100     # maximum number of iterations allowed
, verbose   = FALSE   # set to TRUE to verbose convergence

                      # RETURNS list of class nearcor with components cor, iterations, converged
){
  if (!(is.numeric(R) && is.matrix(R) && identical(R,t(R))))
    stop('Error: Input matrix R must be square and symmetric')

  # Inf norm
  inorm <- function(x)max(rowSums(abs(x)))
  # Froebenius norm
  fnorm <- function(x)sqrt(sum(diag(t(x) %*% x)))

  n <- ncol(R)
  U <- matrix(0, n, n)
  Y <- R
  iter <- 0

  while (TRUE){
      T <- Y - U

      # PROJECT ONTO PSD MATRICES
      e <- eigen(Y, symmetric=TRUE)
      Q <- e$vectors
      d <- e$values
      D <- diag(d)

      # create mask from relative positive eigenvalues
      p <- (d>eig.tol*d[1]);

      # use p mask to only compute 'positive' part
      X <- Q[,p,drop=FALSE] %*% D[p,p,drop=FALSE] %*% t(Q[,p,drop=FALSE])

      # UPDATE DYKSTRA'S CORRECTION
      U <- X - T

      # PROJECT ONTO UNIT DIAG MATRICES
      X <- (X + t(X))/2
      diag(X) <- 1

      conv <- inorm(Y-X) / inorm(Y)
      iter <- iter + 1
      if (verbose)
        cat("iter=", iter, "  conv=", conv, "\n", sep="")

      if (conv <= conv.tol){
        converged <- TRUE
        break
      }else if (iter==maxits){
        warning(paste("nearcor did not converge in", iter, "iterations"))
        converged <- FALSE
        break
      }
      Y <- X
  }
  X <- (X + t(X))/2
  # begin from posdefify(sfsmisc)
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
      d[d < Eps] <- Eps
      Q <- e$vectors
      o.diag <- diag(X)
      X <- Q %*% (d * t(Q))
      D <- sqrt(pmax(Eps, o.diag)/diag(X))
      X[] <- D * X * rep(D, each = n)
  }
  # end from posdefify(sfsmisc)
  # force symmetry
  X <- (X + t(X))/2
  diag(X) <- 1
  ret <- list(cor=X, fnorm=fnorm(R-X), iterations=iter, converged=converged)
  class(ret) <- "nearcor"
  ret
}
