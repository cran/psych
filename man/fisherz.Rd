\name{fisherz}
\alias{fisherz}
\alias{fisherz2r}
\alias{r.con}\alias{r2c}
\alias{r2t}
\alias{t2r}
\alias{g2r}
\alias{chi2r}
\alias{r2chi}
\alias{cor2cov}
\alias{r2p}

\title{Transformations of r, d, and t  including Fisher r to z and z to r and confidence intervals}
\description{Convert a correlation to a z or t, or d, or chi or covariance matrix
  or z to r using the Fisher transformation or find the confidence intervals for a specified correlation.  r2d converts a correlation to an effect size (Cohen's d) and d2r converts a d into an r. g2r converts Hedge's g to a correlation.   t2r converts a t test to r, r2t converts a correlation to a t-test value. chi2r converts a chi square to r, r2chi converts it back.  r2c and cor2cov convert a correlation matrix to a covariance matrix. d2t and t2d convert cohen's d into a t and a t into a cohen d.  See \code{\link{cohen.d}} for other conversions.
}
\usage{
fisherz(rho)
fisherz2r(z)
r.con(rho,n,p=.95,twotailed=TRUE)
r2t(rho,n)
t2r(t,df)
g2r(g,df,n)
chi2r(chi2,n)
r2chi(rho,n)
r2c(rho,sigma)
cor2cov(rho,sigma)
r2p(rho,n)




}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{ a Pearson r }
  \item{z}{A Fisher z}
  \item{n}{Sample size for confidence intervals}
  \item{df}{degrees of freedom for t, or g}
  \item{p}{Confidence interval}
  \item{twotailed}{Treat p as twotailed p}
  \item{g}{An effect size (Hedge's g)}
  \item{t}{A student's t value}
  \item{chi2}{A chi square}
  \item{sigma}{a vector of standard deviations to be used to convert a correlation matrix to a covariance matrix}
}

\value{ 
\item{z}{ value corresponding to r  (fisherz)} 
\item{r}{r corresponding to z (fisherz2r)} 
\item{r.con}{lower and upper p confidence intervals (r.con)}
\item{t}{t with n-2 df (r2t)}
\item{r}{r corresponding to effect size d or d corresponding to r.}
\item{r2c}{r2c  is the reverse of the cor2con function of base R.  It just converts a correlation matrix to the corresponding covariance matrix given a vector of standard deviations.}
}


\author{
Maintainer: William Revelle  \email{revelle@northwestern.edu }
}

\examples{

 n <- 30
 r <- seq(0,.9,.1)
 d <- r2d(r)
 rc <- matrix(r.con(r,n),ncol=2)
 t <- r*sqrt(n-2)/sqrt(1-r^2)
 p <- (1-pt(t,n-2))*2
 r1 <- t2r(t,(n-2))
 r2 <- d2r(d)
 chi <- r2chi(r,n)
 r3 <- chi2r(chi,n)
 r.rc <- data.frame(r=r,z=fisherz(r),lower=rc[,1],upper=rc[,2],t=t,p=p,d=d,
     chi2=chi,d2r=r2,t2r=r1,chi2r=r3)
 round(r.rc,2)


}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
