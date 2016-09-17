\name{fisherz}
\alias{fisherz}
\alias{fisherz2r}
\alias{r.con}
\alias{r2t}
\alias{r2d}
\alias{d2r}
\alias{t2r}
\alias{g2r}
\alias{chi2r}
\alias{r2chi}

\title{Transformations of r including Fisher r to z and z to r and confidence intervals}
\description{Convert a correlation to a z or t, or d, or chi
 score or z to r using the Fisher transformation or find the confidence intervals for a specified correlation.  r2d converts a correlation to an effect size (Cohen's d) and d2r converts a d into an r. g2r converts Hedge's g to a correlation.   t2r converts a t test to r, r2t converts a correlation to a t-test value. chi2r converts a chi square to r, r2chi converts it back. 
}
\usage{
fisherz(rho)
fisherz2r(z)
r.con(rho,n,p=.95,twotailed=TRUE)
r2t(rho,n)
r2d(rho)
d2r(d)
t2r(t,df)
g2r(g,df,n)
chi2r(chi,n)
r2chi(rho,n)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{ a Pearson r }
  \item{z}{A Fisher z}
  \item{n}{Sample size for confidence intervals}
  \item{df}{degrees of freedom for t, or g}
  \item{p}{Confidence interval}
  \item{twotailed}{Treat p as twotailed p}
  \item{d}{an effect size (Cohen's d)}
  \item{g}{An effect size (Hedge's g)}
  \item{t}{A student's t value}
  \item{chi}{A chi square}
}

\value{ z value corresponding to r  (fisherz) \\
r corresponding to z (fisherz2r) \\
lower and upper p confidence intervals (r.con) \\
t with n-2 df (r2t)
r corresponding to effect size d or d corresponding to r.
}


\author{
Maintainer: William Revelle  \email{revelle@northwestern.edu }
}

\examples{

cors <- seq(-.9,.9,.1)
zs <- fisherz(cors)
rs <- fisherz2r(zs)
round(zs,2)
 n <- 30
 r <- seq(0,.9,.1)
 rc <- matrix(r.con(r,n),ncol=2)
 t <- r*sqrt(n-2)/sqrt(1-r^2)
 p <- (1-pt(t,n-2))/2
 r.rc <- data.frame(r=r,z=fisherz(r),lower=rc[,1],upper=rc[,2],t=t,p=p)
 round(r.rc,2)

}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
