\name{fisherz}
\alias{fisherz}
\alias{fisherz2r}
\alias{r.con}
\alias{r2t} 

\title{ Fisher r to z and z to r and confidence intervals}
\description{convert a correlation to a z score or z to r using the Fisher transformation or find the confidence intervals for a specified correlation
}
\usage{
fisherz(rho)
fisherz2r(z)
r.con(rho,n,p=.95,twotailed=TRUE)
r2t(rho,n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{ a Pearson r }
  \item{z}{A Fisher z}
  \item{n}{Sample size for confidence intervals}
  \item{p}{Confidence interval}
  \item{twotailed}{Treat p as twotailed p}
}

\value{ z value corresponding to r  (fisherz) \\
r corresponding to z (fisherz2r) \\
lower and upper p confidence intervals (r.con) \\
t with n-2 df (r2t)
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
