\name{unidim}
\alias{unidim}
\title{Several indices of the unidimensionality of a set of variables.}

\description{There are a variety of ways of assessing whether a set of items measures one latent trait.  \code{\link{unidim}} is just one more way.  If a one factor model holds in the data, then the factor analytic decomposition F implies that FF' should reproduce the correlations with communalities along the diagonal. In this case, the fit FF' should be identical to the correlation matrix minus the uniquenesses.  unidim is just the ratio of these two estimates.  The higher it is, the more the evidence for unidimensionality.

}
\usage{
unidim(x, keys.list = NULL, flip = FALSE)
}

\arguments{
  \item{x}{An input matrix or data frame.  If x is not a correlation matrix, then the correlations are found.}
  \item{keys.list}{If specified,  then a number of scales can be tested at once. (See \code{\link{scoreItems}} for a similar procedure.)}

  \item{flip}{If TRUE, then items will be keyed based upon their loadings on the first factor.  Automatically done if key.list is NULL.}

}
\details{
This is an exploratory index that is still under development.  A number of test cases suggest that it provides high values when the data are in fact unidimensional, low values when they are not.

The logic is deceptively simple:  Unidimensionality implies that a one factor model of the data fits the covariances of the data.  If this is the case, then factor model implies R = FF' + U2 will have residuals of 0.  Similarly, this also implies that the observed correlations will equal the model.  Thus, the sum of the observed corelations (with the diagonal replaced by the communalities) should match the factor model.  Compare these two models:  R - U2  versus FF'.   
}
\value{
  \item{uni.orig }{The raw value of the unidimensional criterion}
 \item{uni.adj }{The unidimensional criterion when items are keyed in positive direction.}
\item{fit1}{The off diagonal fit from  \code{\link{fa}}  }
\item{alpha}{Standardized alpha of the keyed items (after appropriate reversals)}
\item{av.r}{The average interitem correlation of the keyed items.}
\item{raw.model}{The ratio of the FF' model to the sum(R)}
\item{adj.model}{The ratio of the FF' model to the sum(R) when items are flipped.}
item{raw.total}{The ratio of the sum(R - uniqueness)/sum(R)}
item{adj.total}{Same ratio with flipped items}

}

\author{William Revelle}

\note{A perhaps interesting idea but still an exploratory statistic. Treat with appropriate caution. 
}


\seealso{ \code{\link{fa}} for factor analysis, \code{\link{omega}} for reliability.

}
\examples{
#test the unidimensionality of the five factors of the bfi data set.

keys.list <-
  list(agree=c("-A1","A2","A3","A4","A5"),conscientious=c("C1","C2","C3","-C4","-C5"),
extraversion=c("-E1","-E2","E3","E4","E5"),neuroticism=c("N1","N2","N3","N4","N5"),
openness = c("O1","-O2","O3","O4","-O5"), all = 
c("-A1","A2","A3","A4","A5","C1","C2","C3","-C4","-C5","-E1","-E2","E3","E4","E5","N1" 
,"N2","N3","N4","N5","O1","-O2","O3","O4","-O5") )
unidim(psychTools::bfi,keys.list) 

#Try a known 3 factor structure
x <- sim.minor(nfact=3,bipolar=FALSE)
unidim(x$model) 
keys.list <- list(first =c(1:4),second = 5:8,third=9:12,all=1:12)
unidim(x$model,keys.list)

x <- sim.minor(nfact=3)
unidim(x$model,keys.list,flip=TRUE)

}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line