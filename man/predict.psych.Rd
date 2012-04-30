\name{predict.psych}
\alias{predict.psych}
\title{Prediction function for factor analysis or principal components
}
\description{Finds predicted factor/component scores from a factor analysis or components analysis of data set A predicted to data set B.  Predicted factor scores use the weights matrix used to find estimated factor scores, predicted components use the loadings matrix. Scores are either standardized with respect to the prediction sample or based upon the original data.
}
\usage{
\method{predict}{psych}(object, data,old.data,...)
}

\arguments{
  \item{object}{the result of a factor analysis or principal components analysis of data set A}
  \item{data}{Data set B, of the same number of variables as data set A.}
  \item{old.data}{if specified, the data set B will be standardized in terms of values from the old data.  This is probably the preferred option. }
  \item{...}{More options to pass to predictions }
}

\value{
Predicted factor/components scores.  The scores are based upon standardized items where the standardization is either that of the original data (old.data) or of the prediction set.  This latter case can lead to confusion if just a small number of predicted scores are found.  
}

\author{William Revelle
}
\note{Thanks to Reinhold Hatzinger for the suggestion and request}
\seealso{ \code{\link{fa}}, \code{\link{principal}}
}
\examples{
set.seed(42)
x <- sim.item(12,500)
f2 <- fa(x[1:250,],2,scores="regression")  # a two factor solution
p2 <- principal(x[1:250,],2,scores=TRUE)  # a two component solution
round(cor(f2$scores,p2$scores),2)  #correlate the components and factors from the A set
#find the predicted scores (The B set)
pf2 <- predict(f2,x[251:500,],x[1:250,])  #use the original data for standardization values 
pp2 <- predict(p2,x[251:500,],x[1:250,])  #standardized based upon the first set 
round(cor(pf2,pp2),2)   #find the correlations in the B set
#test how well these predicted scores match the factor scores from the second set
fp2 <- fa(x[251:500,],2,scores=TRUE)
round(cor(fp2$scores,pf2),2)

pf2.n <- predict(f2,x[251:500,])  #Standardized based upon the new data set
round(cor(fp2$scores,pf2.n))   #predict factors of set two from factors of set 1, factor order is arbitrary


#note that the signs of the factors in the second set are arbitrary
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ models }
