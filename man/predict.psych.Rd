\name{predict.psych}
\alias{predict.psych}
\title{Prediction function for factor analysis, principal components (pca), bestScales
}
\description{Finds predicted factor/component scores from a factor analysis or principal components analysis (pca) of data set A predicted to data set B.  Predicted factor scores use the weights matrix used to find estimated factor scores, predicted components use the loadings matrix. Scores are either standardized with respect to the prediction sample or based upon the original data. Predicted scores from a bestScales model are based upon the statistics from the original sample.
}
\usage{
\method{predict}{psych}(object, data,old.data,options=NULL,missing=FALSE,impute="none",...)
}

\arguments{
  \item{object}{the result of a factor analysis, principal components analysis (pca) or bestScales of data set A}
  \item{data}{Data set B, of the same number of variables as data set A.}
  \item{old.data}{if specified, the data set B will be standardized in terms of values from the old data.  This is probably the preferred option. This is done automatically if object is from \code{\link{bestScales}} }
  \item{options}{scoring options for bestScales objects 
       ("best.keys","weights","optimal.keys","optimal.weights")}
\item{missing}{If missing=FALSE, cases with missing data are given NA scores, otherwise they are given the values based upon the wts x complete data }
\item{impute}{Should missing cases be replaced by "means", "medians" or treated as missing 
("none" is the default}
  \item{...}{More options to pass to predictions }
}

\value{
Predicted factor/components/criteria scores.  If predicting from either \code{\link{fa}} or \code{\link{pca}},the scores are based upon standardized items where the standardization is either that of the original data (old.data) or of the prediction set.  This latter case can lead to confusion if just a small number of predicted scores are found.   

If the object is from \code{\link{bestScales}}, unit weighted scales are found (by default) using the best.keys and the predicted scores are then put into the metric of the means and standard deviations of the derivation sample. Other scoring key options may be specified using the "options" parameter.  Possible values are best.keys","weights","optimal.keys","optimal.weights".  See \code{\link{bestScales}} for details.

By default, predicted scores are found by the matrix product of the standardized data with the factor or regression weights.  If missing is TRUE, then the predicted scores are the mean of the standardized data x weights for those data points that are not NA.

}

\author{William Revelle
}
\note{Thanks to Reinhold Hatzinger for the suggestion and request and to Sarah McDougald for the bestScales prediction.}
\seealso{ \code{\link{fa}}, \code{\link{principal}}, \code{\link{bestScales}}
}
\examples{
set.seed(42)
x <- sim.item(12,500)
f2 <- fa(x[1:250,],2,scores="regression")  # a two factor solution
p2 <- principal(x[1:250,],2,scores=TRUE)  # a two component solution
round(cor(f2$scores,p2$scores),2)  #correlate the components and factors from the A set
#find the predicted scores (The B set)
pf2 <- predict(f2,x[251:500,],x[1:250,]) 

  #use the original data for standardization values 
pp2 <- predict(p2,x[251:500,],x[1:250,]) 
 #standardized based upon the first set 
round(cor(pf2,pp2),2)   #find the correlations in the B set
#test how well these predicted scores match the factor scores from the second set
fp2 <- fa(x[251:500,],2,scores=TRUE)
round(cor(fp2$scores,pf2),2)

pf2.n <- predict(f2,x[251:500,])  #Standardized based upon the new data set
round(cor(fp2$scores,pf2.n))   
   #predict factors of set two from factors of set 1, factor order is arbitrary


#note that the signs of the factors in the second set are arbitrary
\donttest{
#predictions from bestScales
#the derivation sample
bs <- bestScales(psychTools::bfi[1:1400,], cs(gender,education,age),folds=10,p.keyed=.5) 
pred <- predict(bs,psychTools::bfi[1401:2800,]) #The prediction sample
cor2(pred,psychTools::bfi[1401:2800,26:28] ) #the validity of the prediction
summary(bs) #compare with bestScales cross validations
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ models }
