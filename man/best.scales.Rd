\name{bestScales}
\alias{bestItems}
\alias{bestScales}
\title{A bootstrap aggregation (bagging) function for choosing most predictive items} 

\description{\code{\link{bestScales}} forms scales from the items/scales most correlated with a particular criterion and then cross validates on a hold out sample using unit weighted scales.  This may be repeated n.iter times using either basic bootstrap aggregation (bagging) techniques or K-fold cross validation. Thus, the technique is known as BISCUIT (Best Items Scales that are Cross validated, Unit weighted, Informative, and Transparent).  Given a dictionary of item content, \code{\link{bestScales}} will sort by criteria correlations and display the item content. Options for bagging (bootstrap aggregation) are included.
}
\usage{
bestScales(x,criteria,cut=.1,n.item =10, overlap=FALSE, dictionary=NULL, check=FALSE, 
impute="none", n.iter =1, folds=1,p.keyed=.9,log.p=FALSE,digits=2)
bestItems(x,criteria=1,cut=.3, abs=TRUE, dictionary=NULL,check=FALSE,digits=2) 
}
\arguments{
\item{x}{A data matrix or data frame depending upon the function.}
\item{criteria}{Which variables (by name or location) should be the empirical target for bestScales and bestItems.  May be a separate object. }
\item{cut}{Return all values in abs(x[,c1]) > cut.}
\item{abs}{if TRUE, sort by absolute value in bestItems}
\item{dictionary}{a data.frame with rownames corresponding to rownames in the f$loadings matrix or colnames of the data matrix or correlation matrix, and entries (may be multiple columns) of item content.}
\item{check}{if TRUE, delete items with no variance}
\item{n.item}{How many items make up an empirical scale. }
\item{overlap}{Are the correlations with other criteria fair game for bestScales}
\item{impute}{When finding the best scales, and thus the correlations with the criteria,  how should we handle missing data?  The default is to drop missing items. (That is to say, to use pairwise complete correlations.)}
\item{n.iter}{How many times to perform a bootstrap estimate. Replicate the best item function n.iter times, sampling roughly 1-1/e  of the cases each time, and validating on the remaining 1/e of the cases for each iteration.}
\item{folds}{If folds > 1, this is k-folds validation.  Note, set  n.iter > 1  to do bootstrap aggregation, or set folds > 1  to do k-folds.   }
\item{p.keyed}{The proportion of replications needed to include items in the final best keys.}
\item{log.p}{Select items based upon the log of the probability of the correlations.  This will only have an effect if the number of pairwise cases differs drastically from pair to pair.  }
\item{digits}{round to digits when showing output.}
  
}
\details{
\code{\link{bestScales}} will find up to n.items per criterion that have absolute correlations with a criterion greater than cut.  If the overlap option is FALSE (default) the other criteria are not used.  This is an example of ``dust bowl empiricism" in that there is no latent construct being measured, just those items that most correlate with a set of criteria. The empirically identified items are then formed into scales (ignoring concepts of internal consistency) which are then correlated with the criteria.  

Clearly, \code{\link{bestScales}} is capitalizing on chance associations.  Thus, we should validate the empirical scales by deriving them on a fraction of the total number of subjects, and cross validating on the remaining subjects. (This is known both as K-fold cross validation and bagging.  Both may be done).  If folds > 1, then a k-fold cross validation is done.  This removes 1/k (a fold) from the sample for the derivation sample and validates on that remaining fold. This is done k-folds times.  Traditional cross validation would thus be a k-fold with k =2.  More modern applications seem to prefer k =10 to have 90\% derivation sample and a 10\% cross validation sample.

The alternative, known as 'bagging' is to do a bootstrap sample (which because it is sampling with replacement will typically extract 1- 1/e = 63.2\% of the sample) for the derivation sample (the bag) and then validate it on the remaining 1/e of the sample (the out of bag).  This is done n.iter times. This should be repeated multiple times (n.iter > 1, e.g. n.iter=1000) to get stable cross validations.

One can compare the validity of these two approaches by  trying each.  The average predictability of the n.iter samples are shown as are the average validity of the cross validations.   This can only be done if x is a data matrix/ data.frame, not a correlation matrix.  For very large data sets (e.g., those from SAPA) these scales seem very stable. 

\code{\link{bestScales}} is effectively a straight forward application of 'bagging' (bootstrap aggregation) and machine learning as well as k-fold validation. 

The criteria can be the colnames of elements of x, or can be a separate data.frame. 

\code{\link{bestItems}} and \code{\link{lookup}} are simple helper functions to summarize correlation matrices or factor loading matrices.  \code{\link{bestItems}} will sort the specified column (criteria) of x on the basis of the (absolute) value of the column.  The return as a default is just the rowname of the variable with those absolute values > cut.   If there is a dictionary of item content and item names, then include the contents as a two column matrix with rownames corresponding to the item name and then as many fields as desired for item content. (See the example dictionary \code{\link{bfi.dictionary}}).
}

\value{
\code{\link{bestScales}} returns the correlation of the empirically constructed scale with each criteria and the items used in the scale.  If a dictionary is specified, it also returns a list (value) that shows the item content. Also returns the keys.list so that scales can be found using \code{\link{cluster.cor}} or \code{\link{scoreItems}}.  If using replications (bagging or kfold) then it also returns the best.keys , a list suitable for scoring.

The best.keys object is a list of items (with keying information) that may be used in subsequent analyses.  These ``best.keys" are formed into scales for the ``final.valid" object which reports how well the best.keys work on the entire sample.  This is, of course, not cross validated.  
 
\code{\link{bestItems}} returns a sorted list of factor loadings or correlations with the labels as provided in the dictionary.

The stats object can be used to create \code{\link{error.dots}} plots to show the mean estimate and the standard error of the estimates.   See  the examples.

}

\note{

To create a dictionary, create an object with row names as the item numbers, and the columns as the item content.  See the \code{link{bfi.dictionary}} as an example.
}
\references{
Revelle, W.  (in preparation) An introduction to psychometric theory with applications in {R}. Springer.  (Available online at \url{https://personality-project.org/r/book}). 

}
\author{William Revelle}
\note{Although empirical scale construction is appealing, it has the basic problem of capitalizing on chance.  Thus, be careful of over interpreting the results unless working with large samples.  Iteration and bootstrapping aggregation (bagging) gives information on the stability of the solutions.
}

\seealso{
\code{\link{fa}}, \code{\link{iclust}},\code{\link{principal}}, \code{\link{error.dots}}
}
\examples{
#This is an example of 'bagging' (bootstrap aggregation)
bestboot <- bestScales(bfi,criteria=cs(gender,age,education), 
 n.iter=10,dictionary=bfi.dictionary[1:3])
bestboot
#compare with 10 fold cross validation 
tenfold <- bestScales(bfi,criteria=cs(gender,age,education),fold=10,dictionary=bfi.dictionary[1:3])
tenfold
#Then, to display the results graphically
#Note that we scale the two graphs with the same x.lim values
error.dots(bestboot,eyes=TRUE,xlim=c(0,.4))
error.dots(tenfold,add=TRUE,pch=16,xlim=c(0,.4))

}
\keyword{ models }
\keyword{multivariate }
\keyword{tree }
