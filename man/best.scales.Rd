\name{bestScales}
\alias{bestItems}
\alias{bestScales}
\title{A bootstrap aggregation (bagging) function for choosing most predictive items} 

\description{\code{\link{bestScales}} forms scales from the items/scales most correlated with a particular criterion and then cross validates on a hold out sample. This is repeated n.iter times using basic bootstrap aggregation (bagging) techniques.   Given a dictionary of item content, \code{\link{bestScales}} will sort by criteria correlations and display the item content. Options for bagging (bootstrap aggregation) are included.
}
\usage{
bestScales(x,criteria,cut=.1,n.item =10,overlap=FALSE,dictionary=NULL,check=FALSE,
   impute="none", n.iter =1,frac=.9,p.keyed,digits=2)
bestItems(x,criteria=1,cut=.3, abs=TRUE, dictionary=NULL,check=FALSE,digits=2) 
}
\arguments{
\item{x}{A data matrix or data frame depending upon the function.}
\item{criteria}{Which variables (by name or location) should be the empirical target for bestScales and bestItems.  May be a separate object. }
\item{cut}{Return all values in abs(x[,c1]) > cut.}
\item{abs}{if TRUE, sort by absolute value in bestItems}
\item{dictionary}{a data.frame with rownames corresponding to rownames in the f$loadings matrix or colnames of the data matrix or correlation matrix, and entries (may be multiple columns) of item content.}
\item{check}{if TRUE, delete items with no variance}
\item{n.item}{How many items make up an empirical scale}
\item{overlap}{Are the correlations with other criteria fair game for bestScales}
\item{impute}{When finding the best scales, and thus the correlations with the criteria,  how should we handle missing data?  The default is to drop missing items.}
\item{n.iter}{Replicate the best item function n.iter times, sampling frac of the cases each time, and validating on 1 - frac of the cases.}
\item{frac}{What fraction of the subjects should be used for the derivation sample.  frac is set to one if n.iter = 1}
\item{p.keyed}{The proportion of replications needed to include items in the best keys}
\item{digits}{round to digits}
  
}
\details{
\code{\link{bestScales}} will find up to n.items that have absolute correlations with a criterion greater than cut.  If the overlap option is FALSE (default) the other criteria are not used.  This is an example of ``dust bowl empiricism" in that there is no latent construct being measured, just those items that most correlate with a set of criteria. The empirically identified items are then formed into scales (ignoring concepts of internal consistency) which are then correlated with the criteria.  

Clearly, \code{\link{bestScales}} is capitalizing on chance associations.  Thus, we can validate the empirical scales by deriving them on a fraction of the total number of subjects, and cross validating on the remaining subjects. (K-fold cross validation).  This is done n.iter times to show the variability of these empirically derived scale validities.  This can only be done if x is a data matrix, not a correlation matrix.  For very large data sets (e.g., those from SAPA) these scales seem very stable. 

This is effectively a straight forward application of 'bagging' (bootstrap aggregation) and machine learning.  

The criteria can be the colnames of elements of x, or can be a separate data.frame. 

\code{\link{bestItems}} and \code{\link{lookup}} are simple helper functions to summarize correlation matrices or factor loading matrices.  \code{\link{bestItems}} will sort the specified column (criteria) of x on the basis of the (absolute) value of the column.  The return as a default is just the rowname of the variable with those absolute values > cut.   If there is a dictionary of item content and item names, then include the contents as a two column matrix with rownames corresponding to the item name and then as many fields as desired for item content. (See the example dictionary \code{\link{bfi.dictionary}}).
}

\value{
\code{\link{bestScales}} returns the correlation of the empirically constructed scale with each criteria and the items used in the scale.  If a dictionary is specified, it also returns a list (value) that shows the item content. Also returns the keys.list so that scales can be found using \code{\link{cluster.cor}} or \code{\link{scoreItems}}.  If using replications (bagging) then it also returns the best.keys , a list suitable for scoring.

The best.keys object is a list of items (with keying information) that may be used in subsequent analyses.
 
\code{\link{bestItems}} returns a sorted list of factor loadings or correlations with the labels as provided in the dictionary.

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
\code{\link{fa}}, \code{\link{iclust}},\code{\link{principal}}
}
\examples{
bs <- bestScales(bfi,criteria=cs(gender,age,education),n.iter=10,dictionary=bfi.dictionary[1:3])
bs


}
\keyword{ models }
\keyword{multivariate }
\keyword{tree }
