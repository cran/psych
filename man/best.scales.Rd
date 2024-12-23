\name{bestScales}
\alias{bestItems}
\alias{bestScales}
\alias{BISCUIT}
\alias{biscuit}
\alias{BISCWIT}
\alias{biscwit}


\title{A bootstrap aggregation function for choosing most predictive unit weighted items} 

\description{\code{\link{bestScales}} forms scales from the items/scales most correlated with a particular criterion and then cross validates on a hold out sample using unit weighted scales.  This may be repeated n.iter times using either basic bootstrap aggregation (bagging) techniques or K-fold cross validation. Thus, the technique is known as \code{\link{BISCUIT}} (Best Items Scales that are Cross validated, Unit weighted, Informative, and Transparent).  Given a dictionary of item content, \code{\link{bestScales}} will sort by criteria correlations and display the item content. Options for bagging (bootstrap aggregation) are included. An alternative to unit weighting is to weight items by their zero order correlations (cross validated) with the criteria. This weighted version is called \code{\link{BISCWIT}} and is an optional output. 
}
\usage{
bestScales(x,criteria,min.item=NULL,max.item=NULL, delta = 0,
           cut=.1, n.item =10, wtd.cut=0, wtd.n=10, 
          n.iter =1, folds=1, p.keyed=.9,
          overlap=FALSE, dictionary=NULL, check=TRUE,labels=NULL,
           impute="none",log.p=FALSE,digits=2)

bestItems(x,criteria=1,cut=.1, n.item=10, abs=TRUE, 
   dictionary=NULL,check=FALSE,digits=2,use="pairwise",method="pearson") 
}
\arguments{
\item{x}{A data matrix or data frame depending upon the function.}
\item{criteria}{Which variables (by name or location) should be the empirical target for bestScales and bestItems.  May be a separate object. }
\item{min.item}{Find unit weighted and correlation weighted scales from min.item to max.item}
\item{max.item}{These are all summarized in the final.multi.valid object}
\item{delta}{Return items where the predicted r + delta * se of r < max value}
\item{cut}{Return all values in abs(x[,c1]) > cut.}
\item{wtd.cut}{When finding the weighted scales, use all items with zero 
    order correlations > wtd.cut}
\item{wtd.n}{When finding the weighted scales, use the  wtd.n items that
   are > than wtd.cut}

\item{abs}{if TRUE, sort by absolute value in bestItems}
\item{dictionary}{a data.frame with rownames corresponding to rownames in the f$loadings 
matrix or colnames of the data matrix or correlation matrix, and entries (may be multiple
 columns) of item content.}
\item{check}{if TRUE, delete items with no variance}
\item{labels}{labels for criteria if desired}
\item{n.item}{How many items make up an empirical scale, or (bestItems, show the best
 n.items) }
\item{overlap}{Are the correlations with other criteria fair game for bestScales}
\item{impute}{When finding the best scales, and thus the correlations with the criteria,  how should we handle missing data?  The default is to drop missing items. (That is to say, to use pairwise complete correlations.)}
\item{n.iter}{How many times to perform a bootstrap estimate. Replicate the best item function n.iter times, sampling roughly 1-1/e  of the cases each time, and validating on the remaining 1/e of the cases for each iteration.}
\item{folds}{If folds > 1, this is k-folds validation.  Note, set  n.iter > 1  to do bootstrap aggregation, or set folds > 1  to do k-folds.   }
\item{p.keyed}{The proportion of replications needed to include items in the final best keys.}
\item{log.p}{Select items based upon the log of the probability of the correlations.  This will only have an effect if the number of pairwise cases differs drastically from pair to pair.  }

\item{digits}{round to digits when showing output.}
\item{use}{How to handle missing data.  Defaults to "pairwise"}
\item{method}{Which correlation to find.  Defaults to "pearson"}
  
}
\details{
There are a number of procedures that can be used for predicting criteria from a set of predictors.  The generic term for this is "machine learning" or "statistical learning".  The basic logic of these procedures is to find a set of items that best predict a criteria according to some fit statistic and then cross validate these items numerous times.  "lasso" regression (least absolute shrinkage and selection) is one such example. \code{\link{bestScales}} differs from these procedures by unit weighting items chosen from their zero order correlations with the criteria rather than weighting the partial correlations ala regression.  This is an admittedly simple procedure that takes into account the well known finding (Wilks, 1938;  Wainer, 1976; Dawes, 1979; Waller, 2008) that although regression weights are optimal for any particular data set, unit weights are almost as good (fungible) and more robust across sample variation.  

Following some suggestions, we have added the ability to find scales where items are weighted by their zero order correlations with the criteria.   This is effectively a comprimise between unit weighting and regression weights (where the weights are the zero order correlations times the inverse of the correlation matrix). This weighted version may be thought of as \code{\link{BISCWIT}} in contrast to the unit weighted version \code{\link{BISCUIT}}.

To be comparable to other ML algorithms, we now consider multiple solutions (for number of items >= min.item to max.item).  The final scale consists of the number  items which maximize  the validity or at least are within delta * standard error of r of the maximum.  
 
Thus, \code{\link{bestScales}} will find up to n.items per criterion that have absolute correlations with a criterion greater than cut.  If the overlap option is FALSE (default) the other criteria are not used.  This is an example of ``dust bowl empiricism" in that there is no latent construct being measured, just those items that most correlate with a set of criteria. The empirically identified items are then formed into scales (ignoring concepts of internal consistency) which are then correlated with the criteria.  

Clearly, \code{\link{bestScales}} is capitalizing on chance associations.  Thus, we should validate the empirical scales by deriving them on a fraction of the total number of subjects, and cross validating on the remaining subjects. (This is known both as K-fold cross validation and bagging.  Both may be done).  If folds > 1, then a k-fold cross validation is done.  This removes 1/k (a fold) from the sample for the derivation sample and validates on that remaining fold. This is done k-folds times.  Traditional cross validation would thus be a k-fold with k =2.  More modern applications seem to prefer k =10 to have 90\% derivation sample and a 10\% cross validation sample.

The alternative, known as 'bagging' is to do a bootstrap sample (which because it is sampling with replacement will typically extract 1- 1/e = 63.2\% of the sample) for the derivation sample (the bag) and then validate it on the remaining 1/e of the sample (the out of bag).  This is done n.iter times. This should be repeated multiple times (n.iter > 1, e.g. n.iter=1000) to get stable cross validations.

One can compare the validity of these two approaches by  trying each.  The average predictability of the n.iter samples are shown as are the average validity of the cross validations.   This can only be done if x is a data matrix/ data.frame, not a correlation matrix.  For very large data sets (e.g., those from SAPA) these scales seem very stable. 

\code{\link{bestScales}} is effectively a straight forward application of 'bagging' (bootstrap aggregation) and machine learning as well as k-fold validation. 

The criteria can be the colnames of elements of x, or can be a separate data.frame. 

\code{\link{bestItems}} and \code{\link{lookup}} are simple helper functions to summarize correlation matrices or factor loading matrices.  \code{\link{bestItems}} will sort the specified column (criteria) of x on the basis of the (absolute) value of the column.  The return as a default is just the rowname of the variable with those absolute values > cut.   If there is a dictionary of item content and item names, then include the contents as a two column matrix with rownames corresponding to the item name and then as many fields as desired for item content. (See the example dictionary \code{\link[psychTools]{bfi.dictionary}}).

The derived model can be further validated against yet another hold out sample using the \code{\link{predict.psych}} function if given the best scale object and the new data set. 
}

\value{
\code{\link{bestScales}} returns the correlation of the empirically constructed scale with each criteria and the items used in the scale.  If a dictionary is specified, it also returns a list (value) that shows the item content. Also returns the keys.list so that scales can be found using \code{\link{cluster.cor}} or \code{\link{scoreItems}}.  If using replications (bagging or kfold) then it also returns the best.keys, a list suitable for scoring.

There are actually four keys lists reported.

best.keys are all the items used to form unit weighted scales with the restriction of n.item.  

weights may be used in the \code{\link{scoreWtd}} function to find scales based upon the raw correlation weights.

If the min.item and max.item options are used, then two more sets of weights are provided.

optimal.keys are a subset of the best.keys, taking just those items that increase the cross validation values up to the delta * se of the maximum.  This is a slightly more parsimonious set.

optimal.weights is analogous to optimal keys, but for supplies weights for just those items that are used to predict cross validation values up to delta * se of the maximum.

The best.keys object is a list of items (with keying information) that may be used in subsequent analyses.  These ``best.keys" are formed into scale scores for the ``final.valid" object which reports how well the best.keys work on the entire sample.  This is, of course, not cross validated.  Further cross validation can be done using the \code{\link{predict.psych}} function.

\item{scores}{ Are the unit weighted scores from the original items}
\item{best.keys}{ A key list of those items that were used in the unit weighting.}

\item{wtd.scores}{Are the zero-order correlation based scores.}

\item{weights}{ the scoring weights used}
\item{final.multi.valid}{An object with the unit weighted and correlation weighted correlations from low.step to high.step}
 


The print and summary output list a number of summary statistics for each criteria.  This is given for the default case (number of items fixed) and then if requested, the optimal values chosen from min.item to max.item:

The default statistics:

 \describe{
 \item{derivation mean}{Mean correlation of fixed length scale with the criteria, derivation sample}
 \item{derivation.sd}{The standard deviation of these estimates}
 \item{validation.m}{The mean cross validated correlations with the criteria}
 \item{validation.sd}{The standard deviations of these estimates}
 \item{final.valid}{The correlation of the pooled models with all the subjects}
 \item{final.wtd}{The correlation of the pooled weighted model with all subjects}
 \item{N.wtd}{Number of items used in the final weighted model}
 }

The optimal number of items statistics:

 \describe{
 \item{n}{The mean number of items meeting the criteria}
 \item{unit}{The mean derivation predictive valididy}
 \item{n.wtd}{the mean number of items used in the wtd scales}
 \item{wtd}{The mean derivation wtd correlaton}
 \item{valid.n}{the mean number of items in the cross validation sample}
 \item{valid.unit}{The mean cross validated unit weighted correlations}
 \item{valid.wtd.n}{The mean number of items used in the cross validated correlated weighed scale}
 \item{valid.wtd}{The mean cross validated weighted correlation with criteria}
 \item{n.final}{The optimal number of items on the final cross validation sample}
 \item{n.wtd.final}{The optimal number of weighted items on the final cross validation.}
\item{derviation.mean}{}
 }

\code{\link{bestItems}} returns a sorted list of factor loadings or correlations with the labels as provided in the dictionary. If given raw data, then the correlations with the criteria variables are found first according to "use" and "method".

The stats object can be used to create \code{\link{error.dots}} plots to show the mean estimate and the standard error of the estimates.   See  the examples.

The resulting best scales can be cross validated on a different hold out sample using \code{\link{crossValidation}}.  See the last example. 

}

\note{

Although \code{\link{bestScales}} was designed to form the best unit weighted scales, for large data sets, there seems to be some additional information in weighting by the average zero-order correlation. 

To create a dictionary, create an object with row names as the item numbers, and the columns as the item content.  See the \code{link{bfi.dictionary}} as an example.

This is a very computationally intensive function which can be speeded up considerably by using multiple cores and using the parallel package. The number of cores to use when doing bestScales may be specified using the options command. The greatest step in speed is going from 1 core to 2. This is about a 50% savings. Going to 4 cores seems to have about at 66% savings, and 8 a 75% savings. The number of parallel processes defaults to 2 but can be modified by using the options command: options("mc.cores"=4) will set the number of cores to 4.

}
\references{

Dawes, R.M. (1979) The robust beauty of improper linear models in decision making, American Psychologist, 34, 571-582.

Elleman, L. G., McDougald, S. K., Condon, D. M., & Revelle, W. 2020 (in press). That takes the BISCUIT: A comparative study of predictive accuracy and parsimony of four statistical learning techniques in personality data, with data missingness conditions.  European Journal of Psychological Assessment. (Preprint available at https://psyarxiv.com/tuqap/)

Revelle, W.  (in preparation) An introduction to psychometric theory with applications in {R}. Springer.  (Available online at \url{https://personality-project.org/r/book/}). 


Wainer, H. (1979) Estimating coefficients in linear models: It don't make no nevermind. Psychological Buletin, 83, 213-217.

Waller, N.G. (2008), Fungible weights in multiple regression.  Psychometrica, 73, 691-703. 

Wilks, S. S. (1938), Weighting systems for linear functions of correlated variables when there is no dependent variable. Psychometrika. 3. 23-40.	

}
\author{William Revelle}

\note{Although empirical scale construction is appealing, it has the basic problem of capitalizing on chance.  Thus, be careful of over interpreting the results unless working with large samples.  Iteration and bootstrapping aggregation (bagging) gives information on the stability of the solutions.

It is also important to realize that the variables that are most related to the criterion might be because of other, trivial reasons (e.g. height predicts gender)
}

\seealso{
\code{\link{fa}}, \code{\link{iclust}},\code{\link{principal}}, \code{\link{error.dots}}
}
\examples{
#This is an example of 'bagging' (bootstrap aggregation)
#not run in order to pass the timing tests for Debian at CRAN
#bestboot <- bestScales(bfi,criteria=cs(gender,age,education), 
# n.iter=10,dictionary=bfi.dictionary[1:3])
#bestboot
#compare with 10 fold cross validation 
\donttest{#don't test for purposes of speed in installation to pass the debian CRAN test
tenfold <- bestScales(bfi,criteria=cs(gender,age,education),fold=10,
           dictionary= bfi.dictionary[1:3])
tenfold

#Then, to display the results graphically
#Note that we scale the two graphs with the same x.lim values

 error.dots(tenfold,eyes=TRUE,pch=16,xlim=c(0,.4))
 # error.dots(bestboot,add=TRUE,xlim=c(0,.4))
#do this again, but this time display the scale fits from 1 to 15 Items
# tenfold <- bestScales(bfi,criteria=cs(gender,age,education),fold=10,
# dictionary= bfi.dictionary[1:3],min.item=1,max.item=15)
# matplot(tenfold$multi.validities$wtd.deriv,typ="b",
# xlab="Number of Items",main="Fit as a function of number of items")   #the wtd weights
# matpoints(tenfold$multi.validities$unit.deriv,typ="b",lwd=2) #the unit weights 

#an example of using crossValidation with bestScales
set.seed(42)
ss <- sample(1:2800,1400)
model <- bestScales(bfi[ss,],criteria=cs(gender,education,age))
original.fit <- crossValidation(model,bfi[ss,]) #the derivation set
cross.fit <- crossValidation(model,bfi[-ss,])  #the cross validation set
summary(original.fit)
summary(cross.fit)
}
}
\keyword{ models }
\keyword{multivariate }
\keyword{tree }
