\name{scoreIrt}
\alias{scoreIrt}

\alias{scoreIrt.1pl}
\alias{scoreIrt.2pl}
\alias{score.irt}
\alias{score.irt.poly}
\alias{score.irt.2}
\alias{irt.stats.like}
\alias{irt.tau}
\alias{irt.se}
\title{Find Item Response Theory (IRT) based scores for dichotomous or polytomous items}

\description{\code{\link{irt.fa}} finds Item Response Theory (IRT) parameters through factor analysis of the tetrachoric or polychoric correlations of dichtomous or polytomous items. \code{\link{scoreIrt}} uses these parameter estimates of discrimination and location to find IRT based scores for the responses. As many factors as found for the correlation matrix will be scored. \code{\link{scoreIrt.2pl}} will score lists of scales.
}
\usage{
scoreIrt(stats=NULL, items, keys=NULL,cut = 0.3,bounds=c(-4,4),mod="logistic") 
scoreIrt.1pl(keys.list,items,correct=.5,messages=FALSE,cut=.3,bounds=c(-4,4),
     mod="logistic")  #Rasch like scaling
scoreIrt.2pl(itemLists,items,correct=.5,messages=FALSE,cut=.3,bounds=c(-4,4),
   mod="logistic")  #2 pl scoring
#the next is an alias for scoreIrt both of which are wrappers for 
#     score.irt.2 and score.irt.poly
score.irt(stats=NULL, items, keys=NULL,cut = 0.3,bounds=c(-4,4),mod="logistic") 
 #the higher order call just calls one of the next two
  #for dichotomous items 
score.irt.2(stats, items,keys=NULL,cut = 0.3,bounds=c(-4,4),mod="logistic") 
  #for polytomous items
score.irt.poly(stats, items, keys=NULL, cut = 0.3,bounds=c(-4,4),mod="logistic")
    #to create irt like statistics for plotting
irt.stats.like(items,stats,keys=NULL,cut=.3)
   
irt.tau(x)    #find the tau values for the x object
irt.se(stats,scores=0,D=1.702)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stats}{Output from irt.fa is used for parameter estimates of location and discrimination. Stats may also be the  output from a normal factor analysis (fa) }
  \item{items}{The raw data, may be either dichotomous or polytomous.}
  \item{itemLists}{a list of items to be factored and scored for each scale, can be a keys.list as used in scoreItems or scoreIrt.1pl}
  \item{keys.list}{A list of items to be scored with keying direction  (see example)}
  \item{keys}{A keys matrix of which items should be scored for each factor}
  \item{cut}{Only items with discrimination values > cut will be used for scoring.}
  \item{x}{The raw data to be used to find the tau parameter in irt.tau}
  \item{bounds}{The lower and upper estimates for the fitting function}
  \item{mod}{Should a logistic or normal model be used in estimating the scores?}
  \item{correct}{What value should be used for continuity correction when finding the
     tetrachoric or polychoric correlations when using irt.fa}
  \item{messages}{Should messages be suppressed when running multiple scales?}
  \item{scores}{A single score or a vector of scores to find standard errors}
  \item{D}{The scaling function for the test information statistic used in irt.se}
}
\details{Although there are more elegant ways of finding subject scores given a set of item locations (difficulties) and discriminations, simply finding that value of theta \eqn{\theta} that best fits the equation \eqn{P(x|\theta) = 1/(1+exp(\beta(\delta - \theta) )} for a score vector X, and location \eqn{\delta} and discrimination \eqn{\beta} provides more information than just total scores.  With complete data, total scores and irt estimates are almost perfectly correlated.  However, the irt estimates provide much more information in the case of missing data.

The bounds parameter sets the lower and upper limits to the estimate.  This is relevant for the case of a subject who gives just the lowest score on every item, or just the top score on every item. Formerly (prior to 1.6.12) this was done by  estimating these taial scores  by finding the probability of missing every item taken, converting this to a quantile score based upon the normal distribution, and then assigning a z value equivalent to 1/2 of that quantile.  Similarly, if a person gets all the items they take correct, their score is defined as the quantile of the z equivalent to the probability of getting all of the items correct, and then moving up the distribution half way.  If these estimates exceed either the upper or lower bounds, they are adjusted to those boundaries. 

As of 1.6.9, the procedure is very different.  We now assume that all items are bounded with one passed item that is easier than all items given, and one failed item that is harder than any item given.  This produces much cleaner results.

There are several more elegant packages in R that provide Full Information Maximum Likeliood IRT based estimates. In particular, the MIRT package seems especially good.  The ltm package give equivalent estimates to MIRT for dichotomous data but produces unstable estimates for polytomous data and should be avoided.   

Although the scoreIrt estimates are are not FIML based they seem to correlated with  the MIRT estiamtes with values exceeding .99.  Indeed, based upon very limited simulations there are some small hints that the solutions match the true score estimates  slightly better than do the MIRT estimates.  \code{\link{scoreIrt}} seems to do a good job of recovering the basic structure.

The two wrapper functions \code{\link{scoreIrt.1pl}} and \code{\link{scoreIrt.2pl}} are very fast and are meant for scoring one or many scales at a time with a one factor model (\code{\link{scoreIrt.2pl}}) or just Rasch like scoring.  Just specify the scoring direction for a number of scales  (\code{\link{scoreIrt.1pl}}) or just items to score  for a number of scales \code{\link{scoreIrt.2pl}}.  \code{\link{scoreIrt.2pl}} will then apply \code{\link{irt.fa}} to the items for each scale separately, and then find the 2pl scores. 

The keys.list is a list of items to score for each scale.  Preceding the item name with a negative sign will reverse score that item (relevant for \code{\link{scoreIrt.1pl}}.  Alternatively, a keys matrix can be created using \code{\link{make.keys}}.  The keys matrix is a matrix of 1s, 0s, and -1s reflecting whether an item should be scored or not scored for a particular factor.  See \code{\link{scoreItems}} or \code{\link{make.keys}} for details.  The default case is to score all items with absolute discriminations > cut.


If one wants to score scales taking advantage of differences in item location but not do a full IRT analysis, then find the item difficulties from the raw data using \code{\link{irt.tau}} or combine this information with a scoring keys matrix (see \code{\link{scoreItems}} and \code{\link{make.keys}} and create quasi-IRT statistics using \code{\link{irt.stats.like}}.   This is the equivalent of doing a quasi-Rasch model, in that all items are assumed to be equally discriminating.  In this case, tau values may be found first (using \code{\link{irt.tau}} or just found before doing the scoring.  This is all done for you inside of \code{\link{scoreIrt.1pl}}. 

Such irt based scores are particularly useful if finding scales based upon massively missing data (e.g., the SAPA data sets).  Even without doing the full irt analysis, we can take into account different item difficulties. 

David Condon has added a very nice function to do 2PL analysis for a number of scales at one time.  \code{\link{scoreIrt.2pl}} takes the raw data file and a list of items to score for each of multiple scales.  These are then factored (currently just one factor for each scale) and the loadings and difficulties are used for scoring.  

There are conventionally two different metrics and models that are used.  The logistic metric and model and the normal metric and model. These are chosen using the mod parameter.

\code{\link{irt.se}} finds the standard errors for scores with a particular value.  These are based upon the information curves calculated by \code{\link{irt.fa}} and are not based upon the particular score of a particular subject.
}
\value{
 \item{scores}{A data frame of theta estimates, total scores based upon raw sums, and estimates of fit.}
 \item{tau}{Returned by irt.tau: A data frame of the tau values for an object of dichotomous or polytomous items.  Found without bothering to find the correlations.}
}
\references{
Kamata, Akihito and Bauer, Daniel J. (2008) A Note on the Relation Between Factor Analytic and Item Response Theory Models
Structural Equation Modeling, 15 (1) 136-153.
	
McDonald, Roderick P. (1999) Test theory: A unified treatment. L. Erlbaum Associates.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

}
\author{William Revelle, David Condon
}
\note{
It is very important to note that when using \code{\link{irt.fa}} to find the discriminations, to set the sort option to be FALSE. This is now the default.  Otherwise, the discriminations will not match the item order.  

Always under development.  Suggestions for improvement are most appreciated.

scoreIrt is just a wrapper to score.irt.poly and score.irt.2.  The previous version had score.irt which is now deprecated as I try to move to camelCase.

scoreIrt.2pl is a wrapper for irt.fa and scoreIrt.  It was originally developed by David Condon.
}


\seealso{
\code{\link{irt.fa}} for finding the parameters. For more conventional scoring algorithms see \code{\link{scoreItems}}. \code{\link{irt.responses}} will plot the empirical response patterns for the alternative response choices for  multiple choice items. For more conventional IRT estimations, see the ltm package.
}
\examples{
if(FALSE) {  #not run in the interest of time, but worth doing
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
scores <- scoreIrt(test,d9$items)
scores.df <- data.frame(scores,true=d9$theta) #combine the estimates with the true thetas.
pairs.panels(scores.df,pch=".",
main="Comparing IRT and classical with complete data") 
#now show how to do this with a quasi-Rasch model
tau <- irt.tau(d9$items)
scores.rasch <- scoreIrt(tau,d9$items,key=rep(1,9))
scores.dfr<- data.frame(scores.df,scores.rasch) #almost identical to 2PL model!
pairs.panels(scores.dfr)
#with all the data, why bother ?

#now delete some of the data
d9$items[1:333,1:3] <- NA
d9$items[334:666,4:6] <- NA
d9$items[667:1000,7:9] <- NA
scores <- scoreIrt(test,d9$items)
scores.df <- data.frame(scores,true=d9$theta) #combine the estimates with the true thetas.
pairs.panels(scores.df, pch=".",
main="Comparing IRT and classical with random missing data")
 #with missing data, the theta estimates are noticably better.
#now show how to do this with a quasi-Rasch model
tau <- irt.tau(d9$items)
scores.rasch <- scoreIrt(tau,d9$items,key=rep(1,9))
scores.dfr <- data.frame(scores.df,rasch = scores.rasch)
pairs.panels(scores.dfr)  #rasch is actually better!



v9 <- sim.irt(9,1000,-2.,2.,mod="normal") #dichotomous items
items <- v9$items
test <- irt.fa(items)
total <- rowSums(items)
ord <- order(total)
items <- items[ord,]


#now delete some of the data - note that they are ordered by score
items[1:333,5:9] <- NA
items[334:666,3:7] <- NA
items[667:1000,1:4] <- NA
items[990:995,1:9] <- NA   #the case of terrible data
items[996:998,] <- 0   #all wrong
items[999:1000] <- 1   #all right
scores <- scoreIrt(test,items)
unitweighted <- scoreIrt(items=items,keys=rep(1,9)) #each item has a discrimination of 1
#combine the estimates with the true thetas.
scores.df <- data.frame(v9$theta[ord],scores,unitweighted) 
   
colnames(scores.df) <- c("True theta","irt theta","total","fit","rasch","total","fit")
pairs.panels(scores.df,pch=".",main="Comparing IRT and classical with missing data") 
 #with missing data, the theta estimates are noticably better estimates 
 #of the generating theta than using the empirically derived factor loading weights

#now show the ability to score multiple scales using keys
ab.tau <- irt.tau(ability)  #first find the tau values
ab.keys <- make.keys(ability,list(g=1:16,reason=1:4,letter=5:8,matrix=9:12,rotate=13:16))
ab.scores <- scoreIrt(stats=ab.tau, items = ability, keys = ab.keys)

#and now do it for polytomous items 
bfi.tau <- irt.tau(bfi[1:25])
bfi.keys <- make.keys(bfi[1:25],list(agree=c(-1,2:5),conscientious=c(6:8,-9,-10),
extra=c(11:12,-13,-14,-15),neur = 16:20,open =c(21,-22,23,24,-25)))
bfi.scores <- scoreIrt(bfi.tau,bfi[1:25],bfi.keys)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}
