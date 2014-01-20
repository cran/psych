\name{score.irt}
\alias{score.irt}
\alias{score.irt.poly}
\alias{score.irt.2}
\alias{irt.stats.like}
\alias{irt.tau}
\title{Find Item Response Theory (IRT) based scores for dichotomous or polytomous items}

\description{\code{\link{irt.fa}} finds Item Response Theory (IRT) parameters through factor analysis of the tetrachoric or polychoric correlations of dichtomous or polytomous items. \code{\link{score.irt}} uses these parameter estimates of discrimination and location to find IRT based scores for the responses. As many factors as found for the correlation matrix will be scored. 
}
\usage{
score.irt(stats=NULL, items, keys=NULL,cut = 0.3,bounds=c(-5,5),mod="logistic") 
 #the higher order call just calls one of the next two
  #for dichotomous items 
score.irt.2(stats, items,keys=NULL,cut = 0.3,bounds=c(-5,5),mod="logistic") 
  #for polytomous items
score.irt.poly(stats, items, keys=NULL, cut = 0.3,bounds=c(-5,5))
    #to create irt like statistics for plotting
irt.stats.like(items,stats,keys=NULL,cut=.3)
   
irt.tau(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stats}{Output from irt.fa is used for parameter estimates of location and discrimination. Stats may also be the  output from a normal factor analysis (fa) }
  \item{items}{The raw data, may be either dichotomous or polytomous.}
  \item{keys}{A keys matrix of which items should be scored for each factor}
  \item{cut}{Only items with discrimination values > cut will be used for scoring.}
  \item{x}{The raw data to be used to find the tau parameter in irt.tau}
  \item{bounds}{The lower and upper estimates for the fitting function}
  \item{mod}{Should a logistic or normal model be used in estimating the scores?}
}
\details{Although there are more elegant ways of finding subject scores given a set of item locations (difficulties) and discriminations, simply finding that value of theta \eqn{\theta} that best fits the equation \eqn{P(x|\theta) = 1/(1+exp(\beta(\delta - \theta) )} for a score vector X, and location \eqn{\delta} and discrimination \eqn{\beta} provides more information than just total scores.  With complete data, total scores and irt estimates are almost perfectly correlated.  However, the irt estimates provide much more information in the case of missing data.

The bounds parameter sets the lower and upper limits to the estimate.  This is relevant for the case of a subject who gives just the lowest score on every item, or just the top score on every item. In this case, the scores are estimated by finding the probability of missing every item taken, converting this to a quantile score based upon the normal distribution, and then assigning a z value equivalent to 1/2 of that quantile.  Similarly, if a person gets all the items they take correct, their score is defined as the quantile of the z equivalent to the probability of getting all of the items correct, and then moving up the distribution half way.  If these estimates exceed either the upper or lower bounds, they are adjusted to those boundaries. 

There are several more elegant packages in R that provide Full Information Maximum Likeliood IRT based estimates.  The estimates from score.irt do not do so.  However, the score.irt seems to do a good job of recovering the basic structure.

The keys matrix is a matrix of 1s, 0s, and -1s reflecting whether an item should be scored or not scored for a particular factor.  See \code{\link{score.items}} or \code{\link{make.keys}} for details.  The default case is to score all items with absolute discriminations > cut.

If one wants to score scales taking advantage of differences in item location but not do a full irt analysis, then find the item difficulties from the raw data using \code{\link{irt.tau}} or combine this information with a scoring keys matrix (see \code{\link{score.items}} and code{\link{make.keys}} and create quasi-irt statistics using \code{\link{irt.stats.like}}. 

There are conventionally two different metrics and models that are used.  The logistic metric and model and the normal metric and model. These are chosen using the mod parameter.
}
\value{
 \item{scores}{A data frame of theta estimates, total scores based upon raw sums, and estimates of fit.}
}
\references{
Kamata, Akihito and Bauer, Daniel J. (2008) A Note on the Relation Between Factor Analytic and Item Response Theory Models
Structural Equation Modeling, 15 (1) 136-153.
	
McDonald, Roderick P. (1999) Test theory: A unified treatment. L. Erlbaum Associates.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 

}
\author{William Revelle
}
\note{
Still under development.  Suggestions for improvement are most appreciated.

score.irt is just a wrapper to score.irt.poly and score.irt.2
}


\seealso{
\code{\link{irt.fa}} for finding the parameters. For more conventional scoring algorithms see \code{\link{score.items}}. \code{\link{irt.responses}} will plot the empirical response patterns for the alternative response choices for  multiple choice items. For more conventional IRT estimations, see the ltm package.
}
\examples{
if(FALSE) {  #not run in the interest of time, but worth doing
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
scores <- score.irt(test,d9$items)
scores.df <- data.frame(scores,true=d9$theta) #combine the estimates with the true thetas.
pairs.panels(scores.df,pch=".",
main="Comparing IRT and classical with complete data") 
#with all the data, why bother ?

#now delete some of the data
d9$items[1:333,1:3] <- NA
d9$items[334:666,4:6] <- NA
d9$items[667:1000,7:9] <- NA
scores <- score.irt(test,d9$items)
scores.df <- data.frame(scores,true=d9$theta) #combine the estimates with the true thetas.
pairs.panels(scores.df, pch=".",
main="Comparing IRT and classical with random missing data")
 #with missing data, the theta estimates are noticably better.
}

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
scores <- score.irt(test,items)
unitweighted <- score.irt(items=items,keys=rep(1,9)) #each item has a discrimination of 1
#combine the estimates with the true thetas.
scores.df <- data.frame(v9$theta[ord],scores,unitweighted) 
   
colnames(scores.df) <- c("True theta","irt theta","total","fit","rasch","total","fit")
pairs.panels(scores.df,pch=".",main="Comparing IRT and classical with missing data") 
 #with missing data, the theta estimates are noticably better estimates 
 #of the generating theta than calling them all equal
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}
