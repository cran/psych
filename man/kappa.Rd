\name{cohen.kappa}
\alias{wkappa}
\alias{cohen.kappa}

\title{Find Cohen's kappa and weighted kappa coefficients for correlation of two raters}
\description{Cohen's kappa (Cohen, 1960) and weighted kappa (Cohen, 1968) may be used to find the agreement of two raters when using nominal scores.  Light's kappa is just the average cohen.kappa if using more than 2 raters. 

weighted.kappa is (probability of observed matches - probability of expected matches)/(1 - probability of expected matches).  Kappa just considers the matches on the main diagonal.  Weighted kappa considers off diagonal elements as well.
}
\usage{
cohen.kappa(x, w=NULL,n.obs=NULL,alpha=.05,levels=NULL,w.exp=2)  
wkappa(x, w = NULL)    #deprectated
}

\arguments{
  \item{x}{Either a two by n data with categorical values from 1 to p or a p x p table.  If a data array, a table will be found.}
  \item{w}{A p x p matrix of weights.  If not specified, they are set to be 0 (on the diagonal) and (distance from diagonal) off the diagonal)^2. }
  \item{n.obs}{Number of observations (if input is a square matrix.}
  \item{alpha}{Probability level for confidence intervals}
  \item{levels}{Specify the levels if some levels of x or y are completely missing.  See Examples}
  \item{w.exp}{Exponent to apply to weights matrix -- see examples}
  
}
\details{When cateogorical judgments are made with two cateories, a measure of relationship is the phi coefficient.  However, some categorical judgments are made using more than two outcomes.  For example, two diagnosticians might be asked to categorize patients three ways (e.g., Personality disorder, Neurosis, Psychosis) or to categorize the stages of a disease.  Just as base rates affect observed cell frequencies in a two by two table, they need to be considered in the n-way table (Cohen, 1960). 

Kappa considers the matches on the main diagonal.  A penalty function (weight) may be applied to the off diagonal matches.  If the weights increase by the square of the distance from the diagonal, weighted kappa is similar to an Intra Class Correlation (\code{\link{ICC}}).

Derivations of weighted kappa are sometimes expressed in terms of similarities, and sometimes in terms of dissimilarities. In the latter case, the weights on the diagonal are 1 and the weights off the diagonal are less than one. In this  case, if the weights are 1 - squared distance from the diagonal / k, then the result is similar to the ICC (for any positive k).  

cohen.kappa may use either similarity weighting (diagonal = 0) or dissimilarity weighting (diagonal = 1) in order to match various published examples. 

The input may be a two column data.frame or matrix with columns representing the two judges and rows the subjects being rated. Alternatively, the input may be a square n x n matrix of counts or proportion of matches.  If proportions are used, it is necessary to specify the number of observations (n.obs) in order to correctly find the confidence intervals.

The confidence intervals are based upon the variance estimates discussed by Fleiss, Cohen, and Everitt who corrected the formulae of Cohen (1968) and Blashfield.

Some data sets will include data with numeric categories with some category values missing completely.  In the sense that kappa is a measure of category relationship, this should not matter.  But when finding weighted kappa, the number of categories weighted will be less than the number of categories potentially in the data.  This can be remedied by specifying the levels parameter.  This is a vector of the levels potentially in the data (even if some are missing).   See the examples.

If there are more than 2 raters, then the average of all raters is known as Light's kappa. (Conger, 1980). 

Following a request from Julius Pfadt, the weights matrix may be formed by squared distances (the default) or linear distances (w.exp=1)
}

\value{
  \item{kappa }{Unweighted kappa}
  \item{weighted.kappa }{The default weights are quadratric.}
  \item{var.kappa}{Variance of kappa}
  \item{var.weighted}{Variance of weighted kappa}
  \item{n.obs}{number of observations}
  \item{weight}{The weights used in the estimation of weighted kappa}
  \item{confid}{The alpha/2 confidence intervals for unweighted and weighted kappa}
  \item{plevel}{The alpha level used in determining the confidence limits}
}
\references{
Banerjee, M., Capozzoli, M., McSweeney, L and Sinha, D. (1999) Beyond Kappa: A review of interrater agreement measures The Canadian Journal of Statistics / La Revue Canadienne de Statistique, 27, 3-23

Cohen, J. (1960). A coefficient of agreement for nominal scales. Educational and Psychological Measurement, 20 37-46

Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for scaled disagreement or partial credit. Psychological Bulletin, 70, 213-220.

Conger, A. J. (1980) Integration and generalization of kappas for multiple raters, Psychological Bulletin,, 88, 322-328.

Fleiss, J. L., Cohen, J.  and Everitt, B.S. (1969) Large sample standard errors of kappa and weighted kappa. Psychological Bulletin, 72, 332-327.

Light, R. J. (12971) Measures of response agreement for qualitative data: Some generalizations and alternatives, Psychological Bulletin, 76, 365-377.

Zwick, R.  (1988) Another look at interrater agreement. Psychological Bulletin, 103, 374 - 378.

 }
 
\author{William Revelle }
\note{As is true of many R functions, there are alternatives in other packages.  The Kappa function in the vcd package estimates unweighted and weighted kappa and reports the variance of the estimate. The input is a square matrix. The ckappa and wkappa functions in the psy package take raw data matrices. The kappam.light function from the irr package finds Light's average kappa.

To avoid confusion with Kappa (from vcd) or the kappa function from base, the function was originally named wkappa. With additional features modified from psy::ckappa to allow input with a different number of categories, the function has been renamed cohen.kappa.  

Unfortunately, to make it more confusing, the weights described by Cohen are a function of the reciprocals of those discucssed by Fleiss and Cohen. The cohen.kappa function uses the appropriate formula for Cohen or Fleiss-Cohen weights. 

There are some cases where the large sample size approximation of Fleiss et al. will produce confidence intervals exceeding +/- 1.  Clearly, for these cases, the upper (or lower for negative values) should be set to 1.  Boot strap resampling shows the problem is that the values are not symmetric.  See the last (unrun) example.

It is also possible to have more than 2 raters.  In this case, cohen.kappa is reported for all pairs of raters (e.g. R1 and R2, R1 and R3,  ... R3 and R4).  To see the confidence intervals for these cohen.kappas, use the print command with the all=TRUE option. (See the exmaple of multiple raters.)

}


\examples{

#rating data (with thanks to Tim Bates)
rater1 = c(1,2,3,4,5,6,7,8,9) # rater one's ratings
rater2 = c(1,3,1,6,1,5,5,6,7) # rater two's ratings
cohen.kappa(x=cbind(rater1,rater2))

#data matrix taken from Cohen
cohen <- matrix(c(
0.44, 0.07, 0.09,
0.05, 0.20, 0.05,
0.01, 0.03, 0.06),ncol=3,byrow=TRUE)

#cohen.weights  weight differences
cohen.weights <- matrix(c(
0,1,3,
1,0,6,
3,6,0),ncol=3)


cohen.kappa(cohen,cohen.weights,n.obs=200)
#cohen reports .492 and .348 

#another set of weights
#what if the weights are non-symmetric
wc <- matrix(c(
0,1,4,
1,0,6,
2,2,0),ncol=3,byrow=TRUE)
cohen.kappa(cohen,wc)  
#Cohen reports kw = .353

cohen.kappa(cohen,n.obs=200)  #this uses the squared weights

fleiss.cohen <- 1 - cohen.weights/9
cohen.kappa(cohen,fleiss.cohen,n.obs=200)

#however, Fleiss, Cohen and Everitt weight similarities
fleiss <- matrix(c(
106, 10,4,
22,28, 10,
2, 12,  6),ncol=3,byrow=TRUE)

#Fleiss weights the similarities
weights <- matrix(c(
 1.0000, 0.0000, 0.4444,
 0.0000, 1.0000, 0.6667,
 0.4444, 0.6667, 1.0000),ncol=3)
 
 cohen.kappa(fleiss,weights,n.obs=200)
 
 #another example is comparing the scores of two sets of twins
 #data may be a 2 column matrix
 #compare weighted and unweighted
 #also look at the ICC for this data set.
 twins <- matrix(c(
    1, 2, 
    2, 3,
    3, 4,
    5, 6,
    6, 7), ncol=2,byrow=TRUE)
  cohen.kappa(twins)
  
#data may be explicitly categorical
x <- c("red","yellow","blue","red")
y <- c("red",  "blue", "blue" ,"red") 
xy.df <- data.frame(x,y)
ck <- cohen.kappa(xy.df)
ck
ck$agree

#Example for specifying levels
#The problem of missing categories (from Amy Finnegan)
#We need to specify all the categories possible using the levels option
numbers <- data.frame(rater1=c(6,3,7,8,7),
                      rater2=c(6,1,8,5,10))
cohen.kappa(numbers)  #compare with the next analysis
cohen.kappa(numbers,levels=1:10)  #specify the number of levels 
              #   these leads to slightly higher weighted kappa
  
#finally, input can be a data.frame of ratings from more than two raters
ratings <- matrix(rep(1:5,4),ncol=4)
ratings[1,2] <- ratings[2,3] <- ratings[3,4] <- NA
ratings[2,1] <- ratings[3,2] <- ratings[4,3] <- 1
ck <- cohen.kappa(ratings)
ck  #just show the raw and weighted kappas
print(ck, all=TRUE)  #show the confidence intervals as well

 
 #In the case of confidence intervals being artificially truncated to +/- 1, it is 
 #helpful to compare the results of a boot strap resample
 #ck.boot <-function(x,s=1:nrow(x)) {cohen.kappa(x[s,])$kappa}
 #library(boot)
 #ckb <- boot(x,ck.boot,R=1000)
 #hist(ckb$t)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

