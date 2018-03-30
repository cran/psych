\name{splitHalf}
\alias{splitHalf}
\alias{guttman}
\alias{tenberge}
\alias{glb}
\alias{glb.fa}
 

\title{Alternative estimates of test reliabiity }
\description{Eight alternative estimates of test reliability include the six discussed by Guttman (1945), four discussed by ten Berge and Zergers (1978) (\eqn{\mu_0 \dots \mu_3)} as well as \eqn{\beta} (the worst split half, Revelle, 1979),  the glb (greatest lowest bound) discussed by Bentler and Woodward (1980), and \eqn{\omega_h} and \eqn{\omega_t} (McDonald, 1999; Zinbarg et al., 2005). Greatest and lowest split-half values are found by brute force or sampling. 
}
\usage{
splitHalf(r,raw=FALSE,brute=FALSE,n.sample=10000,covar=FALSE,check.keys=TRUE,
           key=NULL,ci=.05,use="pairwise")
guttman(r,key=NULL) 
tenberge(r)
glb(r,key=NULL)
glb.fa(r,key=NULL)
}

\arguments{
  \item{r}{A correlation or covariance matrix or raw data matrix.}
  \item{raw}{return a vector of split half reliabilities}
  \item{brute}{Use brute force to try all combinations of n take n/2. }
  \item{n.sample}{if brute is false, how many samples of split halves should be tried?}
  \item{covar}{Should the covariances or correlations be used for reliability calculations}
  \item{check.keys}{If TRUE, any item with a negative loading on the first factor will be flipped in sign}
 \item{key}{a vector of -1, 0, 1 to select or reverse key items.  See if the key vector is less than the number of variables, then item numbers to be reverse can be specified.}
 \item{use}{Should we find the correlations using "pairwise" or "complete" (see ?cor)}
 \item{ci}{The alpha level to use for the confidence intervals of the split half estimates}
}
\details{Surprisingly, more than a century after Spearman (1904) introduced the concept of reliability to psychologists, there are still multiple approaches for measuring it. Although very popular, Cronbach's \eqn{\alpha}  (1951) underestimates the reliability of a test and over estimates the first factor saturation. Using \code{\link{splitHalf}} for tests with 16 or fewer items, all possible splits may be found fairly easily.  For tests with 17 or more items, n.sample splits are randomly found. Thus, for 16 or fewer items, the upper and lower bounds are precise.  For 17 or more items, they are close but will probably slightly underestimate the highest and overestimate the lowest reliabilities.  

 The guttman function includes the six estimates discussed by Guttman (1945), four of ten Berge and Zergers (1978), as well as Revelle's \eqn{\beta} (1979) using \code{\link{splitHalf}}. The companion function, \code{\link{omega}} calculates omega hierarchical (\eqn{\omega_h})  and omega total (\eqn{\omega_t}). 

Guttman's first estimate \eqn{\lambda_1} assumes that all the variance of an item is error:
\deqn{
\lambda_1 = 1 - \frac{tr(\vec{V_x})}{V_x} = \frac{V_x - tr(\vec{V}_x)}{V_x}
}{lambda 1= 1-tr(Vx)/Vx}
This is a clear underestimate.

The second bound, \eqn{\lambda_2}{\lambda_2}, replaces the diagonal with a function of the square root of the  sums of squares of the off diagonal elements.  Let \eqn{C_2 = \vec{1}( \vec{V}-diag(\vec{V})^2 \vec{1}' }, then 
\deqn{
\lambda_2 = \lambda_1 + \frac{\sqrt{\frac{n}{n-1}C_2}}{V_x} = \frac{V_x - tr(\vec{V}_x) + \sqrt{\frac{n}{n-1}C_2} }{V_x}}{\lambda_2= \lambda_1 + sqrt(n *(n-1)C_2)/V_x)}
Effectively, this is replacing the diagonal with  n * the square root of the average squared off diagonal element.  

Guttman's 3rd lower bound, \eqn{\lambda_3}, also  modifies \eqn{\lambda_1} and estimates the true variance of each item as the average covariance between items and is, of course, the same as Cronbach's \eqn{\alpha}. 
\deqn{
\lambda_3 = \lambda_1 + \frac{\frac{V_X - tr(\vec{V}_X)}{n (n-1)}}{V_X} = \frac{n \lambda_1}{n-1} = \frac{n}{n-1}\Bigl(1 - \frac{tr(\vec{V})_x}{V_x}\Bigr) = \frac{n}{n-1} \frac{V_x - tr(\vec{V}_x)}{V_x} = \alpha
}{\lambda 3 = ((n)/(n-1))(1-tr(Vx)/(Vx)  = ((n)/(n-1))(Vx-tr(Vx)/Vx  = \alpha}

This is just replacing the diagonal elements with the average off diagonal elements.  \eqn{\lambda_2 \geq \lambda_3}{\lambda_2 \ge \lambda_3} with  \eqn{\lambda_2 > \lambda_3} if the covariances are not identical.

 \eqn{\lambda_3} and \eqn{\lambda_2} are both corrections to \eqn{\lambda_1} and this correction may  be  generalized  as an infinite set of successive improvements. (Ten Berge and Zegers, 1978) 
\deqn{
\mu_r = \frac{1}{V_x} \bigl( p_o + (p_1 + (p_2 + \dots (p_{r-1} +( p_r)^{1/2})^{1/2} \dots )^{1/2})^{1/2}    \bigr), r = 0, 1, 2, \dots
}{(1/(Vx))(po + p1 = (p2 + ... (pr1) + pr^.5 )^.5^ ... .5)}
where 
\deqn{
p_h = \sum_{i\ne j}\sigma_{ij}^{2h}, h = 0, 1, 2, \dots r-1
}{p_h = sum(\sigma^2h,  h = 0, 1, 2, ... r-1  }
and
\deqn{
p_h = \frac{n}{n-1}\sigma_{ij}^{2h}, h = r 
}{p_h = n/((n-1) \sigma^2h) }
tenberge and Zegers (1978).  Clearly \eqn{\mu_0 = \lambda_3 = \alpha} and \eqn{ \mu_1 = \lambda_2}.  \eqn{\mu_r \geq \mu_{r-1} \geq \dots \mu_1 \geq \mu_0}{\mu_r \ge \mu_{r-1} \ge \dots \mu_1 \ge \mu_0}, although the series does not improve much after the first two steps.

Guttman's fourth lower bound, \eqn{\lambda_4} was originally proposed as any spit half reliability  but has been interpreted as the greatest split half reliability. If \eqn{\vec{X}} is split into  two parts, \eqn{\vec{X}_a} and \eqn{\vec{X}_b}, with correlation \eqn{r_{ab}} then
\deqn{
\lambda_4 = 2\Bigl(1 - \frac{V_{X_a} + V_{X_b}}{V_X} \Bigr) =  \frac{4 r_{ab}}{V_x} = \frac{4 r_{ab}}{V_{X_a} + V_{X_b}+ 2r_{ab}V_{X_a}  V_{X_b}}
}{\lambda 4 = 4rab/(Va + Vb + 2rabVaVb)}
which is just the normal split half reliability, but in this case, of the most similar splits. For 16 or fewer items, this is found by trying all possible splits.  For 17 or more items, this is estimated by taking n.sample random splits. 


\eqn{\lambda_5}, Guttman's fifth lower bound, replaces the diagonal values with twice the square root of the maximum (across items) of the sums of squared interitem covariances
\deqn{
\lambda_5 =  \lambda_1 + \frac{2 \sqrt{\bar{C_2}}}{V_X}.
}{\lambda_5 = \lambda_1 +2/sqrt(average(C_2)/V_X.) }
Although superior to \eqn{\lambda_1}, \eqn{\lambda_5} underestimates the correction to the diagonal.  A better estimate would be analogous to the correction used in \eqn{\lambda_3}:
\deqn{
\lambda_{5+} =  \lambda_1 + \frac{n}{n-1}\frac{2 \sqrt{\bar{C_2}}}{V_X}.
}{\lambda 5+ = \lambda 1 + ((n/(n-1))2/sqrt(av covariance 12)/Vx}

\eqn{\lambda_6},Guttman's final bound considers the amount of variance in each item that can be accounted for the linear regression of all of the other items (the squared multiple correlation or smc), or more precisely, the variance of the errors, \eqn{e_j^2},  and is
\deqn{\lambda_6 = 1 - \frac{\sum e_j^2}{V_x} = 1 - \frac{\sum(1-r_{smc}^2)}{V_x}
}{\lambda 6 = 1 - sum(e^2)/Vx = 1-sum(1-r^2(smc))/Vx}.

The smc is found from all the items.  A modification to Guttman \eqn{\lambda_6}, \eqn{\lambda_6*} reported by the \code{\link{score.items}} function is to find the smc from the entire pool of items given, not just the items on the selected scale.  



Guttman's \eqn{\lambda_4} is the greatest split half reliability.  Although originally found here by combining the output from three different approaches,this has now been replaced by using \code{\link{splitHalf}} to find the maximum value by brute force (for 16 or fewer items) or by taking a substantial number of random splits.

The algorithms that had been tried before included:

a) Do an ICLUST of the reversed correlation matrix.  ICLUST normally forms the most distinct clusters.  By reversing the correlations, it will tend to find the most related clusters.  Truly a weird approach but tends to work.

b) Alternatively, a kmeans clustering of the correlations (with the diagonal replaced with 0 to make pseudo distances) can produce 2 similar clusters.

c) Clusters identified by assigning items to two clusters based upon their order on the first principal factor.  (Highest to cluster 1, next 2 to cluster 2, etc.)

These three procedures will produce keys vectors for assigning items to the two splits.  The maximum split half reliability is found by taking the maximum of these three approaches.  This is not elegant but is fast.

The brute force and the sampling procedures seem to provide more stable and larger estimates. 

Yet another procedure, implemented in \code{\link{splitHalf}} is actually form all possible (for n items <= 16) or sample 10,000 (or more) split halfs corrected for test length.  This function returns the best and worst splits as item keys that can be used for scoring purposes, if desired.  Can do up to 24 items in reasonable time, but gets much slower for more than about 24 items.  

There are three greatest lower bound functions.  One, glb finds the greatest split half reliability, \eqn{\lambda_4}. This considers the test as set of items and examines how best to partition the items into splits. The other two, \code{\link{glb.fa}} and \code{\link{glb.algebraic}}, are alternative ways of weighting the diagonal of the matrix. 

\code{\link{glb.fa}} estimates the communalities of the variables from a factor model where the number of factors is the number with positive eigen values.  Then reliability is found by 
\deqn{
glb = 1 - \frac{\sum e_j^2}{V_x} = 1 - \frac{\sum(1- h^2)}{V_x}
}{glb = 1 - sum(e^2)/Vx = 1-sum(1-h^2)/Vx}

This estimate will differ slightly from that found by  \code{\link{glb.algebraic}}, written by Andreas Moeltner which uses calls to  csdp in the Rcsdp package. His algorithm, which more closely matches the description of the glb by Jackson and Woodhouse, seems to have a positive bias (i.e., will over estimate the reliability of some items; they are said to be = 1) for small sample sizes.  More exploration of these two algorithms is underway. 

Compared to \code{\link{glb.algebraic}}, \code{\link{glb.fa}} seems to have less (positive) bias for smallish sample sizes (n < 500) but larger for large (> 1000) sample sizes. This interacts with the number of variables so that equal bias sample size differs as a function of the number of variables.  The differences are, however small. As samples sizes grow,  \code{\link{glb.algebraic}} seems to converge on the population value while glb.fa has a positive bias. 
}

\value{
\item{beta}{The worst split half reliability. This is an estimate of the general factor saturation.}
\item{maxrb}{The maximimum split half reliability.  This is Guttman's lambda 4}
\item{alpha}{Also known as Guttman's Lambda 3}
\item{ci}{The 2.5\%, 50\%, and 97.5\%  values of the raw or sampled split half}
\item{tenberge$mu1}{tenBerge mu 1 is functionally alpha}
 \item{tenberge$mu2}{one of the sequence of estimates mu1 ... mu3}
\item{glb}{glb found from factor analysis}

}
\references{
Cronbach, L.J. (1951) Coefficient alpha and the internal strucuture of tests.  Psychometrika, 16, 297-334.

Guttman, L. (1945). A basis for analyzing test-retest reliability. Psychometrika, 10 (4), 255-282. 

Revelle, W. (1979). Hierarchical cluster-analysis and the internal structure of tests. Multivariate Behavioral Research, 14 (1), 57-74. 

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma. Psychometrika, 2009. 

Ten Berge, J. M. F., & Zegers, F. E. (1978). A series of lower bounds to the reliability of a test. Psychometrika, 43 (4), 575-579. 

Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's \eqn{\alpha} , Revelle's \eqn{\beta} , and McDonald's \eqn{\omega_h} ): Their relations with each other and two alternative conceptualizations of reliability. Psychometrika, 70 (1), 123-133.

}
\author{ William Revelle }

\seealso{ \code{\link{alpha}}, \code{\link{omega}}, \code{\link{ICLUST}}, \code{\link{glb.algebraic}}
}
\examples{
data(attitude)
splitHalf(attitude)
splitHalf(attitude,covar=TRUE) #do it on the covariances
glb(attitude)
glb.fa(attitude)
if(require(Rcsdp)) {glb.algebraic(cor(attitude)) }
guttman(attitude)

#to show the histogram of all possible splits for the ability test
#sp <- splitHalf(ability,raw=TRUE)  #this saves the results
#hist(sp$raw,breaks=101,ylab="SplitHalf reliability",main="SplitHalf 
#    reliabilities of a test with 16 ability items")
sp <- splitHalf(bfi[1:10],key=c(1,9,10))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

