\name{guttman}
\alias{tenberge}
\alias{glb}
\alias{guttman} 

\title{Alternative estimates of test reliabiity }
\description{Eight alternative estimates of test reliability include the six discussed by Guttman (1945), four discussed by ten Berge and Zergers (1978) (\eqn{\mu_0 \dots \mu_3)} as well as \eqn{\beta} (the worst split half, Revelle, 1979),  the glb (greatest lowest bound) discussed by Bentler and Woodward (1980), and \eqn{\omega_h} and \eqn{\omega_t} (McDonald, 1999; Zinbarg et al., 2005). 
}
\usage{
guttman(r,key=NULL,digits=2) 
tenberge(r,digits=2)
glb(r,key=NULL,digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{A correlation matrix or raw data matrix.}
  \item{key}{a vector of -1, 0, 1 to select or reverse items}
  \item{digits}{How many digits of accuracy in the output?}
}
\details{Surprisingly, 104 years after Spearman (1904) introduced the concept of reliability to psychologists, there are still multiple approaches for measuring it. Although very popular, Cronbach's \eqn{\alpha}  (1951) underestimates the reliability of a test and over estimates the first factor saturation. The guttman function includes the six estimates discussed by Guttman (1945), four of ten Berge and Zergers (1978), as well as Revelle's \eqn{\beta} (1979) using \code{\link{ICLUST}}. The companion function, \code{\link{omega}} calculates omega hierarchical (\eqn{\omega_h})  and omega total (\eqn{\omega_t}). 

\deqn{
\lambda_1 = 1 - \frac{tr(\vec{V_x})}{V_x} = \frac{V_x - tr(\vec{V}_x)}{V_x}
}
The second bound, \eqn{\lambda_2} replaces the diagonal with a function of the square root of the  sums of squares of the off diagonal elements.  Let \eqn{C_2 = \vec{1}( \vec{V}-diag(\vec{V}))^2 \vec{1}' }, then 
\deqn{
\lambda_2 = \lambda_1 + \frac{\sqrt{\frac{n}{n-1}C_2}}{V_x} = \frac{V_x - tr(\vec{V}_x) + \sqrt{\frac{n}{n-1}C_2} }{V_x}.
}
Effectively, this is replacing the diagonal with  n * the square root of the average squared off diagonal element.  

Guttman's 3rd lower bound, \eqn{\lambda_3}, also  modifies \eqn{\lambda_1} and estimates the true variance of each item as the average covariance between items and is, of course, the same as Cronbach's \eqn{\alpha}. 
\deqn{
\lambda_3 = \lambda_1 + \frac{\frac{V_X - tr(\vec{V}_X)}{n (n-1)}}{V_X} = \frac{n \lambda_1}{n-1} = \frac{n}{n-1}\Bigl(1 - \frac{tr(\vec{V})_x}{V_x}\Bigr) = \frac{n}{n-1} \frac{V_x - tr(\vec{V}_x)}{V_x} = \alpha
}
This is just replacing the diagonal elements with the average off diagonal elements.  \eqn{\lambda_2 \geq \lambda_3} with  \eqn{\lambda_2 > \lambda_3} if the covariances are not identical.

 \eqn{\lambda_3} and \eqn{\lambda_2} are both corrections to \eqn{\lambda_1} and this correction may  be  generalized  as an infinite set of successive improvements. (Ten Berge and Zegers, 1978) 
\deqn{
\mu_r = \frac{1}{V_x} \bigl( p_o + (p_1 + (p_2 + \dots (p_{r-1} +( p_r)^{1/2})^{1/2} \dots )^{1/2})^{1/2}    \bigr), r = 0, 1, 2, \dots
}
where 
\deqn{
p_h = \sum_{i\ne j}\sigma_{ij}^{2h}, h = 0, 1, 2, \dots r-1
}
and
\deqn{
p_h = \frac{n}{n-1}\sigma_{ij}^{2h}, h = r 
} 
tenberge \& Zegers (1978).  Clearly \eqn{\mu_0 = \lambda_3 = \alpha} and \eqn{ \mu_1 = \lambda_2}.  \eqn{\mu_r \geq \mu_{r-1} \geq \dots \mu_1 \geq \mu_0}, although the series does not improve much after the first two steps.

Guttman's fourth lower bound, \eqn{\lambda_4} was originally proposed as any spit half reliability  but has been interpreted as the greatest split half reliability. If \eqn{\vec{X}} is split into  two parts, \eqn{\vec{X}_a} and \eqn{\vec{X}_b}, with correlation \eqn{r_{ab}} then
\deqn{
\lambda_4 = 2\Bigl(1 - \frac{V_{X_a} + V_{X_b}}{V_X} \Bigr) =  \frac{4 r_{ab}}{V_x} = \frac{4 r_{ab}}{V_{X_a} + V_{X_b}+ 2r_{ab}V_{X_a}  V_{X_b}}
}
which is just the normal split half reliability, but in this case, of the most similar splits.

\eqn{\lambda_5}, Guttman's fifth lower bound, replaces the diagonal values with twice the square root of the maximum (across items) of the sums of squared interitem covariances
\deqn{
\lambda_5 =  \lambda_1 + \frac{2 \sqrt{\bar{C_2}}}{V_X}.
}
Although superior to \eqn{\lambda_1}, \eqn{\lambda_5} underestimates the correction to the diagonal.  A better estimate would be analogous to the correction used in \eqn{\lambda_3}:
\deqn{
\lambda_{5+} =  \lambda_1 + \frac{n}{n-1}\frac{2 \sqrt{\bar{C_2}}}{V_X}.
}

Guttman's final bound considers the amount of variance in each item that can be accounted for the linear regression of all of the other items (the squared multiple correlation or smc), or more precisely, the variance of the errors, \eqn{e_j^2},  and is
\deqn{
\lambda_6 = 1 - \frac{\sum e_j^2}{V_x} = 1 - \frac{\sum(1-r_{smc}^2)}{V_x}
}


Guttman's \eqn{\lambda_4} is the greatest split half reliability.  This is found by combining the output from three different approaches, and seems to work for all test cases yet tried. Lambda 4 is reported as the max of these three algorithms.


The algorithms are 

a) Do an ICLUST of the reversed correlation matrix.  ICLUST normally forms the most distinct clusters.  By reversing the correlations, it will tend to find the most related cluster.  Truly a weird approach but tends to work.

b) Alternatively, a kmeans clustering of the correlations (with the diagonal replaced with 0 to make pseudo distances) can produce 2 similar clusters.

c) Clusters identified by assigning items to two clusters based upon their order on the first principal factor.  (Highest to cluster 1, next 2 to cluster 2, etc.)


}

\value{
\item{beta}{The normal beta estimate of cluster similarity from ICLUST.  This is an estimate of the general factor saturation.}
\item{tenberge$mu1}{tenBerge mu 1 is functionally alpha}
 \item{tenberge$mu2}{one of the sequence of estimates mu1 ... mu3}
\item{beta.factor}{For experimental purposes, what is the split half based upon the two factor solution?}
\item{glb.IC}{Greatest split half based upon ICLUST of reversed correlations}
\item{glb.Km}{Greatest split half based upon a kmeans clustering.}
\item{glb.Fa}{Greatest split half based upon the items assigned by factor analysis.}
\item{glb.max}{max of the above estimates}
\item{keys}{scoring keys from each of the alternative methods of forming best splits}
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

\seealso{  \code{\link{omega}}, \code{\link{ICLUST}}, }
\examples{
data(attitude)
glb(attitude)
guttman(attitude)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

