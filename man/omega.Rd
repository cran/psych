\name{omega}
\alias{omega}
\title{  Calculate McDonald's  omega estimates of general and total factor saturation }
\description{McDonald has proposed coefficient omega as an estimate of the general factor saturation of a test.  One way to find omega is to do a factor analysis of the original data set, rotate the factors obliquely, do a Schmid Leiman transformation, and then find omega. This function estimates omega as suggested by McDonald by using hierarchical factor analysis (following Jensen).
}
\usage{
omega(m, nfactors=3, fm = "minres",key = NULL, flip=TRUE, digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A correlation matrix or a data.frame/matrix of data }
  \item{nfactors}{Number of factors believed to be group factors}
  \item{fm}{factor method  fm="pa" for principal axes, fm="minres" for a minimum residual (OLS) solution, fm="pc" for principal components, fm="mle" for maximum likelihood.}
  \item{key}{a vector of +/- 1s to specify the direction of scoring of items.  The default is to assume all items are positively keyed, but if some items are reversed scored, then key should be specified.}
  \item{flip}{If flip is TRUE, then items are automatically flipped to have positive correlations on the general factor. Items that have been reversed are shown with a - sign.}
  \item{digits}{if specified, round the output to digits}
  \item{title}{Title for this analysis}
  \item{sl}{If plotting the results, should the Schmid Leiman solution be shown or should the hierarchical solution be shown? (default sl=TRUE)}
  \item{labels}{If plotting, what labels should be applied to the variables? If not specified, will default to the column names.}
  \item{plot}{plot=TRUE (default) calls omega.graph, plot =FALSE does not}
  \item{n.obs}{Number of observations - used for goodness of fit statistic}
  \item{rotate}{What rotation to apply? The default is oblimin, the alternatives include simplimax, Promax,  cluster and target. target will rotate to an optional keys matrix (See \code{\link{target.rot}})}
  
  \item{...}{Allows additional parameters to be passed through to the factor routines.  See in particular the option parameter in  \code{\link{schmid}} for treating the case of two group factors.}
}

\details{``Many scales are assumed by their developers and users to be primarily a measure of one latent variable. When it is also assumed that the scale conforms to the effect indicator model of measurement (as is almost always the case in psychological assessment), it is important to support such an 
interpretation with evidence regarding the internal structure of that scale. In particular, it is important to examine two related properties pertaining to the internal structure of such a scale. The first property relates to whether all the indicators forming the scale measure a latent variable in common. 

The second internal structural property pertains to the proportion of variance in the scale scores (derived from summing or averaging the indicators) accounted for by this latent variable that is common to all the indicators (Cronbach, 1951; McDonald, 1999; Revelle, 1979). That is, if an effect indicator scale is primarily a measure of one latent variable common to all the indicators forming the scale, then that latent variable should account for the majority of the variance in the scale scores. Put differently, this variance ratio provides important information about the sampling fluctuations when estimating individuals' standing on a latent variable common to all the indicators arising from the sampling of indicators (i.e., when dealing with either Type 2 or Type 12 sampling, to use the terminology of Lord, 1956). That is, this variance proportion can be interpreted as the square of the correlation between the scale score and the latent variable common to all the indicators in the infinite universe of indicators of which the scale indicators are a subset. Put yet another way, this variance ratio is important both as reliability and a validity coefficient. This is a reliability issue as the larger this variance ratio is, the more accurately one can predict an individual's relative standing on the latent variable common to all the scale's indicators based on his or her 
observed scale score. At the same time, this variance ratio also bears on the construct validity of the scale given that construct validity encompasses the internal structure of a scale." (Zinbarg, Yovel, Revelle, and McDonald, 2006).

McDonald has proposed coefficient omega (hierarchical (\eqn{\omega_h}) as an estimate of the general factor saturation of a test.  Zinbarg, Revelle, Yovel and Li (2005) 
\url{http://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf} compare McDonald's \eqn{\omega_h} to Cronbach's \eqn{\alpha} and Revelle's \eqn{\beta}  They conclude that \eqn{\omega_h} is the best estimate. (See also Zinbarg et al., 2006 and Revelle and Zinbarg (2009)).   

One way to find \eqn{\omega_h} is to do a factor analysis of the original data set, rotate the factors obliquely, factor that correlation matrix, do a   Schmid-Leiman (\link{schmid}) transformation to find general factor loadings, and then find \eqn{\omega_h}.  Here we present code to do that.  

\eqn{\omega_h} differs as a function of how the factors are estimated.  Three options are available, pc="pa"  does a principle axes factor analysis (\link{factor.pa}), pc="mle" uses the factanal function, and pc="pc" does a principal components analysis (\link{principal}).  

For ability items, it is typically the case that all items will have positive loadings on the general factor.  However, for non-cognitive items it is frequently the case that some items are to be scored positively, and some negatively.  Although probably better to specify which directions the items are to be scored by specifying a key vector, if flip =TRUE (the default), items will be reversed so that they have positive loadings on the general factor.  The keys are reported so that scores can be found using the \code{\link{score.items}} function.  Arbitrarly reversing items this way can overestimate the general factor. (See the example with a simulated circumplex).


Output from omega is shown graphically using the  \code{\link{omega.graph}} function.  This requires Rgraphviz to be installed.  If Rgraphviz is not available, select plot=FALSE. 

\eqn{\beta}, an alternative to \eqn{\omega}, is defined as the worst split half reliability.  It can be estimated by using \link{ICLUST} (a hierarchical clustering algorithm originally developed for main frames and written in Fortran and that is now available in R.  (For a very complimentary review of why the ICLUST algorithm is useful in scale construction, see Cooksey and Soutar, 2005). 

The \code{\link{omega}} function uses exploratory factor analysis to estimate the \eqn{\omega_h}{omega_h} coefficient.  It is important to remember that  ``A recommendation that should be heeded, regardless of the method chosen to estimate \eqn{\omega_h}, is to always examine the pattern of the estimated general factor loadings prior to estimating \eqn{\omega_h}. Such an examination constitutes an informal test of the assumption that there is a latent variable common to all of the scale's indicators that can be conducted even in the context of EFA. If the loadings were salient for only a relatively small subset of the indicators, this would suggest that there is no true general factor underlying the covariance matrix. Just such an informal assumption test would have afforded a great deal of protection against the possibility of misinterpreting the misleading \eqn{\omega_h} estimates occasionally produced in the simulations reported here." (Zinbarg et al., 2006, p 137).

A simple demonstration of the problem of an omega estimate reflecting just one of two group factors can be found in the last example.

Although omega is uniquely defined only for cases where 3 or more subfactors are extracted, it is sometimes desired to have a two factor solution.  By default this is done by forcing the schmid extraction to treat the two subfactors as having equal loadings.  

There are three possible options for this condition: setting the general factor loadings between the two lower order factors to be "equal" which will be the sqrt(oblique correlations between the factors) or to "first" or "second" in which case the general factor is equated with either the first or second group factor. A  message is issued suggesting that the model is not really well defined. This solution discussed in Zinbarg et al., 2007.  To do this in omega, add the option="first" or option="second" to the call.

Although obviously not meaningful for a 1 factor solution, it is of course possible to find the sum of the loadings on the first (and only) factor, square them, and compare them to the overall matrix variance.  This is done, with appropriate complaints.


In addition to \eqn{\omega_h}, another of McDonald's coefficients is \eqn{\omega_t}.  This is an estimate of the total reliability of a test. 

McDonald's \eqn{\omega_t}, which is similar to Guttman's \eqn{\lambda_6}, \code{\link{guttman}} but uses the estimates of uniqueness (\eqn{u^2} from factor analysis to find \eqn{e_j^2}. This is based on a decomposition of the variance of a test score, \eqn{V_x}  into four parts: that due to a general factor, \eqn{\vec{g}}, that due to a set of group factors, \eqn{\vec{f}},  (factors common to some but not all of the items), specific factors, \eqn{\vec{s}} unique to each item, and \eqn{\vec{e}}, random error.  (Because specific variance can not be distinguished from random error unless the test is given at least twice,  some combine these both into error). 

Letting \eqn{\vec{x} =  \vec{cg} + \vec{Af} + \vec {Ds} + \vec{e}}
then the communality of item\eqn{_j}, based upon general as well as group factors,
\eqn{h_j^2 = c_j^2 + \sum{f_{ij}^2}}
and the unique variance for the item
\eqn{u_j^2 = \sigma_j^2 (1-h_j^2)}
may be used to estimate the test reliability.
That is, if \eqn{h_j^2} is the communality of item\eqn{_j}, based upon general as well as group factors,  then for standardized items,  \eqn{e_j^2 = 1 - h_j^2} and
\deqn{
\omega_t = \frac{\vec{1}\vec{cc'}\vec{1} + \vec{1}\vec{AA'}\vec{1}'}{V_x} = 1 - \frac{\sum(1-h_j^2)}{V_x} = 1 - \frac{\sum u^2}{V_x}
}
Because \eqn{h_j^2 \geq r_{smc}^2}, \eqn{\omega_t \geq \lambda_6}.

It is important to distinguish here between the two \eqn{\omega} coefficients of McDonald, 1978 and Equation 6.20a of McDonald, 1999, \eqn{\omega_t} and \eqn{\omega_h}.  While the former is based upon the sum of squared loadings on all the factors, the latter is based upon the sum of the squared loadings on the  general factor. 
\deqn{\omega_h = \frac{ \vec{1}\vec{cc'}\vec{1}}{V_x}}

Another estimate reported is the omega for an infinite lenght test with a structure similar to the observed test.  This is found by 
\deqn{\omega_limit = \frac{ \vec{1}\vec{cc'}\vec{1}}{\vec{1}\vec{cc'}\vec{1} + \vec{1}\vec{AA'}\vec{1}'}}

}
\value{
 
  \item{omega hierarchical}{The \eqn{\omega_h} coefficient}
  \item{omega.lim}{The limit of \eqn{\omega_h} as the test becomes infinitly large}
  \item{omega total}{The \eqn{omega_t} coefficient}
  \item{alpha}{Cronbach's \eqn{\alpha}}
  \item{schmid}{The Schmid Leiman transformed factor matrix and associated matrices}
  \item{schmid$sl}{The g factor loadings as well as the residualized factors}
  \item{schmid$orthog}{Varimax rotated solution of the original factors}
  \item{schmid$oblique}{The oblimin or promax transformed factors}
  \item{schmid$phi}{the correlation matrix of the oblique factors}
  \item{schmid$gloading}{The loadings on the higher order, g, factor of the oblimin factors}
  \item{key}{A vector of -1 or 1 showing which direction the items were scored.}
  \item{model}{a matrix suitable to be given to the sem function for structure equation models}
  \item{various fit statistics}{various fit statistics, see output}
}


\references{ \url{http://personality-project.org/r/r.omega.html} \cr

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma.  Psychometrika (in press).

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{http://personality-project.org/revelle/publications/iclust.pdf})

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{http://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

Zinbarg, R., Yovel, I. & Revelle, W.  (2007).  Estimating omega  for structures containing two group factors:  Perils and prospects.  Applied Psychological Measurement. 31 (2), 135-157.
 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814 \url{http://apm.sagepub.com/cgi/reprint/30/2/121}
}
\author{ 
 \url{http://personality-project.org/revelle.html} \cr
Maintainer: William Revelle  \email{  revelle@northwestern.edu  } 
}
\note{Requires the GPArotation package.

The default rotation uses oblimin from the GPArotation package.  Alternatives include the simplimax function, as well as \code{\link{Promax}}.

If the factor solution leads to an exactly orthogonal solution (probably only for demonstration data sets), then use the rotate="Promax" option to get a solution.}

\seealso{ \code{\link{omega.graph}} \code{\link{ICLUST}}, \code{\link{ICLUST.graph}}, \code{\link{VSS}}, \code{\link{schmid} }, \code{\link{make.hierarchical} }}
\examples{
\dontrun{
test.data <- Harman74.cor$cov
if(!require(GPArotation) {message("Omega requires GPA rotation" )} else {my.omega <- omega(test.data,3)       
print(my.omega,digits=2)
#
         }
#create 9 variables with a hierarchical structure
jen.data <- sim.hierarchical()
#with correlations of
jen.data
#find omega 
if(require(Rgraphviz)) {jen.omega <- omega(jen.data,digits=2)} else {jen.omega <- omega(jen.data,digits=2,plot=FALSE)}
jen.omega

#create 8 items with a two factor solution, showing the use of the flip option
sim2 <- item.sim(8)
omega(sim2)   #an example of misidentification-- remember to look at the loadings matrices.
omega(sim2,2)  #this shows that in fact there is no general factor
omega(sim2,2,option="first") #but, if we define one of the two group factors as a general factor, we get a falsely high omega 
#apply omega to analyze 6 mental ability tests 
data(ability.cov)   #has a covariance matrix
if(require(Rgraphviz)) {omega(ability.cov$cov)} else {omega(ability.cov$cov,plot=FALSE)} 
}
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
