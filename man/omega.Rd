\name{omega}
\alias{omega}
\alias{omegaSem}
\alias{omegaFromSem}
\alias{omegah}
\alias{omegaDirect}
\alias{directSl}
\title{  Calculate McDonald's  omega estimates of general and total factor saturation }
\description{McDonald has proposed coefficient omega as an estimate of the general factor saturation of a test.  One way to find omega is to do a factor analysis of the original data set, rotate the factors obliquely, do a Schmid Leiman transformation, and then find omega. This function estimates omega as suggested by McDonald by using hierarchical factor analysis (following Jensen). A related option is to define the model using omega and then perform a confirmatory (bi-factor) analysis using the sem or lavaan packages.  This is done by omegaSem and omegaFromSem. omegaFromSem will convert appropriate sem/lavaan objects to find omega.  Yet another option is to do the direct Schmid-Leiman of Waller.
}
\usage{
omega(m,nfactors=3,fm="minres",n.iter=1,p=.05,poly=FALSE,key=NULL,
    flip=TRUE,digits=2, title="Omega",sl=TRUE,labels=NULL,
    plot=TRUE,n.obs=NA,rotate="oblimin",Phi=NULL,option="equal",covar=FALSE, ...)
omegaSem(m,nfactors=3,fm="minres",key=NULL,flip=TRUE,digits=2,title="Omega",
  sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin",
  Phi = NULL, option="equal",lavaan=TRUE,...)
  
omegah(m,nfactors=3,fm="minres",key=NULL,flip=TRUE, 
digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,
   n.obs=NA,rotate="oblimin",Phi = NULL,option="equal",covar=FALSE,...) 

omegaFromSem(fit,m=NULL,flip=TRUE,plot=TRUE)
omegaDirect(m,nfactors=3,fm="minres",rotate="oblimin",cut=.3,
   plot=TRUE,main="Direct Schmid Leiman")
directSl(m,nfactors=3,fm="minres",rotate="oblimin",cut=.3)
}


\arguments{
  \item{m}{A correlation matrix, or a data.frame/matrix of data, or (if Phi) is specified, an oblique factor pattern matrix }
  \item{nfactors}{Number of factors believed to be group factors}
 \item{n.iter}{How many replications to do in omega for bootstrapped estimates}
  \item{fm}{factor method (the default is minres)  fm="pa" for principal axes, fm="minres" for a minimum residual (OLS) solution, fm="pc" for principal components (see note), or  fm="ml" for maximum likelihood.}
  \item{poly}{should the correlation matrix be found using polychoric/tetrachoric or normal Pearson correlations}
  \item{key}{a vector of +/- 1s to specify the direction of scoring of items.  The default is to assume all items are positively keyed, but if some items are reversed scored, then key should be specified.}
  \item{flip}{If flip is TRUE, then items are automatically flipped to have positive correlations on the general factor. Items that have been reversed are shown with a - sign.}
 
  \item{p}{probability of two tailed conference boundaries}
  \item{digits}{if specified, round the output to digits}
  \item{title}{Title for this analysis}
  \item{main}{main for this analysis  (directSl)}
  \item{cut}{Loadings greater than cut are used in directSl}
  \item{sl}{If plotting the results, should the Schmid Leiman solution be shown or should the hierarchical solution be shown? (default sl=TRUE)}
  \item{labels}{If plotting, what labels should be applied to the variables? If not specified, will default to the column names.}
  \item{plot}{plot=TRUE (default) calls omega.diagram, plot =FALSE does not.  If Rgraphviz is available, then \code{\link{omega.graph}} may be used separately.}
  \item{n.obs}{Number of observations - used for goodness of fit statistic}
  \item{rotate}{What rotation to apply? The default is oblimin, the alternatives include simplimax, Promax,  cluster and target. target will rotate to an optional keys matrix (See \code{\link{target.rot}})}
  \item{Phi}{If specified, then omega is found from the pattern matrix (m) and the factor intercorrelation matrix (Phi).}
  \item{option}{In the two factor case (not recommended), should the loadings be equal, emphasize the first factor, or emphasize the second factor. See in particular the option parameter in  \code{\link{schmid}} for treating the case of two group factors.}
  \item{covar}{defaults to FALSE and the correlation matrix is found (standardized variables.)  If TRUE, the do the calculations on the unstandardized variables and use covariances.}
  \item{lavaan}{if FALSE, will use John Fox's sem package to do the omegaSem.  If TRUE, will use Yves Rosseel's lavaan package. }
  \item{fit}{The fitted object from lavaan or sem. For lavaan, this includes the correlation matrix and the variable names and thus m needs not be specified.}
  \item{...}{Allows additional parameters to be passed through to the factor routines.  }
  
}

\details{``Many scales are assumed by their developers and users to be primarily a measure of one latent variable. When it is also assumed that the scale conforms to the effect indicator model of measurement (as is almost always the case in psychological assessment), it is important to support such an 
interpretation with evidence regarding the internal structure of that scale. In particular, it is important to examine two related properties pertaining to the internal structure of such a scale. The first property relates to whether all the indicators forming the scale measure a latent variable in common. 

The second internal structural property pertains to the proportion of variance in the scale scores (derived from summing or averaging the indicators) accounted for by this latent variable that is common to all the indicators (Cronbach, 1951; McDonald, 1999; Revelle, 1979). That is, if an effect indicator scale is primarily a measure of one latent variable common to all the indicators forming the scale, then that latent variable should account for the majority of the variance in the scale scores. Put differently, this variance ratio provides important information about the sampling fluctuations when estimating individuals' standing on a latent variable common to all the indicators arising from the sampling of indicators (i.e., when dealing with either Type 2 or Type 12 sampling, to use the terminology of Lord, 1956). That is, this variance proportion can be interpreted as the square of the correlation between the scale score and the latent variable common to all the indicators in the infinite universe of indicators of which the scale indicators are a subset. Put yet another way, this variance ratio is important both as reliability and a validity coefficient. This is a reliability issue as the larger this variance ratio is, the more accurately one can predict an individual's relative standing on the latent variable common to all the scale's indicators based on his or her 
observed scale score. At the same time, this variance ratio also bears on the construct validity of the scale given that construct validity encompasses the internal structure of a scale." (Zinbarg, Yovel, Revelle, and McDonald, 2006).

McDonald has proposed coefficient omega_hierarchical (\eqn{\omega_h}) as an estimate of the general factor saturation of a test.  Zinbarg, Revelle, Yovel and Li (2005) 
\url{https://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf} compare McDonald's \eqn{\omega_h} to Cronbach's \eqn{\alpha} and Revelle's \eqn{\beta}.  They conclude that \eqn{\omega_h} is the best estimate. (See also Zinbarg et al., 2006 and Revelle and Zinbarg (2009)).   

One way to find \eqn{\omega_h}{omega_h} is to do a factor analysis of the original data set, rotate the factors obliquely, factor that correlation matrix, do a   Schmid-Leiman (\link{schmid}) transformation to find general factor loadings, and then find \eqn{\omega_h}{omega_h}.  Here we present code to do that.  

\eqn{\omega_h}{omega_h} differs as a function of how the factors are estimated.  Four options are available, three use the \code{\link{fa}} function but with different factoring methods: the default does a minres factor solution, fm="pa"  does a principle axes factor analysis  fm="mle" does a maximum likelihood solution; fm="pc" does a principal components analysis using (\code{\link{principal}}).  

For ability items, it is typically the case that all items will have positive loadings on the general factor.  However, for non-cognitive items it is frequently the case that some items are to be scored positively, and some negatively.  Although probably better to specify which directions the items are to be scored by specifying a key vector, if flip =TRUE (the default), items will be reversed so that they have positive loadings on the general factor.  The keys are reported so that scores can be found using the \code{\link{scoreItems}} function.  Arbitrarily reversing items this way can overestimate the general factor. (See the example with a simulated circumplex).



\eqn{\beta}{beta}, an alternative to \eqn{\omega_h}, is defined as the worst split half 
reliability (Revelle, 1979).  It can be estimated by using \code{\link{ICLUST}} (a 
hierarchical clustering algorithm originally developed for main frames and written in 
Fortran and that is now part of the psych package.  (For a very complimentary review of 
why the ICLUST algorithm is useful in scale construction, see Cooksey and Soutar, 2005)). 




The \code{\link{omega}} function uses exploratory factor analysis to estimate the \eqn{\omega_h} coefficient.  It is important to remember that  ``A recommendation that should be heeded, regardless of the method chosen to estimate \eqn{\omega_h}, is to always examine the pattern of the estimated general factor loadings prior to estimating \eqn{\omega_h}. Such an examination constitutes an informal test of the assumption that there is a latent variable common to all of the scale's indicators that can be conducted even in the context of EFA. If the loadings were salient for only a relatively small subset of the indicators, this would suggest that there is no true general factor underlying the covariance matrix. Just such an informal assumption test would have afforded a great deal of protection against the possibility of misinterpreting the misleading \eqn{\omega_h} estimates occasionally produced in the simulations reported here." (Zinbarg et al., 2006, p 137).

A simple demonstration of the problem of an omega estimate reflecting just one of two group factors can be found in the last example.  

Diagnostic statistics that reflect the quality of the omega solution include a comparison of the relative size of the g factor eigen value to the other eigen values, the percent of the common variance for each item that is general factor variance (p2), the mean of p2, and the standard deviation of p2.  Further diagnostics can be done by describing (\link{describe}) the $schmid$sl results.


Although omega_h is uniquely defined only for cases where 3 or more subfactors are extracted, it is sometimes desired to have a two factor solution.  By default this is done by forcing the schmid extraction to treat the two subfactors as having equal loadings.  

There are three possible options for this condition: setting the general factor loadings between the two lower order factors to be "equal" which will be the sqrt(oblique correlations between the factors) or to "first" or "second" in which case the general factor is equated with either the first or second group factor. A  message is issued suggesting that the model is not really well defined. This solution discussed in Zinbarg et al., 2007.  To do this in omega, add the option="first" or option="second" to the call.

Although obviously not meaningful for a 1 factor solution, it is of course possible to find the sum of the loadings on the first (and only) factor, square them, and compare them to the overall matrix variance.  This is done, with appropriate complaints.

In addition to \eqn{\omega_h}, another of McDonald's coefficients is \eqn{\omega_t}.  This is an estimate of the total reliability of a test. 



McDonald's \eqn{\omega_t}, which is similar to Guttman's \eqn{\lambda_6}, \code{\link{guttman}} but uses the estimates of uniqueness (\eqn{u^2}) from factor analysis to find \eqn{e_j^2}. This is based on a decomposition of the variance of a test score, \eqn{V_x}  into four parts: that due to a general factor, \eqn{\vec{g}}, that due to a set of group factors, \eqn{\vec{f}},  (factors common to some but not all of the items), specific factors, \eqn{\vec{s}} unique to each item, and \eqn{\vec{e}}, random error.  (Because specific variance can not be distinguished from random error unless the test is given at least twice,  some combine these both into error). 

Letting \eqn{\vec{x} =  \vec{cg} + \vec{Af} + \vec {Ds} + \vec{e}}{x = cg + Af + Ds + e}
then the communality of item\eqn{_j}, based upon general as well as group factors,
\eqn{h_j^2 = c_j^2 + \sum{f_{ij}^2}}{h_j^2 = c_j^2 + sum(f_ij^2)}
and the unique variance for the item
\eqn{u_j^2 = \sigma_j^2 (1-h_j^2)}
may be used to estimate the test reliability.
That is, if \eqn{h_j^2} is the communality of item\eqn{_j}, based upon general as well as group factors,  then for standardized items,  \eqn{e_j^2 = 1 - h_j^2} and
\deqn{
\omega_t = \frac{\vec{1}\vec{cc'}\vec{1} + \vec{1}\vec{AA'}\vec{1}'}{V_x} = 1 - \frac{\sum(1-h_j^2)}{V_x} = 1 - \frac{\sum u^2}{V_x}}{\omega_t = (1 cc' 1 + 1 AA' 1')/(V_x)}


Because \eqn{h_j^2 \geq r_{smc}^2}, \eqn{\omega_t \geq \lambda_6}.

It is important to distinguish here between the two \eqn{\omega} coefficients of McDonald, 1978 and Equation 6.20a of McDonald, 1999, \eqn{\omega_t} and \eqn{\omega_h}.  While the former is based upon the sum of squared loadings on all the factors, the latter is based upon the sum of the squared loadings on the general factor. 
\deqn{\omega_h = \frac{ \vec{1}\vec{cc'}\vec{1}}{V_x}}{\omega_h = (1 cc' 1')/Vx}

Another estimate reported is the omega for an infinite length test with a structure similar to the observed test (omega H asymptotic).  This is found by 
\deqn{\omega_{limit} = \frac{\vec{1}\vec{cc'}\vec{1}}{\vec{1}\vec{cc'}\vec{1} + \vec{1}\vec{AA'}\vec{1}'}}{\omega_{limit} = (1 cc' 1')/(1 cc' 1' + 1 AA' 1')}. 

Following suggestions by Steve Reise, the Explained Common Variance (ECV) is also reported.  This is the ratio of the general factor eigen value to the sum of all of the eigen values.  As such, it is a better indicator of unidimensionality than of the amount of test variance accounted for by a general factor.


The input to omega may be a correlation matrix or a raw data matrix, or a factor pattern matrix with the factor intercorrelations (Phi) matrix.  

\code{\link{omega}} is an exploratory factor analysis function that uses a Schmid-Leiman transformation.  \code{\link{omegaSem}} first calls \code{\link{omega}} and then takes the Schmid-Leiman solution, converts this to a confirmatory sem model and then calls the sem package to conduct a confirmatory model.  \eqn{\omega_h} is then calculated from the CFA output. Although for well behaved problems, the efa and cfa solutions will be practically identical, the CFA solution will not always agree with the EFA solution. In particular, the estimated   \eqn{R^2} will sometimes exceed 1. (An example of this is the Harman 24 cognitive abilities problem.)


In addition, not all EFA solutions will produce workable CFA solutions.  Model misspecifications will lead to very strange CFA estimates. 

It is also possible to give \code{\link{omega}} a factor pattern matrix and the associated factor intercorrelation.  In this case, the analysis will be done on these matrices.  This is particularly useful if one is not satisfied with the exploratory EFA solutions and rotation options and somehow comes up with an alternative. (For instance, one might want to do a EFA using fm='pa' with a Kaiser normalized Promax solution with a specified m value.)  

\code{\link{omegaFromSem}} takes the output from a sem model and uses it to find 
 \eqn{\omega_h}.  The estimate of factor indeterminacy, found by the multiple \eqn{R^2} of the variables with the factors, will not match that found by the EFA model.  In particular, the estimated   \eqn{R^2} will sometimes exceed 1. (An example of this is the Harman 24 cognitive abilities problem.)
 
The notion of omega may be applied to the individual factors as well as the overall test. A typical use of omega is to identify subscales of a total inventory.  Some of that variability is due to the general factor of the inventory, some to the specific variance of each subscale.  Thus, we can find a number of different omega estimates:   what percentage of the variance of the items identified with each subfactor is actually due to the general factor.  What variance is common but unique to the subfactor, and what is the total reliable variance of each subfactor.   These results are reported in omega.group object and in the last few lines of the normal output.

Finally, and still be tested, is \code{\link{omegaDirect}} adapted from Waller (2017).  This is a direct rotation to a Schmid-Leiman  like solution without doing the hierarchical factoring (\code{\link{directSl}}).  This rotation is then interpreted in terms of omega.  It is included here to allow for comparisons with the alternative procedures \code{\link{omega}} and \code{\link{omegaSem}}.  Preliminary analyses suggests that it produces inappropriate solutions for the case where there is no general factor.

Moral: Finding omega_h is tricky and one should probably compare \code{\link{omega}}, \code{\link{omegaSem}},  \code{\link{omegaDirect}} and even \code{\link{iclust}} solutions to understand the differences.

The summary of the omega object is a reduced set of the most useful output. 

The various objects returned from omega include:
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
  \item{model}{a list of two elements, one suitable to give  to the sem function for structure equation models, the other, to give to the lavaan package. }
  \item{sem}{The output from a sem analysis}
  \item{omega.group}{The summary statistics for the omega total, omega hierarchical (general) and omega within each group.}
  \item{scores}{Factor score estimates are found for the Schmid-Leiman solution.  To get 
  scores for the hierarchical model see the note.} 
  \item{various fit statistics}{various fit statistics, see output}
  \item{OmegaSem}{ is an object that contains the fits for the OmegaSem output.}
  \item{loadings}{The direct SL rotated object (from omegaDirect)}
  \item{orth.f}{The original, unrotated solution from omegaDirect}
  \item{Target}{The cluster based target for rotation in directSl}
}

\references{ \url{https://personality-project.org/r/r.omega.html} \cr

Jensen, Arthur R. and Li-Jen  Weng (1994) What is a good g?  Intelligence, 18, 3, 231-258

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{https://personality-project.org/revelle/publications/iclust.pdf})

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma.  Psychometrika, 74, 1, 145-154. (\url{https://personality-project.org/revelle/publications/rz09.pdf}

Waller, N. G. (2017) Direct {Schmid-Leiman} Transformations and Rank-Deficient Loadings Matrices.  Psychometrika.  DOI: 10.1007/s11336-017-9599-0

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{https://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

Zinbarg, R., Yovel, I. & Revelle, W.  (2007).  Estimating omega  for structures containing two group factors:  Perils and prospects.  Applied Psychological Measurement. 31 (2), 135-157.
 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814 \url{https://journals.sagepub.com/doi/10.1177/0146621605278814}
}
\author{ 
 \url{https://personality-project.org/revelle.html} \cr
Maintainer: William Revelle  \email{  revelle@northwestern.edu  } 
}
\note{Requires the GPArotation package.

The default rotation uses oblimin from the GPArotation package.  Alternatives include the simplimax function, as well as \code{\link{Promax}} or the \code{\link{promax}} rotations.  promax will do a Kaiser normalization before applying Promax rotation.

If the factor solution leads to an exactly orthogonal solution (probably only for demonstration data sets), then use the rotate="Promax" option to get a solution.

\code{\link{omegaSem}} requires the sem or lavaan packages.  \code{\link{omegaFromSem}} uses the output from the sem or lavaan package.

\code{\link{omega}} may be run on raw data (finding either  Pearson or  tetrachoric/polychoric corrlations, depending upon the poly option) a correlation matrix, a polychoric correlation matrix (found by e.g., \code{\link{polychoric}}), or the output of a previous omega run.  This last case is particularly useful when working with categorical data using the poly=TRUE option.  For in this case, most of the time is spent in finding the correlation matrix.  The matrix is saved as part of the omega output and may be used as input for subsequent runs.  A similar feature is found in \code{\link{irt.fa}} where the output of one analysis can be taken as the input to the subsequent analyses.  

However, simulations based upon tetrachoric and polychoric correlations suggest that although the structure is better defined, that the estimates of omega are inflated over the true general factor saturation.

Omega returns factor scores based upon the Schmid-Leiman transformation.  To get the hierarchical factor scores, it is necessary to do this outside of omega. See the example (not run).

Consider the case of the raw data in an object data.  Then 

f3 <- fa(data,3,scores="tenBerge", oblique.rotation=TRUE
f1 <- fa(f3$scores)
hier.scores <- data.frame(f1$scores,f3$scores)

When doing fm="pc", principal components are done for the original correlation matrix, but minres is used when examining the intercomponent correlations.   A warning is issued that the method was changed to minres for the higher order solution.  omega is a factor model, and finding loadings using principal components will overestimate the resulting solution.  This is particularly problematic for the amount of group saturation, and thus the omega.group statistics are overestimates.  

The last three lines of omega report "Total, General and Subset omega for each subset".  These are available as the omega.group object in the output.

The last of these (omega group)  is effectively what Steve Reise calls omegaS for the subset omega.

The omega general is the amount of variance in the group that is accounted for by the general factor, the omega total is the amount of variance in the group accounted for by general + group.

This is based upon a cluster solution (that is to say, every item is assigned to one group) and this is why for  first column the omega general and group do not add up to omega total.  Some of the variance is found in the cross loadings between groups.

Reise and others like to report the ratio of the second line to the first line (what portion of the reliable variance is general factor) and the third row to the first (what portion of the reliable variance is within group but not general.  This may be found by using the omega.group object that is returned by omega. (See the last example.)

If using the lavaan=TRUE option in \code{\link{omegaSem}} please note that variable names can not start with a digit (e.g. 4.Letter.Words in the  \code{\link{Thurstone}} data set. The leading digit needs to be removed. 

 \code{\link{omegaSem}} will do an exploratory efa and omega, create (and return) the commands for doing either a sem or lavaan analysis.  The commands are returned as the model object.  This can be used for further sem/lavaan analyses.  
 
 Omega can also be found from  an analysis  done using lavaan or sem directly by calling  \code{\link{omegaFromSem}} with the original correlation matrix and the fit of the sem/lavaan model.  See the last (not run) example)
}

\seealso{ \code{\link{omega.graph}} for a dot code graphic.   \code{\link{ICLUST}}, \code{\link{ICLUST.diagram}} for hierarchical cluster analysis, \code{\link{VSS}}, \code{\link{nfactors}}, \code{\link{fa.parallel}} for alternative ways of estimating the appropriate number of factors. \code{\link{schmid} } for the decomposition used in omega.
 \code{\link{make.hierarchical}} to simulate a hierarchical model.

 \code{\link{fa.multi}} for hierarchical factor analysis with an arbitrary number of 2nd  order factors. 
 }
\examples{
\dontrun{
 test.data <- Harman74.cor$cov
# if(!require(GPArotation)) {message("Omega requires GPA rotation" )} else {
      my.omega <- omega(test.data)       
      print(my.omega,digits=2)
#}
 
#create 9 variables with a hierarchical structure
v9 <- sim.hierarchical()  
#with correlations of
round(v9,2)
#find omega 
v9.omega <- omega(v9,digits=2)
v9.omega

#create 8 items with a two factor solution, showing the use of the flip option
sim2 <- item.sim(8)
omega(sim2)   #an example of misidentification-- remember to look at the loadings matrices.
omega(sim2,2)  #this shows that in fact there is no general factor
omega(sim2,2,option="first") #but, if we define one of the two group factors 
     #as a general factor, we get a falsely high omega 
#apply omega to analyze 6 mental ability tests 
data(ability.cov)   #has a covariance matrix
omega(ability.cov$cov)

#om <- omega(Thurstone)
#round(om$omega.group,2)
#round(om$omega.group[2]/om$omega.group[1],2)  #fraction of reliable that is general variance
# round(om$omega.group[3]/om$omega.group[1],2)  #fraction of reliable that is group variance

#To find factor score estimates for the hierarchical model it is necessary to 
#do two extra steps.

#Consider the case of the raw data in an object data.  (An example from simulation)
# set.seed(42)
# gload <- matrix(c(.9,.8,.7),nrow=3)
# fload <- matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,rep(0,9),.7,.6,.4),   ncol=3)
# data <- sim.hierarchical(gload=gload,fload=fload, n=100000, raw=TRUE)
# 
# f3 <- fa(data$observed,3,scores="tenBerge", oblique.scores=TRUE)
# f1 <- fa(f3$scores)

# om <- omega(data$observed,sl=FALSE) #draw the hierarchical figure
# The scores from om are based upon the Schmid-Leiman factors and although the g factor 
# is identical, the group factors are not.
# This is seen in the following correlation matrix
# hier.scores <- cbind(om$scores,f1$scores,f3$scores)
# lowerCor(hier.scores)
#
#this next set of examples require lavaan
#jensen <- sim.hierarchical()   #create a hierarchical structure (same as v9 above)
#om.jen <- omegaSem(jensen,lavaan=TRUE)  #do the exploratory omega with confirmatory as well
#lav.mod <- om.jen$omegaSem$model$lavaan #get the lavaan code or create it yourself
# lav.mod <- 'g =~ +V1+V2+V3+V4+V5+V6+V7+V8+V9
#              F1=~  + V1 + V2 + V3             
#              F2=~  + V4 + V5 + V6 
#              F3=~  + V7 + V8 + V9 '  
#lav.jen <- cfa(lav.mod,sample.cov=jensen,sample.nobs=500,orthogonal=TRUE,std.lv=TRUE)
# omegaFromSem(lav.jen,jensen)
#the directSl solution
#direct.jen <- directSl(jen)
#direct.jen 

#try a one factor solution -- this is not recommended, but sometimes done
#it will just give omega_total
# lav.mod.1 <- 'g =~ +V1+V2+V3+V4+V5+V6+V7+V8+V9 '  
#lav.jen.1<- cfa(lav.mod.1,sample.cov=jensen,sample.nobs=500,orthogonal=TRUE,std.lv=TRUE)
# omegaFromSem(lav.jen.1,jensen)



}
}
\keyword{ multivariate }
\keyword{ models }
