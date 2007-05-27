\name{omega}
\alias{omega}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{   Calculate the omega estimate of factor saturation}
\description{McDonald has proposed coefficient omega as an estimate of the general factor saturation of a test.  One way to find omega is to do a factor analysis of the original data set, rotate the factors obliquely, do a Schmid Leiman transformation, and then find omega. This function estimates omega as suggested by McDonald by using hierarchical factor analysis (following Jensen).
}
\usage{
omega(m, nfactors, pc = "pa",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A correlation matrix or a data.frame/matrix of data }
  \item{nfactors}{Number of factors believed to be group factors}
  \item{pc}{  pc="pa" for principal axes, pc="pc" for principal components, pc="mle" for maximum likelihood .}
  \item{...}{Allows additional parameters to be passed through to the factor routines}
}

\details{``Many scales are assumed by their developers and users to be primarily a measure of one latent variable. When it is also assumed that the scale conforms to the effect indicator model of measurement (as is almost always the case in psychological assessment), it is important to support such an 
interpretation with evidence regarding the internal structure of that scale. In particular, it is important to examine two related properties pertaining to the internal structure of such a scale. The first property relates to whether all the indicators forming the scale measure a latent variable in common. 

The second internal structural property pertains to the proportion of variance in the scale scores (derived from summing or averaging the indicators) accounted for by this latent variable that is common to all the indicators (Cronbach, 1951; McDonald, 1999; Revelle, 1979). That is, if an effect indicator scale is primarily a measure of one latent variable common to all the indicators forming the scale, then that latent variable should account for the majority of the variance in the scale scores. Put differently, this variance ratio provides important information about the sampling fluctuations when estimating individuals' standing on a latent variable common to all the indicators arising from the sampling of indicators (i.e., when dealing with either Type 2 or Type 12 sampling, to use the terminology of Lord, 1956). That is, this variance proportion can be interpreted as the square of the correlation between the scale score and the latent variable common to all the indicators in the infinite universe of indicators of which the scale indicators are a subset. Put yet another way, this variance ratio is important both as reliability and a validity coefficient. This is a reliability issue as the larger this variance ratio is, the more accurately one can predict an individual's relative standing on the latent variable common to all the scale's indicators based on his or her 
observed scale score. At the same time, this variance ratio also bears on the construct validity of the scale given that construct validity encompasses the internal structure of a scale." (Zinbarg, Yovel, Revelle, and McDonald, 2006).


McDonald has proposed coefficient omega as an estimate of the general factor saturation of a test.  Zinbarg, Revelle, Yovel and Li (2005) 
\url{http://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf} compare McDonald's Omega to Cronbach's alpha and Revelle's beta.  They conclude that omega is the best estimate. (See also Zinbarg et al., 2006)   

One way to find omega is to do a factor analysis of the original data set, rotate the factors obliquely, do a   Schmid-Leiman (\link{schmid}) transformation, and then find omega.  Here we present code to do that.  

Omega differs as a function of how the factors are estimated.  Three options are available, pc="pa"  does a principle axes factor analysis (\link{factor.pa}), pc="mle" uses the factanal function, and pc="pc" does a principal components analysis (\link{principal}).  

Beta, an alternative to omega, is defined as the worst split half reliability.  It can be estimated by using \link{ICLUST} (a hierarchical clustering algorithm originally developed for main frames and written in Fortran and that is now available in R.  (For a very complimentary review of why the ICLUST algorithm is useful in scale construction, see Cooksey and Soutar, 2005). 

}
\value{
 
  \item(omega) {The omega coefficient}
  \item{alpha }{Cronbach's alpha}
  \item{schmid }{The Schmid Leiman transformed factor matrix and associated matrices}
  \item{schmid$sl}{The g factor loadings as well as the residualized factors}
  \item{schmid$orthog}{Varimax rotated solution of the original factors}
  \item{schmid$oblique}{The oblimin transformed factors}
  \item{schmid$fcor}{the correlation matrix of the oblique factors}
  \item{schid$gloading}{The loadings on the higher order, g, factor of the oblimin factors}
}
\references{ \url{http://personality-project.org/r/r.omega.html} \cr

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{http://personality-project.org/revelle/publications/iclust.pdf})

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{http://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814 \url{http://apm.sagepub.com/cgi/reprint/30/2/121}
}
\author{ 
 \url{http://personality-project.org/revelle.html} \cr
Maintainer: William Revelle  \email{  revelle@northwestern.edu           } 
}
\note{Requires the GPArotation package }
\seealso{ \code{\link{ICLUST}}, \code{\link{ICLUST.graph}}, \code{\link{VSS}}, \code{\link{schmid} }}
\examples{
\dontrun{
test.data <- Harman74.cor$cov
my.omega <- omega(test.data,3)       
print(my.omega,digits=2)
}
#produces this output


#$omega
#[1] 0.64
#
#$alpha
#[1] 0.91
#
#$schmid
#$schmid$sl
#                       g factor Factor1 Factor2 Factor3    h2   u2
#VisualPerception           0.53   0.018  0.4688 0.02089 0.494 0.51
#Cubes                      0.34   0.022  0.3029 0.04544 0.209 0.79
#PaperFormBoard             0.38   0.033  0.3971 0.18505 0.398 0.60
#Flags                      0.43   0.109  0.3233 0.06148 0.261 0.74
#GeneralInformation         0.57   0.564  0.0078 0.09900 0.606 0.39
#PargraphComprehension      0.57   0.599  0.0244 0.02960 0.671 0.33
#SentenceCompletion         0.56   0.624  0.0352 0.02783 0.730 0.27
#WordClassification         0.56   0.394  0.1343 0.09255 0.341 0.66
#WordMeaning                0.58   0.637  0.0135 0.06034 0.762 0.24
#Addition                   0.35   0.047  0.0939 0.81706 0.858 0.14
#Code                       0.44   0.061  0.1442 0.44377 0.300 0.70
#CountingDots               0.37   0.100  0.1576 0.57732 0.491 0.51
#StraightCurvedCapitals     0.50   0.036  0.2677 0.33168 0.301 0.70
#WordRecognition            0.34   0.127  0.1506 0.10015 0.093 0.91
#NumberRecognition          0.32   0.060  0.2005 0.07972 0.105 0.90
#FigureRecognition          0.44   0.022  0.4080 0.00192 0.374 0.63
#ObjectNumber               0.37   0.064  0.1769 0.22163 0.139 0.86
#NumberFigure               0.43   0.076  0.3290 0.26170 0.339 0.66
#FigureWord                 0.37   0.053  0.2431 0.10447 0.151 0.85
#Deduction                  0.53   0.231  0.2814 0.00299 0.277 0.72
#NumericalPuzzles           0.50   0.025  0.2877 0.30211 0.301 0.70
#ProblemReasoning           0.52   0.222  0.2840 0.00067 0.272 0.73
#SeriesCompletion           0.59   0.198  0.3304 0.07553 0.325 0.68
#ArithmeticProblems         0.52   0.211  0.1106 0.40982 0.320 0.68
#
#$schmid$orthog
#                       Factor1 Factor2  Factor3
#VisualPerception         0.025   0.702 -0.02336
#Cubes                    0.030   0.454 -0.05080
#PaperFormBoard           0.045   0.595 -0.20689
#Flags                    0.149   0.484 -0.06874
#GeneralInformation       0.771  -0.012  0.11068
#PargraphComprehension    0.817   0.037 -0.03310
#SentenceCompletion       0.852  -0.053  0.03111
#WordClassification       0.538   0.201  0.10348
#WordMeaning              0.870   0.020 -0.06746
#Addition                 0.065  -0.141  0.91350
#Code                     0.083   0.216  0.49615
#CountingDots            -0.136   0.236  0.64546
#StraightCurvedCapitals   0.049   0.401  0.37082
#WordRecognition          0.173   0.226  0.11197
#NumberRecognition        0.082   0.300  0.08912
#FigureRecognition       -0.029   0.611  0.00214
#ObjectNumber             0.087   0.265  0.24779
#NumberFigure            -0.104   0.493  0.29259
#FigureWord               0.073   0.364  0.11681
#Deduction                0.315   0.421 -0.00334
#NumericalPuzzles         0.035   0.431  0.33777
#ProblemReasoning         0.303   0.425 -0.00074
#SeriesCompletion         0.270   0.495  0.08445
#ArithmeticProblems       0.288   0.166  0.45820
#
#$schmid$fcor
#     [,1] [,2] [,3]
#[1,] 1.00 0.51 0.30
#[2,] 0.51 1.00 0.33
#[3,] 0.30 0.33 1.00
#
#$schmid$gloading
#
#Loadings:
#     Factor1
#[1,] 0.681  
#[2,] 0.744  
#[3,] 0.447  
#
#               Factor1
#SS loadings      1.218
#Proportion Var   0.406
#
#
# The function is currently defined as
function(m,nfactors=3,pc="pa",...) {
      #m is a correlation matrix
      #nfactors is the number of factors to extract
      require(GPArotation)
      nvar <-dim(m)[2]
      gf<-schmid(m,nfactors,pc,...)
      Vt <- sum(m)   #find the total variance in the scale
      Vitem <-sum(diag(m)) #
      gload <- gf$sl[,1]
      gsq <- (sum(gload))^2
      alpha <- ((Vt-Vitem)/Vt)*(nvar/(nvar-1))
      omega <- list(omega= gsq/Vt,alpha=alpha,schmid=gf)
      }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
