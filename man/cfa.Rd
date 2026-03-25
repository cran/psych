\name{CFA}
\alias{CFA}
\alias{CFA.bifactor}
\alias{irt.CFA}

\title{Confirmatory Factor Analysis Without Iteration}
\description{Factor analysis is a procedure for identifying latent variables  thought to account for the correlations or covariances between observed variables. There are two approaches to factor analysis: Exploratory Factor Analysis (e.g., EFA using the \code{\link{fa}} function)  and Confirmatory Factor Analysis (CFA).  Perhaps the best way to do Confirmatory Factor Analysis is with the laavan package's cfa function.  \code{\link{CFA}} in psych is a simple and more limited version for those who want to stay within the psych package and take advantage of various psych package options.  \code{\link{CFA}} uses the direct approach (the multiple group method) using the Spearman/Guttman approach as discussed by Dhaene and Rosseel, 2025. 
}
\usage{
CFA(model=NULL,x=NULL, all=FALSE, cor = "cor", use ="pairwise", n.obs = NA, 
 orthog=FALSE,  weight=NULL,correct=0, method="regression", 
 missing=FALSE,impute="none",Grice=FALSE,n.iter=100)

CFA.bifactor(model=NULL,x=NULL,all=FALSE,g=FALSE, cor="cor", use="pairwise", n.obs=NA, 
  orthog=FALSE, weight=NULL, correct=0, method="regression",
   missing=FALSE,impute="none",
  Grice=FALSE ,
  n.iter=1)
  
irt.CFA(model=NULL,x=NULL,all=TRUE,n.obs=NA, orthog=FALSE, weight=NULL, correct=0, 
    method="regression", 
	missing=FALSE,impute="none",
	Grice=FALSE,
	plot=TRUE,
	n.iter=1,
	scores=FALSE)
}
\arguments{
\item{model}{If specified, the model can either be in lavaan syntax or a keys list such as those used in \code{\link{scoreItems}} or \code{\link{scoreOverlap}}.  See examples. If the model is missing, then the analysis will be a for one factor model of all of the data. }
  \item{x}{A data matrix or correlation /covariance matrix. (If missing, then assumed to be the values of the first object.) }
  
  \item{all}{if TRUE, then do the analysis for  all the variables in the r matrix.  If FALSE, then select just those variables in the  r matrix defined by the model.}
  \item{g}{if TRUE, find a hierarchical (higher level) model.}
  \item{cor}{How to find the correlations: "cor" is Pearson", "cov" is covariance, 
  "tet" uses \code{\link{tetrachoric}}, "poly" uses \code{\link{polychoric}}, "mixed" uses \code{\link{mixedCor}} for a mixture of tetrachorics, polychorics, Pearsons, biserials, and polyserials, Yuleb is Yulebonett, Yuleq and YuleY are the obvious Yule coefficients as appropriate}
  \item{n.obs}{Number of observations if given a correlation/covariance matrix, defaults to 100}
  \item{orthog}{Should the factors be allowed to correlate (orthog=FALSE) or forced to be orthogonal (orthog=TRUE)}
   \item{use}{How to treat missing data, use="pairwise" is the default".  See cor for other options.}
    \item{weight}{If not NULL, a vector of length n.obs that contains weights for each observation. The NULL case is equivalent to all cases being weighted 1.} 
   \item{correct}{correction value for 0 values in tetrachoric correlation.  See the discussion in \code{\link{tetrachoric}} for alternative values.}
   \item{method}{Correlations are by default found using Pearson.  Alternative methods for the correlation  may be Spearman or Kendall.}
   \item{missing}{if r is a data matrix, and missing=TRUE, then impute missing values using either the median or the mean. Specifying impute="none" will not impute data points and thus will have some missing data in the factor scores.}
  \item{impute}{"median" or "mean" values are used to replace missing values}
  \item{Grice}{If TRUE, use the Grice method for factor indeterminacy. }
  \item{plot}{should irt.CFA display the plots (defaults to TRUE)}
  \item{n.iter}{Number of bootstrap resamples of the CFA solution. Defaults to 100.  }
  \item{scores}{If true, and given raw data, find irt based scores for the CFA based factors.}
  }
  \details{Most EFA and CFA functions use maximum likelihood functions to estimate the coefficients.  However, as Maccallum et al. (2007), and Dhaene and Rosseel (2024) point out, ML approaches are not necessarily optimal for finite (e.g., small) samples.  Maccallum et al. (2007) discuss why ML fails on some problems that minres procedures do not. 

  Confirmatory factor analysis may be done without iteration (and thus not using Maximum Likelihood procedures) by using some very old techniques.  The algorithm follows that of Dhaene and Rosseel (2024) using the ``Spearman" Multiple Group Method to estimate the communalities.  This method was introduced by Guttman (1952) and is discussed by Harman (1967).

\code{\link{CFA}} follows the Spearman approach for communalities discussed by Dhaene and Rosseel (2024) and described as the ``Multiple Group Method".  I use the upper case name (CFA) to avoid conflicts with lavaan's cfa function. Following Harman (1967)	 (Chapter 7, p 115-117) the communality of each variable is estimated by the ratio of the sum of all the correlations to the sum of squared correlations with that  variable.  The square root of the communality is the factor loading.  

Guttman (1952) points out that a weighting matrix of -1, 0, and 1 times is essentially a regression model where the use of differential weights doesn't make much difference. 

\code{\link{CFA.bifactor}} first does a CFA on all of the variables, and then does another CFA using the model matrix or keys list on the residual correlation matrix.  The results are in relatively close agreement with those from lavaan, but are not identical.  

To do a "S-1" solution (Eid, 2017; Li and Savalei, 2025) just specify a model with not all variables defined as group factors. 

Options for \code{\link{CFA.bifactor}} include solving the correlations as simple bi-factor model, or as a hierarchical/higher level model using the g=TRUE option.  Graphical output in the examples shows the difference of the two approaches.

Guttman and Harman's original method seem to be restricted to positive manifolds and finds the communalities based upon the the correlations. 

\eqn{h_i^2 =\frac{(\Sigma r_i)^2 -\Sigma r_i^2}{2(\Sigma r_{i<j}-\Sigma{r_i})}
{h_i^2 =(\Sigma r_i)^2 -\Sigma r_i^2}/{2(\Sigma r_{i<j}-\Sigma{r_i})}
}

 CFA estimates communalities using the absolute values of r when finding the sums. This allows applying the method to personality data sets such as the \code{\link{bfi}} or sapa  data sets as well as mood data as found in the msqR data sets (sapa and msqR are in psychTools).
 
 If the g parameter is set to true, a hierarchical or second order solution is found by first factoring all the variables for a one factor (g) model  and then factoring  the residualized matrix using the model based factors. This is shown in the test.hi example.
 
\code{\link{irt.CFA}} does a confirmatory analysis (CFA) on the polychoric/tetrachoric correlations  and then reports the results in terms of an \code{\link{irt.fa}} (IRT) model. See  \code{\link{irt.fa}} for a discussion of the equivalence of a 2PL IRT model and a (non-linear) factor analysis of items using  polychoric or tetrachoric correlations.

By default, the \code{\link{CFA}} will do n.iter=100 bootstrap resamplings to provide confidence intervals for the solutions.  Because \code{\link{CFA.bifactor}} uses \code{\link{CFA}} to find the general factor and then the residual group factors, confidence intervals are also reported.

The loadings and goodness of fits are based upon the median bootstrapped values.
} 

\value{
\item{loadings }{Factor (Structure) Loadings}
\item{Pattern}{Factor Pattern coefficients}
\item{Phi}{Factor correlations}
\item{communalities}{As estimated using the Spearman/Guttman procedure}
\item{dof}{Degrees of freedom is the number of original corrlations - number of loadings - number of between factor correlations}
\item{stats}{as found by \code{\link{fa.stats}} }
\item{scores}{Factor scores. }
\item{...}{Many other statistics as reported by \code{\link{fa}} }
\item{Call}{echoes the call to the function}
}
\references{Sara Dhaene and Yves Rosseel, 2024, An Evaluation of non-iterative estimators in confirmatory factor analysis. Structural Equation Modeling (31) 1 1-13. 
doi:  10.1080/10705511.2023.2187285

Michael Eid, Christian Geiser, Tobias Koch and Moritz Heene (2017) Anomalous Results in G-Factor Models: Explanations and Alternatives. Psychological Methods, 22, 541-562 

Guttman. L. (1952) Multiple group methods for common-factor analysis:  their basis, computation, and interpretation. Psychometrika, 17, (2) 209-222.

H.H. Harman (1967) Modern Factor Analysis. University of Chicago Press.

Li, Sijia and Savalei, Victoria  (2025), Evaluating Statistical Fit of Confirmatory Bifactor Models: Updated Recommendations and a Review of Current Practice.  Psychological Methods. doi.org/10.1037/met0000730 


MacCallum, Robert C. and Browne, Michael W. and Cai, Li  (2007) Factor analysis models as approximations. In Cudeck, Robert and MacCallum, Robert C. (Eds). Factor analysis at 100: Historical developments and future directions. Lawrence Erlbaum Associates Publishers.
}


\author{William Revelle
}
\note{The examples include a number of comparisons with the output of the lavaan package.  These are not run, but can be examined after loading lavaan.   The general observation is that the results are very similar, but not identical. The loadings are identical for the 9 variable Thurstone problem, but differ slightly for the 24 Holzinger problem. lavaan has several ways of estimating coefficients.  The ULS results match CFA most closely. 

Further note that cross loadings are not allowed.




}

\seealso{ \code{\link{fa}} for exploratory analysis and more discussion of factor analysis in general. \code{\link{omegaStats}} to allow quick comparisons with other functions.
\code{\link{irt.CFA}} calculates and displays the CFA analysis of tetra/polychoric correlations in the IRT metric.

}
\examples{
 #test set from Harman Table 7.1 P 116
har5 <- structure(c(1, 0.485, 0.4, 0.397, 0.295, 0.485, 1, 0.397, 0.397, 
0.247, 0.4, 0.397, 1, 0.335, 0.275, 0.397, 0.397, 0.335, 1, 0.195, 
0.295, 0.247, 0.275, 0.195, 1), dim = c(5L, 5L), dimnames = list(
    c("V1", "V2", "V3", "V4", "V5"), c("V1", "V2", "V3", "V4", 
    "V5")))
#in all of these examples, I am setting n.iter to be very small in order to pass CRAN test
CFA(har5,n.iter=5)   #The Harman example.   Note that the model not necessary for the 1 factor case.

CFA(Harman_5,n.iter=5)  #the Harman example of a Heywood case

v9 <- sim.hierarchical()  #Create a 3 correlated factor model using default values
model <- 'F1=~ V1 + V2 + V3
          F2=~ V4 + V5 + V6
          F3 =~ V7 +V8 + V9'
CFA(model,v9,n.iter=5)


model9 <- 'F1 =~ .9*V1 + .8*V2 + .7*V3
           F2 =~ .8*V4 + .7*V5 +.6*V6
           F3 =~ .7*V7 + .6*V8 +.5*V9
           F1 ~ .6*F2 + .5*F3
           F2 ~  .4*F3'
#An alternative way to create 3 correlated factors
#note that CFA drops the coefficients, the model is for generating the data
 #lavaan does not drop coefficients
v9s <- sim(model9,n=500)
 test <- CFA(model,v9s$observed , n.iter=5)  #do a cfa using Lavaan syntax
 test.bi <- CFA.bifactor(model9,v9,n.iter=5)
 test.hi <- CFA.bifactor(model9,v9,g=TRUE,n.iter=5)

#graphic displays make the output more understandable.
diagram(test)   #show three correlated factors
diagram(test.bi) #show the bifactor solution
diagram(test.hi) #show the hierarchical/higher order solution

#this next example requires psychTools  not run
#for a four factor model using keys

#CFA(psychTools::ability.keys[-1],psychTools::ability,  cor="tet")

CFA(bfi.keys,bfi, n.iter=5)   # a five factor model of the bfi items
#irt.CFA(bfi.keys,bfi[1:25])

colnames(Thurstone) <- rownames(Thurstone) <- paste0("x",1:9  ) #to match lavaan syntax
model <- HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
c3 <- CFA(model,Thurstone,n.obs=213,n.iter=5)  #compare with the lavaan solution 
           #which has a smaller chi^2 

c3  #show the result	
diagram(c3)  #graphically display the result

c3.hi <- CFA.bifactor(model,Thurstone,n.obs=213,n.iter=5)

c5 <- irt.CFA(bfi.keys,bfi[1:25],n.iter=5)

#do not run the next examples, they require lavaan
#They compare lavaan cfa solutions to CFA

if(FALSE) {
#
#The next examples require lavaan and are thus not run
library(lavaan)               
#The basic lavaan example 
fit <- cfa(model,sample.cov=Thurstone,sample.nobs=213,std.lv=TRUE, estimator="ML")
factor.congruence(fit,c3)  #identical loadings to 2 decimals
round(fit@Model@GLIST$lambda-c3$loadings,4)   
#add the g factor
HS.model <- ' general =~  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 
              visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
g.fit <- cfa(HS.model,sample.cov=Thurstone,sample.nobs=213,std.lv=TRUE,orthogonal=TRUE)
fa.congruence(g.fit,c3.hi) #identical congruence to 2 decimals
round(g.fit@Model@GLIST$lambda-c3.hi$loadings,2)  #loadings with ULS are identicla

#All 24 variables from Harman 
harman24 <- psychTools::holzinger.raw[157:301,8:31]
colnames(harman24) <- paste0("v",1:24)
 mod.24<-'
  g=~v1+v2+v3+v4+v5+v6+v7+v8+v9+v10+v11+v12+v13+v14+v15+v16+v17+v18+v19+v20+v21+v22+v23+v24
  spatial =~ v1 + v2 + v3 + v4 
                verbal=~ v5 + v6 + v7 + v8 + v9
                perceptual =~ v10 + v11 + v12 + v13
                recognition =~ v14+v15 + v16  + v17
                memory =~ v18 + v19 + v20
               '
#the lavaan solution (requires lavaan)
lav.har.uls <- cfa(mod.24, data=harman24,std.lv=TRUE,std.ov=TRUE, orthogonal=TRUE,
           estimator="ULS")

lav.har.ml <-cfa(mod.24, data=harman24,std.lv=TRUE,std.ov=TRUE,orthogonal=TRUE)                

 model.har24.5  <- 'spatial =~ v1 + v2 + v3 + v4 
                verbal=~ v5 + v6 + v7 + v8 + v9
                perceptual =~ v10 + v11 + v12 + v13
                recognition =~ v14+v15 + v16  + v17
                memory =~ v18 + v19 + v20'
 cfa.har24 <- CFA(model.har24.5,harman24) 
 cfa.har.bi <- CFA.bifactor(model.har24.5,harman24) 

 factor.congruence(list(lav.har.uls,lav.har.ml,cfa.har.bi)) #g is very good f1-4 very good

round(lav.har@Model@GLIST$lambda-cfa.har.bi$loadings,2)   #not the same
    }           
}

\keyword{ multivariate }
\keyword{ models}


