\name{mediate}
\alias{mediate}
\alias{mediate.diagram}
\alias{moderate.diagram}
\title{Estimate and display direct and indirect effects of mediators and  moderator in path models}
\description{
Find the direct and indirect effects of a predictor in path models of mediation and moderation. Bootstrap confidence intervals for the indirect effects.  Mediation models are just extended regression models making explicit the effect of particular covariates in the model. Moderation is done by multiplication of the predictor variables.   This function supplies basic mediation/moderation analyses for some of the classic problem types. 
}
\usage{
mediate(y, x, m, data, mod = NULL, n.obs = NULL, use = "pairwise", n.iter = 5000, 
       alpha = 0.05, std = FALSE,plot=TRUE)
mediate.diagram(medi,digits=2,ylim=c(3,7),xlim=c(-1,10),show.c=TRUE,
     main="Mediation model",...)
moderate.diagram(medi,digits=2,ylim=c(2,8),main="Moderation model",...)
}

\arguments{
  \item{y}{The dependent variable (or a formula suitable for a linear model)}
  \item{x}{One or more predictor variables}
  \item{m}{One (or more) mediating variables}
  \item{data}{A data frame holding the data or a correlation or covariance  matrix. }
  \item{mod}{A moderating variable, if desired}
  \item{n.obs}{If the data are from a correlation or covariance matrix, 
     how many observations were used. This will lead to simulated data for the bootstrap.}
  \item{use}{use="pairwise" is the default when finding correlations or covariances}
  \item{n.iter}{Number of bootstrap resamplings to conduct}
  \item{alpha}{Set the width of the confidence interval to be 1 - alpha}
  \item{std}{standardize the covariances to find the standardized betas}
  \item{plot}{Plot the resulting paths}
  \item{digits}{The number of digits to report in the mediate.diagram.} 
  \item{medi}{The output from mediate may be imported into mediate.diagram}
  \item{ylim}{The limits for the y axis in the mediate and moderate diagram functions}
  \item{xlim}{The limits for the x axis.  Make the minimum more negative if the x by x correlations do not fit.}
  \item{show.c}{If FALSE, do not draw the c lines, just the partialed (c')  lines}
  \item{main}{The title for the mediate and moderate functions}
  \item{...}{Additional graphical parameters to pass to mediate.diagram}
}
\details{
When doing linear modeling, it is frequently convenient to estimate the direct effect of a predictor controlling for the indirect effect of a mediator.  See Preacher and Hayes (2004) for a very thorough discussion of mediation.  The mediate function will do some basic mediation and moderation models, with bootstrapped confidence intervals for the mediation/moderation effects. 

Functionally, this is just regular linear regression and partial correlation with some different output.

In the case of two predictor variables, X and M, and a criterion variable Y, then the direct effect of X on Y, labeled with the path c, is said to be mediated by the effect of x on M (path a) and the effect of M on Y (path b).  This partial effect (a b) is said to mediate the direct effect of X --c--> Y:     X --a ->  M  --b--> Y with X --c'--> Y where c' = c - ab.

Testing the significance of the ab mediation effect is done through bootstrapping many random resamples (with replacement) of the data.  

For moderation, the moderation effect of Z on the relationship between X -> Y is found by taking the (centered) product of X and Z and then adding this XZ term into the regression.

In the case of being provided just a correlation matrix, the bootstrapped values are based upon bootstrapping from data matching the original covariance/correlation matrix with the addition of normal errors.  This allows us to test the mediation/moderation effect even if not given raw data.

The function has been tested against some of the basic cases and examples in Hayes (2013) and the associated data sets.

Unless there is a temporal component that allows one to directly distinguish causal paths (time does not reverse direction), interpreting mediation models is problematic. Some people find it useful to compare the differences between mediation models where the causal paths (arrows) are reversed.  This is a mistake  and should not be done (Thoemmes, 2015). 

For fine tuning the size of the graphic output, xlim and ylim can be specified in the mediate.diagram function. Otherwise, the graphics produced by mediate and moderate use the default xlim and ylim values.
}
\value{
\item{total}{The total direct effect of x on y  (c)}
\item{direct}{The beta effects of x (c') and m  (b) on y   }
\item{indirect}{The indirect effect of x through m on y (c-ab)}
\item{mean.boot}{mean bootstrapped value of indirect effect}
\item{sd.boot}{Standard deviation of bootstrapped values}
\item{ci.quant}{The upper and lower confidence intervals based upon the quantiles of the bootstrapped distribution.}
\item{boot}{The bootstrapped values themselves.}
\item{a}{The effect of x on m}
\item{b}{The effect of m on y}
\item{b.int}{The interaction of x and mod (if specified)}
}
\references{

Hayes, Andrew F. (2013)  Introduction to mediation, moderation, and conditional process analysis: A regression-based approach.  Guilford Press. 

Preacher, Kristopher J and Hayes, Andrew F (2004) SPSS and SAS procedures for estimating indirect effects in simple mediation models. Behavior Research Methods, Instruments, \& Computers 36, (4) 717-731.

Thoemmes, Felix (2015) Reversing arrows in mediation models does not distinguish plausible models.  Basic and applied social psychology, 27: 226-234. 

Data from  Hayes (2013), Preacher and Hayes (2004), and from Kerchoff (1974)

}
\author{William Revelle
}
\note{ There are a number of other packages that do mediation analysis (e.g., sem and lavaan) and they are probably preferred.  This function is supplied for the more basic cases, with 1..k y variables, 1..n x variables, and 1 ..j mediators.  It will not do two step mediation.
}


\seealso{\code{\link{setCor}} and \code{\link{setCor.diagram}}

}
\examples{
#data from Preacher and Hayes (2004)
sobel <- structure(list(SATIS = c(-0.59, 1.3, 0.02, 0.01, 0.79, -0.35, 
-0.03, 1.75, -0.8, -1.2, -1.27, 0.7, -1.59, 0.68, -0.39, 1.33, 
-1.59, 1.34, 0.1, 0.05, 0.66, 0.56, 0.85, 0.88, 0.14, -0.72, 
0.84, -1.13, -0.13, 0.2), THERAPY = structure(c(0, 1, 1, 0, 1, 
1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 
1, 1, 1, 0), value.labels = structure(c(1, 0), .Names = c("cognitive", 
"standard"))), ATTRIB = c(-1.17, 0.04, 0.58, -0.23, 0.62, -0.26, 
-0.28, 0.52, 0.34, -0.09, -1.09, 1.05, -1.84, -0.95, 0.15, 0.07, 
-0.1, 2.35, 0.75, 0.49, 0.67, 1.21, 0.31, 1.97, -0.94, 0.11, 
-0.54, -0.23, 0.05, -1.07)), .Names = c("SATIS", "THERAPY", "ATTRIB"
), row.names = c(NA, -30L), class = "data.frame", variable.labels = structure(c("Satisfaction", 
"Therapy", "Attributional Positivity"), .Names = c("SATIS", "THERAPY", 
"ATTRIB")))
 #n.iter set to 50 (instead of default of 5000) for speed of example
mediate(1,2,3,sobel,n.iter=50)  #The example in Preacher and Hayes

#the pmi covariance matrix from Hayes. 2013.
#data set from Hayes, 2013 has 123 cases instead of the covariance matrix used here 

C.pmi <- structure(c(0.251232840197254, 0.119718779155005, 0.157470345195255, 
0.124533519925363, 0.03052112488338, 0.0734039717446355, 0.119718779155005, 
1.74573503931761, 0.647207783553245, 0.914575836332134, 0.0133613221378115, 
-0.0379181660669066, 0.157470345195255, 0.647207783553245, 3.01572704251633, 
1.25128282020525, -0.0224576835932294, 0.73973743835799, 0.124533519925363, 
0.914575836332134, 1.25128282020525, 2.40342196454751, -0.0106624017059843, 
-0.752990470478475, 0.03052112488338, 0.0133613221378115, -0.0224576835932294, 
-0.0106624017059843, 0.229241636678662, 0.884479541516727, 0.0734039717446355, 
-0.0379181660669066, 0.73973743835799, -0.752990470478475, 0.884479541516727, 
33.6509729441557), .Dim = c(6L, 6L), .Dimnames = list(c("cond", 
"pmi", "import", "reaction", "gender", "age"), c("cond", "pmi", 
"import", "reaction", "gender", "age")))

 #n.iter set to 50 (instead of default of 5000) for speed of example
mediate(y="reaction",x = "cond",m=c("pmi","import"),data=C.pmi,n.obs=123,n.iter=50)



#Data from sem package taken from Kerckhoff (and in turn, from Lisrel manual)
R.kerch <- structure(list(Intelligence = c(1, -0.1, 0.277, 0.25, 0.572, 
0.489, 0.335), Siblings = c(-0.1, 1, -0.152, -0.108, -0.105, 
-0.213, -0.153), FatherEd = c(0.277, -0.152, 1, 0.611, 0.294, 
0.446, 0.303), FatherOcc = c(0.25, -0.108, 0.611, 1, 0.248, 0.41, 
0.331), Grades = c(0.572, -0.105, 0.294, 0.248, 1, 0.597, 0.478
), EducExp = c(0.489, -0.213, 0.446, 0.41, 0.597, 1, 0.651), 
    OccupAsp = c(0.335, -0.153, 0.303, 0.331, 0.478, 0.651, 1
    )), .Names = c("Intelligence", "Siblings", "FatherEd", "FatherOcc", 
"Grades", "EducExp", "OccupAsp"), class = "data.frame", row.names = c("Intelligence", 
"Siblings", "FatherEd", "FatherOcc", "Grades", "EducExp", "OccupAsp"
))

 #n.iter set to 50 (instead of default of 5000) for speed of demo
mod.k <- mediate("OccupAsp","Intelligence",m= c(2:5),data=R.kerch,n.obs=767,n.iter=50)
mediate.diagram(mod.k) 
#print the path values 
mod.k

#Compare the following solution to the path coefficients found by the sem package

mod.k2 <- mediate(y="OccupAsp",x=c("Intelligence","Siblings","FatherEd","FatherOcc"),
     m= c(5:6),data=R.kerch,n.obs=767,n.iter=50)
mediate.diagram(mod.k2,show.c=FALSE) #simpler output 

#print the path values
mod.k2

}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
