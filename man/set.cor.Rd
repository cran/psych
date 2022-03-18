\name{setCor}
\alias{setCor}
\alias{setCor.diagram}
\alias{set.cor}
\alias{mat.regress}
\alias{matReg} 
\alias{crossValidation}
\alias{matPlot}
\title{Multiple Regression and Set Correlation from matrix or raw input}
\description{Given a correlation matrix or a  matrix or dataframe of raw data, find the multiple regressions and draw a path diagram relating a set of y variables as a function of a set of x variables.  A set of covariates (z) can be partialled from the x and y sets. Regression diagrams are automatically included.   Model can be specified in conventional formula form, or in terms of x variables and y variables.  Multiplicative models (interactions) and quadratic terms may be specified in the formula mode if using raw data. By default, the data may be  zero centered before finding the interactions.  Will also find Cohen's Set Correlation between a predictor set of variables (x) and a criterion set (y). Also finds the canonical correlations between the x and y sets.
}
\usage{
setCor(y,x,data,z=NULL,n.obs=NULL,use="pairwise",std=TRUE,square=FALSE,
       main="Regression Models",plot=TRUE,show=FALSE,zero=TRUE, alpha = .05,part=FALSE)
setCor.diagram(sc,main="Regression model",digits=2,show=FALSE,cex=1,l.cex=1,...)
#an alias to setCor
set.cor(y,x,data,z=NULL,n.obs=NULL,use="pairwise",std=TRUE,square=FALSE,
 main="Regression Models",plot=TRUE,show=FALSE,zero=TRUE,part=FALSE)   
  
mat.regress(y, x,data, z=NULL,n.obs=NULL,use="pairwise",square=FALSE) #the old form

#does not handle formula input
matReg(x,y,C,m=NULL,z=NULL,n.obs=0,means=NULL,std=FALSE,raw=TRUE,part=FALSE) 


crossValidation(model,data,options=NULL,select=NULL)
matPlot(x, type = "b", minlength=6, xlas=0,legend=NULL,
      lab=NULL,pch=16,col=1:6,lty=NULL,...)
}

\arguments{
  \item{y}{Three options: 'formula' form (similar to lm) or either  the column numbers of the y set  (e.g., c(2,4,6) or the column names of the y set (e.g., c("Flags","Addition"). See notes and examples for each.}
  \item{x}{ either the column numbers of the x set (e.g., c(1,3,5) or the column names of the x set (e.g. c("Cubes","PaperFormBoard").  x and y may also be set by use of the formula style of lm.}
  \item{data}{A matrix or data.frame of correlations or, if not square, of raw data}
  \item{C}{A variance/covariance matrix, or a correlation matrix}
  \item{m}{The column name or numbers of the set of mediating variables (see \code{\link{mediate}}).}
  \item{z}{the column names or numbers of the set of covariates. }
  \item{n.obs}{If specified, then confidence intervals, etc. are calculated, not needed if raw data are given.}

  \item{use}{find the correlations using "pairwise" (default) or just use "complete" cases (to match the lm function)}
  \item{std}{Report standardized betas (based upon the correlations) or raw bs (based upon covariances)}
  \item{part}{z is specified should part (TRUE) or partial (default) correlations be found.}
  \item{raw}{Are data from a correlation matrix or data matrix?}
  \item{means}{A vector of means for the data in matReg if giving matrix input}
  \item{main}{The title for setCor.diagram}
  \item{square}{if FALSE, then square matrices are treated as correlation matrices not as data matrices. In the rare case that one has as many cases as variables, then set square=TRUE.}
  \item{sc}{The output of setCor may be used for drawing diagrams}
  \item{digits}{How many digits should be displayed in the setCor.diagram?}
  \item{show}{Show the unweighted matrix correlation between the x and y sets?}
  \item{zero}{zero center the x variables before finding the interaction terms.}
  \item{alpha}{p value of the confidence intervals for the beta coefficients}
  \item{plot}{By default, setCor makes a plot of the results, set to FALSE to suppress the plot}
  \item{cex}{Text size of boxes displaying the variables in the diagram}
  \item{l.cex}{Text size of numbers in arrows, defaults to cex}
  \item{...}{Additional graphical parameters for setCor.diagram}
  \item{model}{The resulting object from setCor or bestScales}
  \item{options}{If using a bestScales object, which set of keys to use ("best.key","weights","optimal.key","optimal.weights")}
  \item{select}{Not yet implemented option to crossValidation}
  \item{type}{type="b" draws both lines and points}
  \item{xlas}{Orientation of the x labels (3 to make vertical)}
  \item{minlength}{When abbreviating x labels how many characters as a minimum}
  \item{legend}{If not NULL, the draw a legend at the appropriate location}
  \item{pch}{What plot character(s) to use in matPlot}
  \item{col}{colors to use in matPlot}
  \item{lty}{line types to use in matPlot}
  \item{lab}{What labels should be supplied to the lines in matPlot.  (Defaults to colnames of the variables.)}
  
}
\details{
 Although it is more common to calculate multiple regression and canonical correlations from the raw data, it is,  of course, possible to do so from a matrix of correlations or covariances.  In this case, the input to the function is a square covariance or correlation matrix, as well as the column numbers (or names) of the x (predictor),  y (criterion) variables, and if desired z (covariates). The function will find the correlations if given raw data.
 
 Input is either the set of y variables and the set of x variables, this can be written in the standard formula style of lm (see last example).  In this case, pairwise  or higher  interactions (product terms) may also be specified.  By default, when finding product terms,  the predictive variables are zero centered (Cohen, Cohen, West and Aiken, 2003), although this option can be turned off (zero=FALSE) to match the results of \code{\link{lm}} or the results discussed in Hayes (2013).  
 
 Covariates to be removed are specified by a negative sign in the formula input or by using the z variable.  Note that when specifying covariates, the regressions are done as if the regressions were done on the partialled variables.  This means that the degrees of freedom and the R2 reflect the regressions of the partialled variables. (See the last example.)  
 
 If using covariates, should they be removed from the dependent as well as the independent variables (part = FALSE, the default which means partial correlations  or just the independent variables (part=TRUE or part aka semi-partial correlations).  The difference between part correlations and partial correlations is whether the variance of the covariate is removed from both the DV and IVs (a partial correlation) or from just the IVs (a part correlation).  The slopes of the regressions remain the same, but the amount of variance in the DV (and thus the standard errors) is larger when using part correlations.  Using partial correlations for the covariates is equivalent to using the covariate in the regression when interpreting the effects of the other variables. 
  
  The output is a set of multiple correlations, one for each dependent variable in the y set, as well as the set of canonical correlations.
  
An additional output is the R2 found using Cohen's set correlation (Cohen, 1982).  This is a measure of how much variance and the x and y set share.
  
 Cohen (1982) introduced the set correlation, a multivariate generalization of the multiple correlation to measure the overall relationship between two sets of variables. It is an application of canoncial correlation (Hotelling, 1936) and \eqn{1 - \prod(1-\rho_i^2)} where \eqn{\rho_i^2} is the squared canonical correlation.  Set correlation is the amount of shared variance (R2) between two sets of variables.  With the addition of a third, covariate set, set correlation will find multivariate R2, as well as partial and semi partial R2.  (The semi and bipartial options are not yet implemented.) Details on set correlation may be found in Cohen (1982), Cohen (1988) and  Cohen, Cohen, Aiken and West (2003). 
 
 R2 between two sets is just \deqn{R^2 = 1- \frac{\left | R_{yx} \right |}{\left | R_y \right | \left |R_x\right |}  =  1 - \prod(1-\rho_i^2)   }{R2 = 1- |R| /(|Ry| * |Rx|)} where R is the  complete correlation matrix of the x and y variables and Rx and Ry are the two sets involved.

 
 Unfortunately, the R2 is sensitive to one of the canonical correlations being very high.  An alternative, T2, is the proportion of additive variance and is the average of the squared canonicals.  (Cohen et al., 2003), see also Cramer and Nicewander (1979).  This  average, because it includes some very small canonical correlations, will tend to be too small.  Cohen et al. admonition is appropriate: "In the final analysis, however, analysts must be guided by their substantive and methodological conceptions of the problem at hand in their choice of a measure of association." ( p613). 
 
 Yet another measure of the association between two sets is just the simple, unweighted correlation between the two sets. That is, 
 \deqn{R_{uw}  =\frac{ 1 R_{xy} 1' }{(1R_{yy}1')^{.5} (1R_{xx}1')^{.5}}    }{Ruw=1Rxy1' / (sqrt(1Ryy1'* 1Rxx1'))} where Rxy is the matrix of correlations between the two sets.  This is just  the simple (unweighted) sums of the correlations in each matrix. This technique exemplifies the robust beauty of linear models and  is particularly appropriate in the case of one dimension in both x and y, and will be a drastic underestimate in the case of items where the betas differ in sign. 
 
 When finding the unweighted correlations, as is done in \code{\link{alpha}}, items are flipped so that they all are positively signed.  
 
 
  A typical use in the SAPA project is to form item composites by clustering or factoring (see  \code{\link{fa}},\code{\link{ICLUST}}, \code{\link{principal}}), extract the clusters from these results (\code{\link{factor2cluster}}), and then form the composite correlation matrix using \code{\link{cluster.cor}}.  The variables in this reduced matrix may then be used in multiple R procedures using \code{\link{setCor}}.

 Although the overall matrix can have missing correlations, the correlations in the subset of the matrix used for prediction must exist.
 
 If the number of observations is entered, then the conventional confidence intervals, statistical significance, and shrinkage estimates are  reported.
 
 If the input is rectangular (not square), correlations or covariances are found from the data.
 
 The print function reports t and p values for the beta weights, the summary function just reports the beta weights.
 
 The Variance Inflation Factor is reported but should be taken with the normal cautions of interpretation discussed by Guide and Ketokivi.  That is to say, VIF > 10 is not a magic cuttoff to define colinearity.  It is merely 1/(1-smc(R(x)).
 
 The Guide and Ketokivi article is well worth reading for all who want to use various regression models. 
 
 \code{\link{crossValidation}} can be used to take the results from \code{\link{setCor}} or \code{\link{bestScales}} and apply the weights to a different data set.
 
 \code{\link{matPlot}} can be used to plot the crossValidation values (just a call to matplot with the xaxis given labels). matPlot has been improved to draw legends and to allow for specifications of the line types.
 
 \code{\link{setCorLookup}} will sort the beta weights and report them with item contents if given a dictionary. 
 \code{\link{matReg}} is primarily a helper function for \code{\link{mediate}} but is a general multiple regression function given a covariance matrix and the specified x,  y and z variables. Its output includes betas, se, t, p and R2.  The call includes m for mediation variables, but these are only used to adjust the degrees of freedom.
 \code{\link{matReg}} does not work on data matrices, nor does it take formula input.  It is really just a helper function for  \code{\link{mediate}}
}
\value{
  
  \item{beta }{the beta weights for each variable in X for each variable in Y}
  \item{R}{The multiple R for each equation (the amount of change a unit in the predictor set leads to in the criterion set). }
  \item{R2 }{The multiple R2 (\% variance acounted for) for each equation}
  \item{VIF}{The Variance Inflation Factor which is just 1/(1-smc(x))}
  \item{se}{Standard errors of beta weights (if n.obs is specified)}
  \item{t}{t value of beta weights (if n.obs is specified)}
  \item{Probability}{Probability of beta = 0 (if n.obs is specified)}
  \item{shrunkenR2}{Estimated shrunken R2 (if n.obs is specified)}
  \item{setR2}{The multiple R2 of the set correlation between the x and y sets}
  item{residual}{The residual correlation matrix of Y with x and z removed}
  \item{ruw}{The unit weighted multiple correlation for each dependent variable}
  \item{Ruw}{The unit weighted set correlation} 

}

\author{William Revelle \cr
       
Maintainer: William Revelle <revelle@northwestern.edu> }
\references{
J. Cohen (1982)  Set correlation as a general multivariate data-analytic method. Multivariate Behavioral Research, 17(3):301-341.


J. Cohen, P. Cohen, S.G. West, and L.S. Aiken. (2003) Applied multiple regression/correlation analysis for the behavioral sciences. L. Erlbaum Associates, Mahwah, N.J., 3rd ed edition.

H. Hotelling. (1936) Relations between two sets of variates. Biometrika 28(3/4):321-377.


E.Cramer and W. A. Nicewander (1979) Some symmetric, invariant measures of multivariate association. Psychometrika, 44:43-54.

V. Daniel R. Guide Jr. and M. Ketokivim (2015) Notes from the Editors: Redefining some methodological criteria for the journal.  Journal of Operations Management. 37. v-viii.

}

\note{As of April 30, 2011, the order of x and y was swapped in the call to be consistent with the general y ~ x  syntax of the lm and aov functions. In addition, the primary name of the function was switched to setCor from mat.regress to reflect the estimation of the set correlation. 

In October, 2017 I added the ability to specify the input in formula mode and allow for higher level and multiple interactions. 

The denominator degrees of freedom for the set correlation does not match that reported by Cohen et al., 2003 in the example on page 621 but does match the formula on page 615, except for the typo in the estimation of F (see Cohen 1982).   The difference seems to be that they are adding in a correction factor of df 2 = df2 + df1. 

Following a suggestion by Keith Widaman, when zero centering for product terms, the y variables are no longer zero centered. This puts the regression in the units of the DV. (2/22/22). 



}

\seealso{\code{\link{mediate}} for an alternative regression model with 'mediation'.  \code{\link{predict}} to find predicted values given regression weights. 
\code{\link{cluster.cor}}, \code{\link{factor2cluster}},\code{\link{principal}},\code{\link{ICLUST}}, \code{link{cancor}} and cca in the yacca package. 
\code{\link{GSBE}} for further demonstrations of mediation and moderation. }
\examples{
#First compare to lm using data input
summary(lm(rating ~ complaints + privileges, data = attitude))
setCor(rating ~ complaints + privileges, data = attitude, std=FALSE) #do not standardize
z.attitude <- data.frame(scale(attitude))  #standardize the data before doing lm
summary(lm(rating ~ complaints + privileges, data = z.attitude))  #regressions on z scores
setCor(rating ~ complaints + privileges, data = attitude)  #by default we standardize and 
# the results are the same as the standardized lm


R <- cor(attitude) #find the correlations
#Do the regression on the correlations  
#Note that these match the regressions on the standard scores of the data
setCor(rating ~ complaints + privileges, data =R, n.obs=30)

#now, partial out learning and critical
setCor(rating ~ complaints + privileges - learning - critical, data =R, n.obs=30)
#compare with the full regression:
setCor(rating ~ complaints + privileges + learning + critical, data =R, n.obs=30)



#Canonical correlations:

#The first Kelley data set from Hotelling
kelley1 <- structure(c(1, 0.6328, 0.2412, 0.0586, 0.6328, 1, -0.0553, 0.0655, 
0.2412, -0.0553, 1, 0.4248, 0.0586, 0.0655, 0.4248, 1), .Dim = c(4L, 
4L), .Dimnames = list(c("reading.speed", "reading.power", "math.speed", 
"math.power"), c("reading.speed", "reading.power", "math.speed", 
"math.power")))
lowerMat(kelley1)
mod1 <- setCor(y = math.speed + math.power ~ reading.speed + reading.power, 
    data = kelley1, n.obs=140)
mod1$cancor
#Hotelling reports .3945 and .0688  we get  0.39450592 0.06884787

#the second Kelley data from Hotelling
kelley <- structure(list(speed = c(1, 0.4248, 0.042, 0.0215, 0.0573), power = c(0.4248, 
1, 0.1487, 0.2489, 0.2843), words = c(0.042, 0.1487, 1, 0.6693, 
0.4662), symbols = c(0.0215, 0.2489, 0.6693, 1, 0.6915), meaningless = c(0.0573, 
0.2843, 0.4662, 0.6915, 1)), .Names = c("speed", "power", "words", 
"symbols", "meaningless"), class = "data.frame", row.names = c("speed", 
"power", "words", "symbols", "meaningless"))

lowerMat(kelley)

setCor(power + speed ~ words + symbols + meaningless,data=kelley)  #formula mode
#setCor(y= 1:2,x = 3:5,data = kelley) #order of variables input

#Hotelling reports canonical correlations of .3073 and .0583  or squared correlations of
# 0.09443329 and 0.00339889 vs. our values of cancor = 0.3076 0.0593  with squared values
#of  0.0946 0.0035,

setCor(y=c(7:9),x=c(1:6),data=Thurstone,n.obs=213)  #easier to just list variable  
                                              #locations if we have long names
#now try partialling out some variables
set.cor(y=c(7:9),x=c(1:3),z=c(4:6),data=Thurstone) #compare with the previous
#compare complete print out with summary printing 
sc <- setCor(SATV + SATQ ~ gender + education,data=sat.act) # regression from raw data
sc
summary(sc)

setCor(Pedigrees ~ Sentences + Vocabulary - First.Letters - Four.Letter.Words ,
data=Thurstone)  #showing formula input with two covariates

#Do some regressions with real data (rather than correlation matrices)
setCor(reaction ~ cond + pmi + import, data = Tal.Or)

#partial out importance
setCor(reaction ~ cond + pmi - import, data = Tal.Or, main="Partial out importance")

#compare with using lm by partialling
mod1 <- lm(reaction ~ cond + pmi + import, data = Tal.Or)
reaction.import <- lm(reaction~import,data=Tal.Or)$resid
cond.import <-  lm(cond~import,data=Tal.Or)$resid
pmi.import <-   lm(pmi~import,data=Tal.Or)$resid
mod.partial <- lm(reaction.import ~ cond.import + pmi.import)
summary(mod.partial)
#lm uses raw scores, so set std = FALSE for setCor
print(setCor(y = reaction ~ cond + pmi - import, data = Tal.Or,std = FALSE,
 main = "Partial out importance"),digits=4)
#notice that the dfs of the partial approach using lm are 1 more than the setCor dfs 
 
 #Show how to find quadratic terms
sc <- setCor(reaction ~ cond + pmi + I(import^2), data = Tal.Or)
sc
#pairs.panels(sc$data) #show the SPLOM of the data

#Consider an example of a derivation and cross validation sample
set.seed(42)
ss <- sample(1:2800,1400)
model <- setCor(y=26:28,x=1:25,data=bfi[ss,],plot=FALSE)
original.fit <- crossValidation(model,bfi[ss,]) #the derivation set
cross.fit <- crossValidation(model,bfi[-ss,])  #the cross validation set
summary(original.fit)
summary(cross.fit)
predicted <- predict(model,bfi[-ss,])
cor2(predicted,bfi[-ss,26:28])

}

\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
