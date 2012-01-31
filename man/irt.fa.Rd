\name{irt.fa}
\alias{irt.fa}
\alias{irt.select}

\title{Item Response Analysis by Exploratory Factor Analysis of tetrachoric/polychoric correlations}
\description{
Although exploratory factor analysis and Item Response Theory seem to be very different models of binary data, they can provide equivalent parameter estimates of item difficulty and item discrimination.  Tetrachoric or polychoric correlations of a data set of dichotomous or polytomous items may be factor analysed using a minimum residual or maximum likelihood factor analysis and the result loadings transformed to item discrimination parameters.  The tau parameter from the tetrachoric/polychoric correlations combined with the item factor loading may be used to estimate item difficulties. 
}
\usage{
irt.fa(x,nfactors=1,correct=TRUE,plot=TRUE,...)
irt.select(x,y)
}
\arguments{
  \item{x}{A data matrix of dichotomous or discrete items, or the result of \code{\link{tetrachoric}} or \code{\link{polychoric}}   }
  \item{nfactors}{Defaults to 1 factor}
  \item{correct}{If true, then correct the tetrachoric correlations for continuity.  (See \code{\link{tetrachoric}}).  }
  \item{plot}{If TRUE, automatically call the \code{\link{plot.irt}} or \code{\link{plot.poly}} functions.}
  \item{y}{the subset of variables to pick from the rho and tau output of a previous irt.fa analysis to allow for further analysis.}
  \item{...}{Additional parameters to pass to the factor analysis function}
}
\details{
\code{\link{irt.fa}} combines several functions into one to make the process of item response analysis easier.  Correlations are found using either \code{\link{tetrachoric}} or \code{\link{polychoric}}.  Exploratory factor analyeses with all the normal options are then done using \code{\link{fa}}.  The results are then organized to be reported in terms of IRT parameters (difficulties and discriminations) as well as the more conventional factor analysis output. In addition, because the correlation step is somewhat slow, reanalyses may be done using the correlation matrix found in the first step.


The tetrachoric correlation matrix of dichotomous items may be factored using a (e.g.) minimum residual factor analysis function \code{\link{fa}} and the resulting loadings, \eqn{\lambda_i} are transformed to discriminations by
\eqn{\alpha = \frac{\lambda_i}{\sqrt{1-\lambda_i^2}} }{a = \lambda / (sqrt(1-\lambda^2)}.

The difficulty parameter, \eqn{\delta} is found from the \eqn{\tau} parameter of the \code{\link{tetrachoric}} or \code{\link{polychoric}} function.

\eqn{\delta_i = \frac{\tau_i}{\sqrt{1-\lambda_i^2}}}{\delta = \tau / (sqrt(1-\lambda^2)}

Similar analyses may be done with discrete item responses using polychoric correlations and distinct estimates of item difficulty (location)  for each item response.

The results may be shown graphically using \code{link{plot.irt}} (which may be called by plotting the irt.fa output, see the examples).   For plotting there are three options: type = "ICC" will plot the item characteristic response function.  type = "IIC" will plot the item information function, and type= "test" will plot the test information function.  Invisible output from the plot function will return tables of item information as a function of several levels of the trait, as well as the standard error of measurement and the reliability at each of those levels.

The normal input is just the raw data.  If, however, the correlation matrix has already been found using \code{\link{tetrachoric}}, \code{\link{polychoric}}, or a previous analysis using \code{\link{irt.fa}} then that result can be processed directly.  Because  \code{\link{irt.fa}} saves the rho and tau matrices from the analysis, subsequent analyses of the same data set are much faster if the input is the object returned on the first run.  A similar feature is available in \code{\link{omega}}. 

The output is best seen in terms of graphic displays.  Plot the output from irt.fa to see item and test information functions.  

The print function will print the item location and discriminations.  The additional factor analysis output is available as an object in the output and may be printed directly by specifying the $fa object.

The \code{\link{irt.select}} function is a helper function to allow for selecting a subset of a prior analysis for further analysis. First run irt.fa, then select a subset of variables to be analyzed in a subsequent irt.fa analysis.  Perhaps a better approach is to just plot and find the information for selected items.  

The plot function for an irt.fa object will plot ICC (item characteristic curves), IIC (item information curves), or test information curves. In addition, by using the "keys" option,  these three kinds of plots can be done for selected items. This is particularly useful when trying to see the information characteristics of short forms of tests based upon the longer form factor analysis.

The plot function will also return (invisibly) the informaton at multiple levels of the trait, the average information (area under the curve) as well as the location of the peak information for each item.  These may be then printed or printed in sorted order using the sort option in print.
}
\value{
\item{irt }{A list of Item location (difficulty) and discrimination}
\item{fa }{A list of statistics for the factor analyis}
\item{rho}{The tetrachoric/polychoric correlation matrix}
\item{tau}{The tetrachoric/polychoric cut points}

}
\references{
Kamata, Akihito and Bauer, Daniel J. (2008) A Note on the Relation Between Factor Analytic and Item Response Theory Models
Structural Equation Modeling, 15 (1) 136-153.
	
McDonald, Roderick P. (1999) Test theory: A unified treatment. L. Erlbaum Associates.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 
}
\author{William Revelle}
\note{

In comparing irt.fa to the ltm function in the ltm package or to the analysis reported in Kamata and Bauer (2008) the discrimination parameters are not identical, because the irt.fa reports them in units of the normal curve while ltm and Kamata and Bauer report them in logistic units.  In addition, Kamata and Bauer do their factor analysis using a logistic error model.  Their results match the irt.fa results (to the 2nd or 3rd decimal) when examining their analyses using a normal model.  (With thanks to Akihito Kamata for sharing that analysis.) 

\code{\link{irt.fa}} reports parameters in normal units.  To convert them to conventional IRT parameters, multiply by 1.702.  In addition, the location parameter is expressed in terms of difficulty (high positive scores imply lower frequency of response.)

The results of \code{\link{irt.fa}}  can be used by \code{\link{score.irt}} for irt based scoring.  First run \code{\link{irt.fa}} and then score the results using a two parameter model using \code{\link{score.irt}}.  

}

\seealso{\code{\link{fa}}, \code{\link{sim.irt}}, \code{\link{tetrachoric}}, \code{\link{polychoric}} as well as \code{\link{plot.psych}} for plotting the IRT item curves.

See also \code{\link{score.irt}} for scoring items based upon these parameter estimates.

}
\examples{

set.seed(17)
d9 <- sim.irt(9,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
test 
op <- par(mfrow=c(3,1))
plot(test,type="ICC")
plot(test,type="IIC")
plot(test,type="test")
par(op)
set.seed(17)
items <- sim.congeneric(N=500,short=FALSE,categorical=TRUE) #500 responses to 4 discrete items
d4 <- irt.fa(items$observed)  #item response analysis of congeneric measures
d4    #show just the irt output
d4$fa  #show just the factor analysis output


op <- par(mfrow=c(2,2))
plot(d4,type="ICC")
par(op)


#using the iq data set for an example of real items
#first need to convert the responses to tf
data(iqitems)
iq.keys <- c(4,4,3,1,4,3,2,3,1,4,1,3,4,3)
iq.tf <- score.multiple.choice(iq.keys,iqitems,score=FALSE)  #just the responses
iq.irt <- irt.fa(iq.tf)
print(iq.irt,short=FALSE) #show the IRT as well as factor analysis output
p.iq <- plot(iq.irt)  #save the invisible summary table
p.iq  #show the summary table of information by ability level
#select a subset of these variables
small.iq.irt <- irt.select(iq.irt,c(1,5,9,10,11,13))
small.irt <- irt.fa(small.iq.irt)
plot(small.irt)
#find the information for three subset of iq items
keys <- make.keys(14,list(all=1:14,some=c(1,5,9,10,11,13),others=c(1:5)))
plot(iq.irt,keys=keys)

#compare output to the ltm package or Kamata and Bauer   -- these are in logistic units 
ls <- irt.fa(lsat6)
#library(ltm)
# lsat.ltm <- ltm(lsat6~z1)
#  round(coefficients(lsat.ltm)/1.702,3)  #convert to normal (approximation)
#
#   Dffclt Dscrmn
#Q1 -1.974  0.485
#Q2 -0.805  0.425
#Q3 -0.164  0.523
#Q4 -1.096  0.405
#Q5 -1.835  0.386


#Normal results  ("Standardized and Marginal")(from Akihito Kamata )       
#Item       discrim             tau 
#  1       0.4169             -1.5520   
#  2       0.4333             -0.5999 
#  3       0.5373             -0.1512 
#  4       0.4044             -0.7723  
#  5       0.3587             -1.1966
#compare to ls 

  #Normal results  ("Standardized and conditional") (from Akihito Kamata )   
#item            discrim   tau
#  1           0.3848    -1.4325  
#  2           0.3976    -0.5505 
#  3           0.4733    -0.1332 
#  4           0.3749    -0.7159 
#  5           0.3377    -1.1264 
#compare to ls$fa and ls$tau 

#Kamata and Bauer (2008) logistic estimates
#1   0.826    2.773
#2   0.723    0.990
#3   0.891    0.249  
#4   0.688    1.285
#5   0.657    2.053
 
 

 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}
