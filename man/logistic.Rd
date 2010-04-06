\name{logistic}
\alias{logistic}
\alias{logit}
\alias{logistic.grm}

\title{Logistic transform from x to p and logit transform from p to x}

\description{The logistic function (1/(1+exp(-x)) and logit function (log(p/(1-p)) are fundamental to Item Response Theory.  Although just one line functions, they are included here for ease of demonstrations and in drawing IRT models. Also included is the logistic.grm for a graded response model.}
\usage{
logistic(x,d=0, a=1,c=0, z=1)
logit(p)
logistic.grm( x,d=0,a=1.5,c=0,z=1,r=2,s=c(-1.5,-.5,.5,1.5))
}

\arguments{
  \item{x}{Any integer or real value }
  \item{d}{Item difficulty or delta parameter }
  \item{a}{The slope of the curve at x=0  is equivalent to the discrimination parameter in 2PL models or alpha parameter. Is either 1  in 1PL or  1.702 in 1PN approximations.  }
  \item{c}{Lower asymptote  =  guessing parameter in 3PL models or gamma }
  \item{z}{The upper asymptote  --- in 4PL models} 
  \item{p}{Probability to be converted to logit value}
  \item{r}{The response category for the graded response model}
  \item{s}{The response thresholds} 
}
\details{These three functions are provided as simple helper functions for demonstrations of Item Response Theory. The one parameter logistic (1PL) model is also known as the Rasch model.  It assumes items differ only in difficulty.
1PL, 2PL, 3PL and 4PL curves may be drawn by choosing the appropriate d (delta or item difficulty), a (discrimination or slope), c (gamma or  guessing) and z (zeta or upper asymptote).

logit is just the inverse of logistic.

logistic.grm will create the responses for a graded response model for the rth category where cutpoints are in s.
}
\value{
\item{p}{logistic returns the probability associated with x}
\item{x}{logit returns the real number associated with p}
}


\author{ William Revelle}

\examples{
curve(logistic(x,a=1.702),-3,3,ylab="Probability of x",main="Logistic transform of x",xlab="z score units") #logistic with a=1.702 is almost the same as pnorm 
curve(pnorm(x),add=TRUE,lty="dashed")  
curve(logistic(x),add=TRUE)
text(2,.8, expression(alpha ==1))
text(2,1.0,expression(alpha==1.7))
curve(logistic(x),-4,4,ylab="Probability of x", main = "Logistic transform of x in logit units",xlab="logits")
curve(logistic(x,d=-1),add=TRUE)
curve(logistic(x,d=1),add=TRUE)
curve(logistic(x,c=.2),add=TRUE,lty="dashed")
text(1.3,.5,"d=1")
text(.3,.5,"d=0")
text(-1.5,.5,"d=-1")
text(-3,.3,"c=.2")
#demo of graded response model
 curve(logistic.grm(x,r=1),-4,4,ylim=c(0,1),main="Five level response scale",ylab="Probability of endorsement",xlab="Latent attribute on logit scale")
 curve(logistic.grm(x,r=2),add=TRUE)
 curve(logistic.grm(x,r=3),add=TRUE)
 curve(logistic.grm(x,r=4),add=TRUE)
 curve(logistic.grm(x,r=5),add=TRUE)
 
 text(-2.,.5,1)
 text(-1.,.4,2)
 text(0,.4,3)
 text(1.,.4,4)
  text(2.,.5,5)
} 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

