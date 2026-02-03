\name{fparse}
\alias{fparse}
\alias{lavParse}

\title{
Parse and extend  formula input from a model and return the DV(s), IV(s), and associated terms.
}
\description{
Formula input from e.g., lm,  may be extended to include mediators,  quadratic and partial terms using a standard syntax. This is used by (e.g.,) \code{\link{lmCor}}, \code{\link{mediate}}, and  \code{\link{error.bars.by}}.   \code{\link{lavParse}} will take a limited number of lavaan syntax commands and return as a model matrix for \code{\link{CFA}},  \code{\link{CFA.bifactor}}, \code{\link{sim}}, and \code{\link{sim.structure}}
 }
\usage{
fparse(expr)
lavParse(model, phi=FALSE)
}
\arguments{
  \item{expr}{A legitimate expression in the form y ~ x1 ,  etc. (see details)}
  \item{model}{A lavaan style factor specification. (see examples)}
  \item{phi}{If TRUE, the model can include the factor correlations. } 
}
\value{
\item{v.mat}{A matrix of variable weights for each factor}
\item{y}{A list of elements from the left side of the formula}
\item{x}{A list of elements from the right side of the formula}
\item{m}{A list of those elements of the formula included in ()}
\item{prod}{A list of elements separated by a * sign}
\item{ex}{A list of elements marked by I() }
\item{model}{A matrix with cols for the factors, rows for the variable, values specified}
}
\details{The basic formula input given as DV1 + DV2 ~ IV1 + IV2 + (IV3) + I(IV4^2) - IV5 will be parsed to return 2 DVs (1 and 2), two normal IVs (1 and 2), a mediator (IV3) a quadratic (IV4) and a variable to be partialed (IV5). See the various examples in \code{\link{lmCor}} and \code{\link{mediate}}.

To allow rapid comparisons with sem and cfa output from lavaan, it is convenient to use  lavaan style syntax.  

\code{\link{CFA}} and\code{\link{CFA}} can specify the model using lavaan syntax as well as the item keys format.  
}


\author{William Revelle}



\examples{
fparse(DV  ~ IV1 + IV2 * IV2*IV3 + (IV4) + I(IV5^2) )
#somewhat more complicated
fparse(DV1 + DV2 ~ IV1 + IV2 + IV3*IV4 + I(IV5^2) + I(Iv6^2) + (IV7) + (IV8) - IV9)


model <- '	F1 =~ .9*V1 + .8*V2 + .7*V3 
         	F2 =~ .8 * V4 + .7*V5 + .6*V6
        	F3 =~ .9*V7 + .7*V8 + .5*V9
       		F1 ~ .3 * F2 + .6* F3
        	F2 ~ .5*F3'
lavParse(model,phi=TRUE)   #will return a 9 x 3 matrix loadings and 3x3 phi matrix for sim 
 
#by specifying the coefficients, we can run sim more readily
sim.data <-sim(model,n=500)
cx <- CFA(model,sim.data$observed) 

#The CFA function will not use the coefficients nor the phi matrix, just the path model
#So we can use the same model for simulations and CFA
#But, lavaan needs to drop the coefficients
modelg <- 'g =~ .4*V1 +.4*V2+.4*V3 +.4*V4+.4*V5+.4*V6 +.4*V7+.4*V8+.4*V9
F1 =~ .9*V1 + .8*V2 + .7*V3 
         	F2 =~ .8 * V4 + .7*V5 + .6*V6
        	F3 =~ .9*V7 + .7*V8 + .5*V9   '
  
modelg <- 'g=~V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 +V9
    	F1 =~  V1 + V2 + V3 
     	F2 =~  V4 + V5 + V6
        F3 =~ V7 + V8 + V9'
lavParse(modelg)
#model without coefficients to run in lavaan
}
\keyword{ utilities }
