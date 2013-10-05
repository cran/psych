\name{sim.anova}
\Rdversion{1.1}
\alias{sim.anova}

\title{Simulate a 3 way balanced ANOVA or linear model, with or without repeated measures.  }
\description{For teaching basic statistics, it is useful to be able to generate examples suitable for analysis of variance or simple linear models.  sim.anova will generate the design matrix of three independent variables (IV1, IV2, IV3) with an arbitrary number of levels and effect sizes for each main effect and interaction.  IVs can be either continuous or categorical and can have linear or quadratic effects. Either a single dependent variable or multiple (within subject) dependent variables are generated according to the specified model. The repeated measures are assumed to be tau equivalent with a specified reliability.
}
\usage{
sim.anova(es1 = 0, es2 = 0, es3 = 0, es12 = 0, es13 = 0,
    es23 = 0, es123 = 0, es11=0,es22=0, es33=0,n = 2,n1 = 2, n2 = 2, n3 = 2, 
    within=NULL,r=.8,factors=TRUE,center = TRUE,std=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 
  \item{es1}{Effect size of IV1}
  \item{es2}{Effect size of IV2}
  \item{es3}{Effect size of IV3}
  \item{es12}{Effect size of the IV1 x IV2 interaction}
  \item{es13}{Effect size of the IV1 x IV3 interaction}
  \item{es23}{Effect size of the IV2 x IV3 interaction}
  \item{es123}{Effect size of the IV1 x IV2 * IV3  interaction}
  \item{es11}{Effect size of the quadratric term of IV1}
  \item{es22}{Effect size of the quadratric term of IV2}
  \item{es33}{Effect size of the quadratric term of IV3}
  \item{n}{Sample size per cell (if all variables are categorical) or (if at least one variable is continuous), the total sample size}
  \item{n1}{Number of levels of IV1 (0) if continuous}
  \item{n2}{Number of levels of IV2}
  \item{n3}{Number of levels of IV3}
  \item{within}{if not NULL, then within should be a vector of the means of  any repeated measures.}
  \item{r}{the correlation between the repeated measures (if they exist).  This can be thought of as the reliablility of the measures.}
  \item{factors}{report the IVs as factors rather than numeric}
  \item{center}{center=TRUE provides orthogonal contrasts, center=FALSE adds the minimum value + 1 to all contrasts}
  \item{std}{Standardize the effect sizes by standardizing the IVs}
}
\details{A simple simulation for teaching about ANOVA, regression and reliability.  A variety of demonstrations of the relation between anova and lm can be shown.

The default is to produce categorical IVs (factors).  For more than two levels of an IV, this will show the difference between the linear model and anova in terms of the comparisons made.

The within vector can be used to add congenerically equivalent dependent variables. These will have intercorrelations (reliabilities) of r and means as specified as values of within.

To demonstrate the effect of centered versus non-centering, make factors = center=FALSE. The default is to center the IVs. By not centering them, the lower order effects will  be incorrect given the higher order interaction terms. 


}
\value{y.df is a data.frame of the 3 IV values as well as the DV values. 
 \item{IV1 ... IV3}{Independent variables 1 ... 3}
 \item{DV}{If there is a single dependent variable}
 \item{DV.1 ... DV.n}{If within is specified, then the n within subject dependent variables}
}

\author{William Revelle}

\seealso{ The general set of simulation functions in the psych package \code{\link{sim}}

}
\examples{
set.seed(42)
data.df <- sim.anova(es1=1,es2=.5,es13=1)  # one main effect and one interaction
describe(data.df)
pairs.panels(data.df)   #show how the design variables are orthogonal
#
summary(lm(DV~IV1*IV2*IV3,data=data.df))
summary(aov(DV~IV1*IV2*IV3,data=data.df))
set.seed(42)
 #demonstrate the effect of not centering the data on the regression
data.df <- sim.anova(es1=1,es2=.5,es13=1,center=FALSE)  #
describe(data.df)
#
#this one is incorrect, because the IVs are not centered
summary(lm(DV~IV1*IV2*IV3,data=data.df)) 

summary(aov(DV~IV1*IV2*IV3,data=data.df)) #compare with the lm model
#now examine multiple levels and quadratic terms
set.seed(42)
data.df <- sim.anova(es1=1,es13=1,n2=3,n3=4,es22=1)
summary(lm(DV~IV1*IV2*IV3,data=data.df))
summary(aov(DV~IV1*IV2*IV3,data=data.df))
pairs.panels(data.df)
#
data.df <- sim.anova(es1=1,es2=-.5,within=c(-1,0,1),n=10)
pairs.panels(data.df)

}

\keyword{models}
\keyword{multivariate}% __ONLY ONE__ keyword per line
