\name{sim.hierarchical}
\alias{sim.hierarchical}
\alias{make.hierarchical}
\alias{sim.bonds}
\title{Create a population or sample correlation matrix, perhaps with hierarchical structure. }
\description{Create a population orthogonal or hierarchical correlation matrix from a set of factor loadings and factor intercorrelations. Samples of size n may be then be drawn from this population.  Return either the sample data, sample correlations, or population correlations.  This is used to create sample data sets for instruction and demonstration.
}
\usage{
sim.hierarchical(gload=NULL, fload=NULL, n = 0, raw = TRUE,mu = NULL,
    categorical=FALSE, low=-3,high=3,threshold=NULL)
sim.bonds(nvar=9,nf=NULL, loads=c(0,0,.5,.6),validity=.8)

make.hierarchical(gload=NULL, fload=NULL, n = 0, raw = FALSE)  #deprecated

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gload}{ Loadings of group factors on a general factor. Defaults to c(.9,.8,.7)  }
  \item{fload}{ Loadings of items on the group factors. Defaults to matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)}
  \item{n}{ Number of subjects to generate: N=0 => population values }
  \item{raw}{ raw=TRUE, report the raw data, raw=FALSE, report the sample  correlation matrix. }
  \item{mu}{means for the individual variables}
  \item{low}{lower cutoff for categorical data}
  \item{categorical}{If True, then create categorical data}
  \item{threshold}{If categorical is TRUE, and binary output is desired, what is the threshold to convert continuous scores into 0/1.  May be a vector.}
  \item{high}{Upper cuttoff for categorical data}
  \item{nvar}{Number of variables to simulate}
  \item{loads}{A vector of loadings that will be sampled (rowwise) to define the factors}
  \item{validity}{The factor loadings of `pure' measures of the factor.}
  \item{nf}{Number of factors to generate in sim.bonds}
}
\details{Many personality and cognitive tests have a hierarchical factor structure.  For demonstration purposes, it is useful to be able to create such matrices, either with population values, or sample values. 

Given a matrix of item factor loadings (fload) and of loadings of these factors on a general factor (gload), we create a population correlation matrix by using the general factor law  \eqn{ R \approx F' \theta F + U^2} where \eqn{\theta = g'g}  


The default is to return population correlation matrices. Sample correlation matrices are generated if n > 0.  Raw data are returned if raw = TRUE.

The default values for gload and fload create a data matrix discussed by Jensen and Weng, 1994.

In order to properly simulate polytomous items, the categorical option will round continuous scores into integers.  To simulate dichotomous items, the threshold vector may be used to specify the value at which the continuous values are cut into 0/1.  If the length of threshold is less than the number of variables, the vector will be sampled (with replacement) to fill it out.

Although written to create hierarchical structures, if the gload matrix is all 0, then a non-hierarchical structure will be generated.

Yet another model is that of Godfrey H. Thomson (1916) who suggested that independent bonds could produce the same factor structure as a g factor model. This is simulated in \code{\link{sim.bonds}}.  Compare the \code{\link{omega}} solutions for a \code{\link{sim.hierarchical}} with a \code{\link{sim.bonds}} model. Both produce reasonable values of omega, although the one was generated without a general factor.
}
\value{
\item{r}{a matrix of correlations}
\item{model}{The population correlation matrix}
\item{observed}{The simulated data matrix with the defined structure}
\item{theta}{The latent factor scores used to generate the data. Compare how these correlate with the observed data  with the results from \code{\link{omega}}.}
\item{sl}{The Schmid Leiman transformed factor loadings.  These may be used to test factor scoring problem. }
}
\references{ \url{https://personality-project.org/r/r.omega.html }
\cr Jensen, A.R., Weng, L.J. (1994) What is a Good g? Intelligence, 18, 231-258.

 Godfrey H. Thomson (1916) A hierarchy without a general factor, British Journal of Psychology, 8, 271-281. 

}
\author{ William Revelle }

\seealso{ \code{\link{omega}},   \code{\link{schmid}},  \code{\link{ICLUST}},  \code{\link{VSS}} for ways of analyzing these data.  Also see \code{\link{sim.structure}} to simulate a variety of structural models (e.g., multiple correlated factor models).  }

\examples{

gload <-  gload<-matrix(c(.9,.8,.7),nrow=3)    # a higher order factor matrix
fload <-matrix(c(                    #a lower order (oblique) factor matrix
           .8,0,0,
           .7,0,.0,
           .6,0,.0,
            0,.7,.0,
            0,.6,.0,
            0,.5,0,
            0,0,.6,
            0,0,.5,
            0,0,.4),   ncol=3,byrow=TRUE)
            
jensen <- sim.hierarchical(gload,fload)    #the test set used by omega
round(jensen,2)
set.seed(42) #for reproducible results
jensen <-  sim.hierarchical(n=10000,categorical =TRUE, threshold =c(-1,0,1)) 
           #use the same gload and fload values, but produce the data  
#items have three levels of difficulty
#Compare factor scores using the sl model with those that generated the data
lowerCor(jensen$theta) #the correlations of the factors
fs <- factor.scores(jensen$observed, jensen$sl)  #find factor scores from the data
lowerCor(fs$scores) #these are now correlated
cor2(fs$scores,jensen$theta)  #correlation with the generating factors 


#compare this to a simulation of the bonds model
set.seed(42)
R <- sim.bonds()
R$R    

#simulate a non-hierarchical structure
fload <- matrix(c(c(c(.9,.8,.7,.6),rep(0,20)),c(c(.9,.8,.7,.6),rep(0,20)),
    c(c(.9,.8,.7,.6),rep(0,20)),c(c(c(.9,.8,.7,.6),rep(0,20)),c(.9,.8,.7,.6))),ncol=5)
gload <- matrix(rep(0,5))
five.factor <- sim.hierarchical(gload,fload,500,TRUE) #create sample data set
#do it again with a hierachical structure
gload <- matrix(rep(.7,5)  )
five.factor.g <- sim.hierarchical(gload,fload,500,TRUE) #create sample data set
#compare these two with omega
#not run
#om.5 <- omega(five.factor$observed,5)
#om.5g <- omega(five.factor.g$observed,5)
}
\keyword{ multivariate}% at least one, from doc/KEYWORDS
\keyword{models }% __ONLY ONE__ keyword per line
\keyword{datagen}
