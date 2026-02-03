\name{fsi}
\alias{fsi}

\title{
Factor Score Indeterminacy estimates )
}
\description{
Factor analysis either using Exploratory Factor Analysis \code{\link{fa}} or 
Confirmatory Factor analysis \code{\link{CFA}} estimates factor loadings to best fit a latent variable model that is underdetermined.  An infinite number of  factor scores from such a model all will fit equally well (Grice, 2001, Nicewander, 2020, Waller, 2023). \code{\link{fsi}} calculates an estimated R2 of observed factor scores with latent factors using either the Nicewander approach (default) or Grice approach (to match previous versions.)}
\usage{
fsi(f, phi = NULL, r = NULL, Grice = FALSE,short=TRUE)
}
\arguments{
  \item{f}{A factor loading matrix}
  \item{phi}{Factor correlations}
  \item{r}{The correlations of the raw data -- needed for the Grice option}
  \item{Grice}{If TRUE, find the Grice (2001) solution.}
 \item{short}{By default, return a vector of R2 values for each factor. Alternatively report the R2 values and H matrix (see Nicewander, 2020).}
}
\details{
A dirty little secret of factor analysis is that the factor scores are not perfectly correlated with their corresponding factors. In fact, there are an infinite number of factor score solutions that all fit equally well. The factor scores are, in fact, indeterminant.  Although known for almost 100 years, the lack of determinancy is rarely reported. ``for any set of p, imperfectly measured variables that fit Spearman's model, a research can construct an infinite set of factor scores, each of which fits the model perfectly." Waller,2023, p 245. 
 

Reviews of the factor score inderminancy (FSI)  problem and possible solutions include Grice (2001), Nicewander (2020) and Waller (2023).  

The Grice (2001)  approach is highly cited, but recent work by Nicewander (2020) and Waller (2023) suggests an alternative approach, which is more highly grounded in statistical theory.  \code{\link{fsi}} will report both of these solutions. Nicewander approach is used in the \code{\link{CFA}} function and has replaced (as of version 2.6.1) the Grice approach used in \code{\link{fa}}. 

The \code{\link{fsi}} will report the R2 between putative factor scores and the latent factors using either the Nicewander (2020) or the Grice (2001) approach.  Prior versions of \code{\link{fa}} in psych reported just the Grice solution.  As of psych_2.6.1 the default for  \code{\link{fa}} is to report the Nicewander estimate, but with an option to use the Grice estimate.  
}
\value{
%%  If it is a LIST, use
 \item{R2 }{Muliple R2 of the factor model with the factors}
 \item{H}{The first nfactor elments of the H matrix of factor correlations in the nvar space.  See Nicewander, 2020.}

}
\references{
Grice, James W. (2001) Computing and evaluating factors scores. Psychological Methods, 6, p 430-450.  10.1037/1082-989X.6.4.430 

Nicewander, W. Alan  (2020) A pespective on the mathematical and psychometric aspects of factor indeterminacy.  Multivariate Behavioral Research. 55  DOI 10.1080/00273171.2019.1684872) 

Waller, N. G. (2023) Breaking our silence on factor score indeterminacy.  Journal of Educational and Behavioral Statistics. 48, p 244-261.  doi; 10.3102/1076998622112.
 

}
\author{William Revelle}



\seealso{\code{\link{fa}} and \code{\link{CFA}}. }

\examples{
#matrix from Waller example for fsIndeterminacy from fungible
Lambda <- matrix(c(.8,  0,
                   .7,  0,
                   .6,  0,
                    0, .5,
                    0, .4,
                    0, .3), 6, 2, byrow=TRUE)
                       
 fsi(Lambda)  #should be .7675 and .3837
 
 A <- matrix( c(.8, 0,
                .8, 0,
                 0,  .6,
                 0,  .6), 4,2, byrow=TRUE)
                 
 fsi(A)   #should be .781 .529  (from Nicewander Appendix 1)                      
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ multivariate }
\keyword{ models }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
