\name{mardia}
\alias{mardia}
\alias{skew}
\alias{kurtosi}

\title{ Calculate univariate or multivariate (Mardia's test) skew and kurtosis for a vector, matrix, or data.frame}
\description{Find the skew and kurtosis for each variable in a data.frame or matrix.  Unlike skew and kurtosis in e1071, this calculates a different skew for each variable or column of a data.frame/matrix. mardia applies Mardia's tests for multivariate skew and kurtosis
}
\usage{
skew(x, na.rm = TRUE,type=3)
kurtosi(x, na.rm = TRUE,type=3)
mardia(x,na.rm = TRUE,plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or matrix }
  \item{na.rm}{ how to treat missing data }
  \item{type}{See the discussion in describe}.
  \item{plot}{Plot the expected normal distribution values versus the Mahalanobis distance of the subjects.}
}
\details{given a matrix or data.frame x, find the skew or kurtosis for each column (for skew and kurtosis) or the multivariate skew and kurtosis in the case of mardia.


As of version 1.2.3,when finding the skew and the kurtosis, there are three different options available.  These match the choices available in skewness and kurtosis found in the e1071 package (see Joanes and Gill (1998) for the advantages of each one). 

If we define \eqn{m_r = [\sum(X- mx)^r]/n}{m_r = [sum(X- mx)^r]/n} then 

Type 1 finds skewness and kurtosis by \eqn{g_1 = m_3/(m_2)^{3/2} } and \eqn{g_2 = m_4/(m_2)^2 -3}.  

Type 2 is \eqn{G1 = g1 * \sqrt{n *(n-1)}/(n-2)} and \eqn{G2 = (n-1)*[(n+1)g2 +6]/((n-2)(n-3))}.  

Type 3 is \eqn{b1 = [(n-1)/n]^{3/2} m_3/m_2^{3/2}} and \eqn{b2 =  [(n-1)/n]^{3/2} m_4/m_2^2)}. 
 
  For consistency with e1071 and with the Joanes and Gill, the types are now defined as above.
  
However, from revision 1.0.93 to 1.2.3, kurtosi by default gives an unbiased estimate of the kurtosis (DeCarlo, 1997). Prior versions used a different equation which produced a biased estimate.  (See the kurtosis function in the e1071 package for the distinction between these two formulae.  The default, type 1 gave what is called type 2 in e1071.  The other is their type 3.)  For comparison with previous releases, specifying type = 2 will give the old estimate.  These type numbers are now changed.  
 


}
\value{
\item{skew}{if input is a matrix or data.frame, skew is a vector of skews}
\item{kurtosi}{if input is a matrix or data.frame, kurtosi is a vector of kurtosi}
\item{bp1}{Mardia's bp1 estimate of multivariate skew}
\item{bp2}{Mardia's bp2 estimate of multivariate kurtosis}
\item{skew}{Mardia's skew statistic}
\item{small.skew}{Mardia's small sample skew statistic}
\item{p.skew}{Probability of skew}
\item{p.small}{Probability of small.skew}
\item{kurtosis}{Mardia's multivariate kurtosis statistic}
\item{p.kurtosis}{Probability of kurtosis statistic}
\item{D}{Mahalanobis distance of cases from centroid}
}

\references{

Joanes, D.N. and Gill, C.A (1998).  Comparing measures of sample skewness and kurtosis.  The Statistician, 47, 183-189.

L.DeCarlo. 1997) On the meaning and use of kurtosis, Psychological Methods, 2(3):292-307,


K.V. Mardia (1970). Measures of multivariate skewness and kurtosis with applications. Biometrika, 57(3):pp. 519-30, 1970.}

\author{William Revelle }

\note{The mean function supplies means for the columns of a data.frame, but the overall mean for a matrix.  Mean will throw a warning for non-numeric data, but colMeans stops with non-numeric data. Thus, the function uses either mean (for data frames) or colMeans (for matrices).  This is true for skew and kurtosi as well.

}
\note{ 	Probability values less than  10^-300 are set to 0.  }

\seealso{\code{\link{describe}}, \code{\link{describe.by}}, mult.norm in QuantPsyc, Kurt in QuantPsyc}


\examples{
round(skew(attitude),2)   #type 3 (default)
round(kurtosi(attitude),2)  #type 3 (default)
#for the differences between the three types of skew and kurtosis:
round(skew(attitude,type=1),2)  #type 1
round(skew(attitude,type=2),2)  #type 2 
mardia(attitude)
x <- matrix(rnorm(1000),ncol=10)
describe(x)
mardia(x)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
