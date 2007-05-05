\name{correct.cor}
\alias{correct.cor}

\title{ Find dis-attenuated correlations and give alpha reliabilities }
\description{Given a raw correlation matrix and a vector of reliabilities, report the disattenuated correlations above the diagonal.
}
\usage{
correct.cor(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A raw correlation matrix }
  \item{y}{ Vector of reliabilities }
}
\details{Disattenuated correlations may be thought of as  correlations between the latent variables measured by a set of observed variables. That is, what would the correlation be between two (unreliable) variables be if both variables were measured perfectly reliably.

Examples of the output of this function are seen in \code{\link{cluster.loadings}} and \code{\link{cluster.cor}}
}
\value{Raw correlations below the diagonal, reliabilities on the diagonal, disattenuated above the diagonal.
}
\references{ http://personality-project.org/revelle/syllabi/405.syllabus.html }
\author{
Maintainer: William Revelle \email{revelle@northwestern.edu}
}
\seealso{ \code{\link{cluster.loadings}} and \code{\link{cluster.cor}}}
\examples{

# attitude from the datasets package
#example 1 is a rather clunky way of doing things
\dontrun{
a1 <- attitude[,c(1:3)]
a2 <- attitude[,c(4:7)]
x1 <- rowSums(a1)  #find the sum of the first 3 attitudes
x2 <- rowSums(a2)   #find the sum of the last 4 attitudes
alpha1 <- alpha.scale(x1,a1)
alpha2 <- alpha.scale(x2,a2)
x <- matrix(c(x1,x2),ncol=2)
x.cor <- cor(x)
alpha <- c(alpha1,alpha2)
round(correct.cor(x.cor,alpha),2)
#
#much better - although uses standardized alpha 
clusters <- matrix(c(rep(1,3),rep(0,7),rep(1,4)),ncol=2)
cluster.loadings(clusters,cor(attitude))
# or 
clusters <- matrix(c(rep(1,3),rep(0,7),rep(1,4)),ncol=2)
cluster.cor(clusters,cor(attitude))
#
}


## The function is currently defined as
"correct.cor" <-
function(x,y) { n=dim(x)[1]   
        { diag(x) <- y 
        if (n> 1)  {
        for (i in 2:n) {
           k=i-1
           for (j in 1:k) {
              x[j,i] <- x[j,i]/sqrt(y[i]*y[j])  }   #fix the upper triangular part of the matrix
             }}
           return(x)  }}
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
