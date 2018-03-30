\name{factor2cluster}
\alias{factor2cluster}
\title{ Extract cluster definitions from factor loadings }
\description{Given a factor or principal components loading matrix, assign each item to a cluster corresponding to the largest (signed) factor loading for that item.  Essentially, this is a Very Simple Structure approach to cluster definition that corresponds to what most people actually do: highlight the largest loading for each item and ignore the rest.
}
\usage{
factor2cluster(loads, cut = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loads}{either a matrix of loadings, or the result of a factor analysis/principal components analyis with a loading component }
  \item{cut}{Extract items with absolute loadings > cut}
}
\details{A factor/principal components analysis loading matrix is converted to a cluster (-1,0,1) definition matrix where each item is assigned to one and only one cluster.  This is a fast way to extract items that will be unit weighted to form cluster composites.  Use this function in combination with cluster.cor to find the corrleations of these composite scores. 

 A typical use in the SAPA project is to form item composites by clustering or factoring (see \code{\link{ICLUST}}, \code{\link{principal}}), extract the clusters from these results (\code{\link{factor2cluster}}), and then form the composite correlation matrix using \code{\link{cluster.cor}}.  The variables in this reduced matrix may then be used in multiple R procedures using mat.regress.
 
 The input may be a matrix of item loadings, or the output from a factor analysis which includes a loadings matrix.
}
\value{a matrix of -1,0,1 cluster definitions for each item.
}
\references{ \url{https://personality-project.org/r/r.vss.html}  }
\author{ 
 \url{https://personality-project.org/revelle.html} \cr

Maintainer: William Revelle  \email{  revelle@northwestern.edu           }
}

\seealso{\code{\link{cluster.cor}}, \code{\link{factor2cluster}}, \code{\link{fa}}, \code{\link{principal}}, \code{\link{ICLUST}} }
\examples{

\dontrun{
f  <- factanal(x,4,covmat=Harman74.cor$cov)
factor2cluster(f) }
#                       Factor1 Factor2 Factor3 Factor4
#VisualPerception             0       1       0       0
#Cubes                        0       1       0       0
#PaperFormBoard               0       1       0       0
#Flags                        0       1       0       0
#GeneralInformation           1       0       0       0
#PargraphComprehension        1       0       0       0
#SentenceCompletion           1       0       0       0
#WordClassification           1       0       0       0
#WordMeaning                  1       0       0       0
#Addition                     0       0       1       0
#Code                         0       0       1       0
#CountingDots                 0       0       1       0
#StraightCurvedCapitals       0       0       1       0
#WordRecognition              0       0       0       1
#NumberRecognition            0       0       0       1
#FigureRecognition            0       0       0       1
#ObjectNumber                 0       0       0       1
#NumberFigure                 0       0       0       1
#FigureWord                   0       0       0       1
#Deduction                    0       1       0       0
#NumericalPuzzles             0       0       1       0
#ProblemReasoning             0       1       0       0
#SeriesCompletion             0       1       0       0
#ArithmeticProblems           0       0       1       0





}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
