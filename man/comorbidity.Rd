\name{comorbidity}
\alias{comorbidity}
\title{ Convert base rates of two diagnoses and their comorbidity into phi, Yule, and tetrachorics }
\description{In medicine and clinical psychology, diagnoses tend to be categorical (someone is depressed or not, someone has an anxiety disorder or not).  Cooccurrence  of both of these symptoms is called comorbidity.   Diagnostic categories vary in their degree of comorbidity with other diagnostic categories.  From the point of view of correlation, comorbidity is just a name applied to one cell in a four fold table.  It is thus possible to analyze comorbidity rates by considering the probability of the separate diagnoses and the probability of the joint diagnosis.  This gives the two by two table needed for a phi, Yule, or tetrachoric correlation.
}
\usage{
comorbidity(d1, d2, com, labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d1}{Proportion of diagnostic category 1}
  \item{d2}{Proportion of diganostic category 2 }
  \item{com}{Proportion of comorbidity (diagnostic category 1 and 2) }
  \item{labels}{Names of categories 1 and 2}
}

\value{
  \item{twobytwo }{The two by two table implied by the input }
  \item{phi }{Phi coefficient of the two by two table}
  \item{Yule}{Yule coefficient of the two by two table}
  \item{tetra}{Tetrachoric coefficient of the two by two table}
}

\author{William Revelle }


 
\seealso{  \code{\link{phi}},  \code{\link{Yule}} }
\examples{
comorbidity(.2,.15,.1,c("Anxiety","Depression")) 
}

\keyword{multivariate }

