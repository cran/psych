\name{block.random}
\alias{block.random}
\title{Create a block randomized structure for n independent variables}

\description{Random assignment of n subjects with an equal number in all of N conditions may  done by block randomization, where the block size is the number of experimental conditions. The number of Independent Variables and the number of levels in each IV are specified as input. The output is a the block randomized design.
}
\usage{
block.random(n, ncond = NULL)
}

\arguments{
  \item{n}{The number of subjects to randomize.  Must be a multiple of the number of experimental conditions}
\item{ncond}{The number of conditions for each IV.  Defaults to 2 levels for one IV.  If more than one IV, specify as a vector. If names are provided, they are used, otherwise the IVs are labeled as IV1 ... IVn}
}

\value{
\item{blocks }{A matrix of subject numbers, block number, and randomized levels for each IV}
}

\author{William Revelle
}
\note{Prepared for a course on Research Methods in Psychology \url{https://personality-project.org/revelle/syllabi/205/205.syllabus.html}
}


\examples{
br <- block.random(n=24,c(2,3))
pairs.panels(br)
br <- block.random(96,c(time=4,drug=3,sex=2))
pairs.panels(br)
}
\keyword{multivariate }

