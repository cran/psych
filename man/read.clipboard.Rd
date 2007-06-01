\name{read.clipboard}
\alias{read.clipboard}
\alias{read.clipboard.csv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{shortcut for reading from the clipboard}
\description{input from the keyboard is easy but a bit obscure, particularly for Mac users. This is just an easier mnemonic to do so.
}
\usage{
read.clipboard(header = TRUE, ...)

#my.data <- read.clipboad()         #assumes headers and tab delimited
#my.data <- read.clipboard.csv()     #assumes heades and comma delimited
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{Does the first row have variable labels }
  \item{\dots}{ Other parameters to pass to read }
}
\details{A typical session of R might involve data stored in text files, generated on line, etc. Although it is easy to just read from a file (particularly if using file.locate(), copying from the file to the clipboard and then reading from the clipboard is also very convenient (and somewhat more intuitive to the naive user.)

Based upon a suggestion by Ken Knoblauch to the R-help listserve.
}
\value{the contents of the clipboard.
}

\author{ William Revelle}

\examples{
#my.data <- read.clipboad()
#my.data <- read.clipboard.csv()
#my.data <- read.clipboad(header=FALSE)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
