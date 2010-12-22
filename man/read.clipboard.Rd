\name{read.clipboard}
\alias{read.clipboard}
\alias{read.clipboard.csv}
\alias{read.clipboard.tab}
\alias{read.clipboard.lower}
\alias{read.clipboard.upper}
\alias{read.clipboard.fwf}
\title{shortcut for reading from the clipboard}
\description{Input from the clipboard is easy but a bit obscure, particularly for Mac users. This is just an easier way to do so.  Data may be copied to the clipboard from Exel spreadsheets, csv files, or fixed width formatted files and then into a data.frame. Data may also be read from lower  (or upper) triangular matrices and filled out to square matrices.
}
\usage{
read.clipboard(header = TRUE, ...)   #assumes headers and tab or space delimited
read.clipboard.csv(header=TRUE,sep=',',...)   #assumes headers and comma delimited
read.clipboard.tab(header=TRUE,sep='\t',...)   #assumes headers and tab delimited
read.clipboard.lower(diag=TRUE,names=NULL,...)  #read in a matrix given the lower off diagonal
read.clipboard.upper(diag=TRUE,names=NULL,...)
read.clipboard.fwf(header=FALSE,widths=rep(1,10),...)  #read in data using a fixed format width (see read.fwf for instructions)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{Does the first row have variable labels }
  \item{sep}{What is the designated separater between data fields?}
  \item{diag}{for upper or lower triangular matrices, is the diagonal specified or not}
  \item{names}{for read.clipboard.lower or upper, what colnames to use if not defined in the first column}
  \item{widths}{how wide are the columns in fixed width input.  The default is to read 10 columns of size 1. }
  \item{\dots}{ Other parameters to pass to read }
}
\details{A typical session of R might involve data stored in text files, generated online, etc. Although it is easy to just read from a file (particularly if using file.choose(), copying from the file to the clipboard and then reading from the clipboard is also very convenient (and somewhat more intuitive to the naive user). This is particularly convenient when copying from a text book or article and just moving a section of text into R.)

Based upon a suggestion by Ken Knoblauch to the R-help listserve.

If the input file that was copied into the clipboard was an Excel file with blanks for missing data, then read.clipboard.tab() will correctly replace the blanks with NAs.  Similarly for a csv file with blank entries, read.clipboard.csv will replace empty fields with NA.  

read.clipboard.lower and read.clipboard.upper are adapted from John Fox's read.moments function in the sem package.  They will read a lower (or upper) triangular matrix from the clipboard and return a full, symmetric matrix for use by factanal, \code{\link{factor.pa}} , \code{\link{ICLUST}} , etc.  If the diagonal is false, it will be replaced by 1.0s.  These two function were added to allow easy reading of examples from various texts and manuscripts with just triangular output. 

Many articles will report lower triangular matrices with variable labels in the first column.  read.clipboard.lower (or read.clipboard.upper) will handle this case as well.

read.clipboard.fwf will read fixed format files from the clipboard.  It includes a patch to read.fwf which will not read from the clipboard or from remote file.  See read.fwf for documentation of how to specify the widths. }
\value{the contents of the clipboard.
}

\author{ William Revelle}

\examples{
#my.data <- read.clipboad()
#my.data <- read.clipboard.csv()
#my.data <- read.clipboad(header=FALSE)
#my.matrix <- read.clipboard.lower()
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
