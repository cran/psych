\name{read.file}
\alias{read.clipboard}
\alias{read.clipboard.csv}
\alias{read.clipboard.tab}
\alias{read.clipboard.lower}
\alias{read.clipboard.upper}
\alias{read.clipboard.fwf}
\alias{read.file}
\alias{read.file.csv}
\alias{write.file}
\alias{write.file.csv}
\alias{read.https}
\title{Shortcuts for reading from the clipboard or a file}
\description{Input from a variety of sources may be read.  data.frames may be read from files with suffixes of .txt, .text, .dat, .data, .csv, .sav   (i.e., data from SPSS sav files may be easily read.) Data exported by JMP or EXCEL in the csv format are also able to be read. Files saved with writeRDS have suffixes of  .rds or  Rds,  and are read using readRDS.  Files associated with objects with suffixes .rda and .Rda are  loaded. The default values for read.spss are adjusted for more standard input from SPSS files.  
Input from the clipboard is easy but a bit obscure, particularly for Mac users. This is just an easier way to do so.  Data may be copied to the clipboard from Excel spreadsheets, csv files, or fixed width formatted files and then into a data.frame. Data may also be read from lower  (or upper) triangular matrices and filled out to square matrices.
}
\usage{
read.file(f=NULL,header=TRUE,use.value.labels=FALSE,to.data.frame=TRUE,...)
 #for .txt, .text, .csv, .sav, .rds, or .rda files

read.clipboard(header = TRUE, ...)   #assumes headers and tab or space delimited
read.clipboard.csv(header=TRUE,sep=',',...)   #assumes headers and comma delimited
read.clipboard.tab(header=TRUE,sep='\t',...)   #assumes headers and tab delimited
 #read in a matrix given the lower off diagonal
 read.clipboard.lower(diag=TRUE,names=FALSE,...) 
 read.clipboard.upper(diag=TRUE,names=FALSE,...)

#read in data using a fixed format width (see read.fwf for instructions)
read.clipboard.fwf(header=FALSE,widths=rep(1,10),...)  
read.https(filename,header=TRUE)

read.file.csv(f=NULL,header=TRUE,...)
write.file(x,f=NULL,row.names=FALSE,...)
write.file.csv(x,f=NULL,row.names=FALSE,...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{Does the first row have variable labels (generally assumed to be TRUE). }
  \item{sep}{What is the designated separater between data fields?}
  \item{diag}{for upper or lower triangular matrices, is the diagonal specified or not}
  \item{names}{for read.clipboard.lower or upper, are colnames in the  the first column}
  \item{widths}{how wide are the columns in fixed width input.  The default is to read 10 columns of size 1. }
  \item{filename}{Name or address of remote https file to read.}
  \item{\dots}{ Other parameters to pass to read }
  \item{f}{A file name to read.  If omitted, \code{\link{file.choose}} is called to dynamically get the file name.}
  \item{x}{The data frame or matrix to write to f}
  \item{row.names}{Should the output file include the rownames?} \item{to.data.frame}{Should the spss input be converted to a data frame?}
 \item{use.value.labels}{Should the SPSS input values be converted to numeric?}
}
\details{A typical session of R might involve data stored in text files, generated online, etc. Although it is easy to just read from a file (particularly if using \code{\link{file.choose}}, copying from the file to the clipboard and then reading from the clipboard is also very convenient (and somewhat more intuitive to the naive user). This is particularly convenient when copying from a text book or article and just moving a section of text into R.)

The \code{\link{read.file}} function combines the \code{\link{file.choose}} and either \code{\link{read.table}} or \code{\link{read.spss}} (from foreign) or \code{\link{load}} or \code{\link{readRDS}} commands. By examining the file suffix, it chooses the appropriate way to read. For more complicated file structures, see the foreign package. For even more complicated file structures, see the rio package. 

Note that \code{\link{read.file}} assumes by default that the first row has column labels (header =TRUE).  If this is not true, then make sure to specify header = FALSE.  


Currently supported input formats are 
\tabular{ll}{
.sav       \tab      SPSS.sav files\cr
.csv       \tab      A comma separated file (e.g. from Excel or Qualtrics)\cr
.txt       \tab      A typical text file \cr
.text      \tab     A typical text file \cr
.data     \tab       A data file \cr
.dat      \tab       A data file \cr
.rds       \tab      A R data file \cr
.Rds      \tab       A R data file  (created by a write) \cr
.Rda      \tab       A R data structure (created using save) \cr
.rda      \tab       A R data structure (created using save) \cr
.RData     \tab       A R data structure (created using save) \cr
.rdata     \tab       A R data structure (created using save) \cr
.R        \tab       A R data structure created using dput \cr
.r        \tab       A R data structure created using dput \cr
}

The foreign function \code{\link{read.spss}} is used to read SPSS .sav files using the most common options.  Just as  \code{\link{read.spss}}  issues various warnings, so does  \code{\link{read.file}}.   In general, these can be ignored.  For more detailed information about using  \code{\link{read.spss}}, see the help pages in the foreign package.


The \code{\link{write.file}} function combines the \code{\link{file.choose}} and either \code{\link{write.table}} or \code{\link{saveRDS}}. By examining the file suffix, it chooses the appropriate way to write. For more complicated file structures, see the foreign package, or the save function in R Base.  

Currently supported output formats are 
\tabular{ll}{
.csv       \tab      A comma separated file (e.g. for reading into Excel)\cr
.txt       \tab      A typical text file \cr
.text      \tab     A typical text file \cr
.rds       \tab      A R data file \cr
.Rds      \tab       A R data file  (created by a write) \cr
.Rda      \tab       A R data structure (created using save) \cr
.rda      \tab       A R data structure (created using save) \cr
.R        \tab       A R data structure created using dput \cr
.r        \tab       A R data structure created using dput \cr
}



\code{\link{read.clipboard}} was based upon a suggestion by Ken Knoblauch to the R-help listserve.

If the input file that was copied into the clipboard was an Excel file with blanks for missing data, then read.clipboard.tab() will correctly replace the blanks with NAs.  Similarly for a csv file with blank entries, read.clipboard.csv will replace empty fields with NA.  

read.clipboard.lower and read.clipboard.upper are adapted from John Fox's read.moments function in the sem package.  They will read a lower (or upper) triangular matrix from the clipboard and return a full, symmetric matrix for use by factanal, \code{\link{factor.pa}} , \code{\link{ICLUST}} , etc.  If the diagonal is false, it will be replaced by 1.0s.  These two function were added to allow easy reading of examples from various texts and manuscripts with just triangular output. 

Many articles will report lower triangular matrices with variable labels in the first column.  read.clipboard.lower  will handle this case. Names must be in the first column if names=TRUE is specified.

Other articles will report upper triangular matrices with variable labels in the first row.  read.clipboard.upper will handle this.  Note that labels in the first column will not work for read.clipboard.upper.  The names, if present, must be in the first row.

read.clipboard.fwf will read fixed format files from the clipboard.  It includes a patch to read.fwf which will not read from the clipboard or from remote file.  See read.fwf for documentation of how to specify the widths. }
\value{the contents of the file to be read or of the clipboard.
}

\author{ William Revelle}

\examples{
#Because these are dynamic functions, they need to be run interactively and 
# can not be shown as examples.

#my.data <- read.file()   #search the directory for a file and then read it.
#my.data <- read.clipboad()
#my.data <- read.clipboard.csv()
#my.data <- read.clipboad(header=FALSE)
#my.matrix <- read.clipboard.lower()
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ IO }% __ONLY ONE__ keyword per line
