\name{Utility}
\alias{fileScan}
\alias{fileCreate}
\alias{filesList}
\alias{filesInfo}

\alias{Utility}

\title{Useful utility functions for file/directory exploration and manipulation.} 
\description{
Wrappers for  dirname, file.choose, readLines. file.create, file.path to be called directly for listing directories, creating files, showing the files in a directory, and listing the content of files in a directory. \code{\link{fileCreate}} gives the functionality of \code{\link{file.choose}}(new=TRUE). \code{\link{filesList}} combines file.choose, dirname, and list.files to show the files in a directory, \code{\link{fileScan}} extends this and then returns the first few lines of each readable file
}
\usage{
fileScan(f = NULL, nlines = 3, max = NULL, from = 1, filter = NULL)
filesList(f=NULL)
filesInfo(f=NULL,max=NULL)
fileCreate(newName="new.file")
}
\arguments{
  \item{f}{File path to use as base path (will use file.choose() if missing.  If f is a 
  directory, will list the files in that directory, if f is a file, will find the 
  directory for that file and then list all of those files.)
}
  \item{nlines}{How many lines to display}
  \item{max}{maximum number of files to display}
  \item{from}{First file (number) to display}
  \item{filter}{Just display files with "filter" in the name}
  \item{newName}{The name of the file to be created.}
}
\details{
Just a collection of simple wrappers to powerful core R functions.  Allows the user more direct control of what directory to list, to create a file, or to display the content of files.  The functions called  include \code{\link{file.choose}}, \code{\link{file.path}}, \code{\link{file.info}},\code{\link{file.create}}, \code{\link{dirname}}, and \code{\link{dir.exists}}.  All of these are very powerful functions, but not easy to call interactively.
}

\author{William Revelle}
\note{Work arounds for core-R functions for interactive file manipulation }


\seealso{\code{\link{read.file}} to read in data from a file or \code{\link{read.clipboard}} from the clipboard.  \code{\link{dfOrder}} to sort data.frames.

}
\examples{
#all of these require interactive input and thus are not given as examples
#fileCreate("my.new.file.txt")
#filesList()   #choose the items in the directory where a file is displayed
#fileScan() #show the content of the files in a director 
#or, if you have a file in mind
# f <- file.choose()  #go find it
#filesList(f)
#fileScan(f)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
