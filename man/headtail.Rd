\name{headTail}
\alias{headtail}
\alias{headTail}
\alias{topBottom}
\alias{quickView}
\title{Combine calls to head and tail}
\description{A quick way to show the first and last n lines of a data.frame, matrix, or a text object.  Just a pretty call to \code{\link{head}} and \code{\link{tail}} or \code{\link{View}}
}
\usage{
headTail(x, top=4,bottom=4,from=1,to=NULL, digits=2, hlength = 4, tlength =4, 
ellipsis=TRUE) 
headtail(x,hlength=4,tlength=4,digits=2,ellipsis=TRUE,from=1,to=NULL)
topBottom(x, top=4,bottom=4,from=1,to=NULL, digits=2, hlength = 4, tlength = 4) 
quickView(x,top=8,bottom=8,from=1,to=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix or data frame or free text}
  \item{top}{The number of lines at the beginning to show}
   \item{bottom}{The number of lines at the end to show}
  \item{digits}{Round off the data to digits}
  \item{ellipsis}{Separate the head and tail with dots (ellipsis)}
  \item{from}{The first column to show (defaults to 1)}
  \item{to}{The last column to show (defaults to the number of columns)}
  \item{hlength}{The number of lines at the beginning to show (an alias for top)}
   \item{tlength}{The number of lines at the end to show (an alias for bottom)}
  
}

\value{The first top and last bottom lines of a matrix or data frame with an ellipsis in between. If the input is neither a matrix nor data frame, the output will be the first top and last bottom lines. 

For each line, just columns starting at from and going to to will be displayed.  Bt default, from = 1 and to = the last column.

topBottom is just a call to headTail with ellipsis = FALSE and returning a matrix output.

quickView is a call to \code{\link{View}} which opens a viewing window which is scrollable (if needed because the number of lines listed is more than a screen's worth).  View (and therefore quickView) is slower than \code{\link{headTail}} or  \code{\link{topBottom}}.  
  
}


\seealso{ \code{\link{head}} and \code{\link{tail}}}

\examples{

headTail(psychTools::iqitems,4,8,to=6) #the first 4 and last 6 items from 1 to 6
topBottom(psychTools::ability,from =2, to = 6) #the first and last 4 items from 2 to 6
headTail(psychTools::bfi,top=4, bottom=4,from =6,to=10) #the first and last 4 from 6 to 10
#not shown
#quickView(ability,hlength=10,tlength=10)  #this loads a spreadsheet like table  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

