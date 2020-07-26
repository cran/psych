\name{pairwiseCount}
\alias{pairwiseCount}
\alias{count.pairwise}
\alias{pairwiseDescribe}
\alias{pairwiseReport}
\alias{pairwiseImpute}
\alias{pairwisePlot}
\title{Count number of pairwise cases for a data set with missing (NA) data and impute values. }
\description{
  When doing cor(x, use= "pairwise"), it is nice to know the number of cases for each pairwise correlation.  This is particularly useful when doing SAPA type analyses. More importantly, when there are some missing pairs, it is useful to supply imputed values so that further analyses may be done.  This is useful if using the Massively Missing Completely at Random (MMCAR) designs used by the SAPA project. 
}
\usage{
pairwiseCount(x, y = NULL,diagonal=TRUE)
pairwiseDescribe(x,y,diagonal=FALSE,...) 
pairwiseImpute(keys,R,fix=FALSE)
pairwiseReport(x,y=NULL,cut=0,diagonal=FALSE,...) 
pairwisePlot(x,y=NULL,upper=TRUE,diagonal=TRUE,labels=TRUE,show.legend=TRUE,n.legend=10,
colors=FALSE,gr=NULL,minlength=6,xlas=1,ylas=2,
main="Relative Frequencies",count=TRUE,...)
count.pairwise(x, y = NULL,diagonal=TRUE) #deprecated
}

\arguments{
  \item{x}{ An input matrix, typically a data matrix ready to be correlated. }
  \item{y}{ An optional second input matrix }
  \item{diagonal}{if TRUE, then report the diagonal, else fill the diagonals with NA}
  \item{...}{Other parameters to pass to describe}
  \item{keys}{A keys.list specifying which items belong to which scale.}
  \item{R}{A correlation matrix to be described or imputed}
  \item{cut}{Report the item pairs and numbers with cell sizes less than cut}
  \item{fix}{If TRUE, then replace all NA correlations with the mean correlation for that  
  within or between scale}
   \item{upper}{Should the upper off diagonal matrix be drawn, or left blank?}
  \item{labels}{if NULL, use column and row names, otherwise use labels}
  \item{show.legend}{A legend (key) to the colors is shown on the right hand side}
   \item{n.legend}{How many categories should be labelled in the legend?}
    \item{colors}{Defaults to FALSE and  will use a grey scale. colors=TRUE use colors \
    from the colorRampPalette from red through white to blue}
     \item{minlength}{If not NULL, then the maximum number of characters to use in 
    row/column labels}
     \item{xlas}{Orientation of the x axis labels (1 = horizontal, 0, parallel to axis, 2 perpendicular to axis)}
  \item{ylas}{Orientation of the y axis labels (1 = horizontal, 0, parallel to axis, 2 perpendicular to axis)}
   \item{main}{A title. Defaults to "Relative Frequencies"}
    \item{gr}{A color gradient: e.g.,  gr <- colorRampPalette(c("#B52127", "white", "#2171B5"))  will produce slightly more pleasing (to some) colors. See next to last example of \code{\link{corPlot}}.  }
 \item{count}{Should we count the number of pairwise observations using pairwiseCount, or just plot the counts for a matrix?}
    
}

\details{When using Massively Missing Completely at Random (MMCAR) designs used by the SAPA project, it is important to count the number of pairwise observations (\code{\link{pairwiseCount}}).  If there are pairs with 1 or fewer observations, these will produce NA values for correlations making subsequent factor analyses \code{\link{fa}} or reliability analsyes \code{\link{omega}} or \code{\link{scoreOverlap}} impossible.

In order to identify item pairs with counts less than a certain value \code{\link{pairwiseReport}} reports the names of those pairs with fewer than 'cut' observations.  By default, it just reports the number of offending items. With short=FALSE, the print will give the items with n.obs < cut. Even more detail is available in the returned objects.

  To remedy the problem of missing correlations, we impute the missing correlations using \code{\link{pairwiseImpute}}.
The technique takes advantage of the scale based structure of SAPA items.  Items within a scale (e.g. Letter Number Series similar to the \code{\link[psychTools]{ability}} items) are  imputed to correlate with items from another scale (e.g., Matrix Reasoning) at the average of these two between scale inter-item mean correlations.  The average correlations within and between scales are reported by \code{\link{pairwiseImpute}} and if the fix paremeter is specified, the imputed correlation matrix is returned. 

Alternative methods of imputing these correlations are not yet implemented.
 }

\value{
\item{result}{ = matrix of counts of pairwise observations (if pairwiseCount)}
\item{av.r}{The average correlation value of the observed correlations
  within/between scales} 
\item{count}{The numer of observed correlations within/between each scale}
\item{percent}{The percentage of complete data by scale}
\item{imputed}{The original correlation matrix with NA values replaced by the mean  
correlation for items within/between the appropriate scale.}
}

\author{

Maintainer: William Revelle \email{revelle@northwestern.edu}}

\examples{

x <- matrix(rnorm(900),ncol=6)
y <- matrix(rnorm(450),ncol=3)
x[x < 0] <- NA
y[y > 1] <- NA

pairwiseCount(x)
pairwiseCount(y)
pairwiseCount(x,y)
pairwiseCount(x,diagonal=FALSE)
pairwiseDescribe(x,quant=c(.1,.25,.5,.75,.9))

#examine the structure of the ability data set
keys <- list(ICAR16=colnames(psychTools::ability),reasoning =  
  cs(reason.4,reason.16,reason.17,reason.19),
  letters=cs(letter.7, letter.33,letter.34,letter.58, letter.7), 
  matrix=cs(matrix.45,matrix.46,matrix.47,matrix.55), 
  rotate=cs(rotate.3,rotate.4,rotate.6,rotate.8))
 pairwiseImpute(keys,psychTools::ability)

    
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
