\name{cor.plot}
\Rdversion{1.1}
\alias{cor.plot}
\alias{corPlot}
\alias{corPlotUpperLowerCi}
\alias{cor.plot.upperLowerCi}

\title{Create an image plot for a correlation or factor matrix}

\description{Correlation matrices may be shown graphically by using the image function to emphasize structure.  This is a particularly useful tool for showing the structure of  correlation matrices with a clear structure.  Partially meant for the pedagogical value of the graphic for teaching or discussing factor analysis and other multivariate techniques.
}
\usage{
corPlot(r,numbers=TRUE,colors=TRUE,n=51,main=NULL,zlim=c(-1,1),
  show.legend=TRUE, labels=NULL,n.legend=10,keep.par=TRUE,select=NULL, pval=NULL, 
  digits=2, trailing=TRUE, cuts=c(.001,.01),scale=TRUE,cex,MAR,upper=TRUE,diag=TRUE, 
  symmetric=TRUE,stars=FALSE, adjust="holm",xaxis=1, xlas=0, ylas=2, gr=NULL, alpha=.75, 
  min.length=NULL,...)

corPlotUpperLowerCi(R,numbers=TRUE,cuts=c(.001,.01,.05),select=NULL,
      main="Upper and lower confidence intervals of correlations",adjust=FALSE,...)
 
 cor.plot(r,numbers=TRUE,colors=TRUE,n=51,main=NULL,zlim=c(-1,1),
 show.legend=TRUE, labels=NULL,n.legend=10,keep.par=TRUE,select=NULL, pval=NULL, 
  digits=2, trailing=TRUE, cuts=c(.001,.01),scale=TRUE,cex,MAR,upper=TRUE,diag=TRUE, 
  symmetric=TRUE,stars=FALSE, adjust="holm",xaxis=1, xlas=0, ylas=2, gr=NULL, alpha=.75, 
  min.length=NULL,...)     #deprecated
           
cor.plot.upperLowerCi(R,numbers=TRUE,cuts=c(.001,.01,.05),select=NULL,
      main="Upper and lower confidence intervals of correlations",adjust=FALSE,...)
      #deprecated
}
\arguments{
  \item{r}{A correlation matrix or the output of  \code{\link{fa}}, \code{\link{principal}} or \code{\link{omega}}, or a raw data matrix. }
  \item{R}{The object returned from  \code{\link{cor.ci}}  }
  \item{numbers}{Display the numeric value of the correlations. (As of September, 2019)  Defaults to TRUE.}
  \item{colors}{Defaults to TRUE and colors use colors from the colorRampPalette from red through white to blue, but colors=FALSE will use a grey scale}
  \item{n}{The number of levels of shading to use.  Defaults to 51}
  \item{main}{A title. Defaults to "correlation plot"}
  \item{zlim}{The range of values to color -- defaults to -1 to 1.  If specified as NULL, then defaults to min and max observed correlation.}
  \item{show.legend}{A legend (key) to the colors is shown on the right hand side}
  \item{labels}{if NULL, use column and row names, otherwise use labels}
  \item{n.legend}{How many categories should be labelled in the legend?}
  \item{keep.par}{restore the graphic parameters when exiting}
  \item{pval}{scale the numbers by their pvals, categorizing them based upon the values of cuts}
  \item{digits}{Round off to digits.  Defaults to 2.}
  \item{trailing}{Show trailing zeros.  }
  \item{cuts}{Scale the numbers by the categories defined by pval < cuts}
  \item{scale}{Should the size of the numbers be scaled by the significance level?}
  \item{select}{Select the subset of variables to plot}
  \item{cex}{Character size.  Should be reduced a bit for large numbers of variables.}
  \item{MAR}{Allows for adjustment of the margins if using really long labels or big fonts}
  \item{upper}{Should the upper off diagonal matrix be drawn, or left blank?}
  \item{diag}{Should we show the diagonal?}
  \item{symmetric}{By default, if given a non-symmetric matrix, we find the correlations using pair.wise complete and then show them.  If wanting to display a non-symmetric matrix, then specify that symmetric is FALSE}
  \item{stars}{For those people who like to show the 'significance' of correlations by using magic astricks, set stars=TRUE}
  \item{adjust}{If showing significance, should we adjust for multiple tests?  The default is to show zero order probabilities below the diagonal and adjust these using the 'holm' correction above the diagonal. Use adjust = "none" if no adjustment is desired.
  adjust is also used in corPlotUpperLowerCI to show the nominal alpha confidence intervals (adjust =FALSE) or the Bonferonni adjusted confidence intervals (adjust=TRUE).}
  \item{xlas}{Orientation of the x axis labels (1 = horizontal, 0, parallel to axis, 2 perpendicular to axis)}
  \item{ylas}{Orientation of the y axis labels (1 = horizontal, 0, parallel to axis, 2 perpendicular to axis)}
  \item{xaxis}{By default, draw this below the figure.  If xaxis=3, then it wil be drawn above the figure}
  \item{gr}{A color gradient: e.g.,  gr <- colorRampPalette(c("#B52127", "white", "#2171B5"))  will produce slightly more pleasing (to some) colors. See next to last example.}
  \item{alpha}{The degree of transparency (0 = completely, 1= not). Default value of .75 makes somewhat moreor pleasing plots when using numbers.}
  \item{min.length}{If not NULL, then the maximum number of characters to use in 
    row/column labels}
  \item{...}{Other parameters for axis (e.g., cex.axis to change the font size, srt to   
  rotate the numbers in the plot)}
}
\details{When summarizing the correlations of large data bases or when teaching about factor analysis or cluster analysis, it is useful to graphically display the structure of correlation matrices.  This is a simple graphical display using the image function. 

The difference between mat.plot with a regular image plot is that the primary diagonal goes from the top left to the lower right. 
zlim defines how to treat the range of possible values. -1 to 1 and the color choice is more reasonable.  Setting it as c(0,1) will lead to negative correlations  treated as zero.  This is advantageous when showing general factor structures, because it makes the 0 white.  

There is an interesting case when plotting correlations corrected for attenuation.  Some of these might exceed 1.  In this case, either set zlim = NULL (to use the observed maximum and minimum values) or all values above 1 will be given a slightly darker shade than 1, but do not differ.

The default shows a legend for the color coding on the right hand side of the figure.

Inspired, in part, by a paper by S. Dray (2008)  on the number of components problem. 

Modified following suggestions by David Condon and Josh Wilt to use a more meaningful color choice ranging from dark red (-1) through white (0) to dark blue (1). Further modified to allow for color choices using the gr option (suggested by Lorien Elleman). Further modified to include the numerical value of the correlation.  (Inspired by the corrplot package).  These values may be scaled according the the probability values found in \code{\link{cor.ci}} or \code{\link{corr.test}}.

Unless specified, the font size is dynamically scaled to have a cex =  10/max(nrow(r),ncol(r).  This can produce fairly small fonts for large problems.  
The font size of the labels may be adjusted using cex.axis which defaults to one.  

By default \code{\link{cor.ci}} calls corPlotUpperLowerCi and scales the correlations based upon "significance" values.  The correlations plotted are the upper and lower confidence boundaries.  To show the correlations themselves, call corPlot directly.



If using the output of \code{\link{corr.test}}, the upper off diagonal will be scaled by the corrected probability, the lower off diagonal the scaling is the uncorrected probabilities.


If given raw data or  correlation matrix, \code{\link{corPlotUpperLowerCi}} will automatically call \code{\link{corr.test}} or \code{\link{cor.ci}}.  


If using the output of \code{\link{corr.test}} or \code{\link{cor.ci}} as input to \code{\link{corPlotUpperLowerCi}}, the upper off diagonal will be the upper bounds and the lower off diagonal the lower bounds of the confidence intervals.  If adjust=TRUE, these will use the Holm or Bonferroni adjusted values (depending upon corr.test). 

To compare the elements of two correlation matrices, \code{\link{corPlot}} the results from \code{\link{lowerUpper}}.  

To do multiple \code{\link{corPlot}} on the same plot, specify that show.legend=FALSE and keep.par=FALSE.  See the last examples.  

Care should be taken when selecting rows and columns from a non-symmetric matrix (e.g., the corrected correlations from \code{\link{scoreItems}} or \code{\link{scoreOverlap}}).  

To show a factor loading matrix (or any non-symmetric matrix), set symmetric=FALSE.  Otherwise the input will be treated as raw data and correlations will be found.

To sort the matrix use \code{\link{mat.sort}} first.  To find correlations other than Pearson, plot the output from e.g., \code{\link{mixed.cor}}. 
}
\author{William Revelle
}

\references{
Dray, Stephane (2008) On the number of principal components: A test of dimensionality based on measurements of similarity between matrices. Computational Statistics & Data Analysis. 52, 4, 2228-2237.
}

\seealso{ \code{\link{fa}}, \code{\link{mat.sort}}, \code{\link{cor.ci}}, \code{\link{corr.test}} \code{\link{lowerUpper}}. 
}
\examples{
corPlot(Thurstone,main="9 cognitive variables from Thurstone") 
#just blue implies positive manifold
#select just some variables to plot
corPlot(Thurstone, zlim=c(0,1),main="9 cognitive variables from Thurstone",select=c(1:3,7:9))
#now show a non-symmetric plot
corPlot(Thurstone[4:9,1:3], zlim=c(0,1),main="9 cognitive variables
 from Thurstone",numbers=TRUE,symmetric=FALSE)

#Two ways of including stars to show significance
#From the raw data
corPlot(sat.act,numbers=TRUE,stars=TRUE)
#from a correlation matrix with pvals
cp <- corr.test(sat.act)  #find the correlations and pvals
r<- cp$r
p <- cp$p
corPlot(r,numbers=TRUE,diag=FALSE,stars=TRUE, pval = p,main="Correlation plot
 with Holm corrected 'significance'")

#now red means less than .5
corPlot(mat.sort(Thurstone),TRUE,zlim=c(0,1), 
       main="9 cognitive variables from Thurstone (sorted by factor loading) ")
simp <- sim.circ(24)
corPlot(cor(simp),main="24 variables in a circumplex")

#scale by raw and adjusted probabilities
rs <- corr.test(sat.act[1:200,] ) #find the probabilities of the correlations
corPlot(r=rs$r,numbers=TRUE,pval=rs$p,main="Correlations scaled by probability values") 
 #Show the upper and lower confidence intervals
corPlotUpperLowerCi(R=rs,numbers=TRUE) 

#now do this again, but with lighter colors
gr <- colorRampPalette(c("#B52127", "white", "#2171B5"))
corPlot(r=rs$r,numbers=TRUE,pval=rs$p,main="Correlations scaled by probability values",gr=gr) 

corPlotUpperLowerCi(R=rs,numbers=TRUE,gr=gr) 




#do multiple plots 
#Also show the xaxis option
op <- par(mfrow=c(2,2))
corPlot(psychTools::ability,show.legend=FALSE,keep.par=FALSE,upper=FALSE)
f4 <- fa(psychTools::ability,4)
corPlot(f4,show.legend=FALSE,keep.par=FALSE,numbers=TRUE,xlas=3)
om <- omega(psychTools::ability,4)
corPlot(om,show.legend=FALSE,keep.par=FALSE,numbers=TRUE,xaxis=3)
par(op)


corPlotUpperLowerCi(rs,adjust=TRUE,main="Holm adjusted confidence intervals",gr=gr)

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line

