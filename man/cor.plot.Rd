\name{cor.plot}
\Rdversion{1.1}
\alias{cor.plot}
\alias{corPlot}
\alias{cor.plot.upperLowerCi}
\title{Create an image plot for a correlation or factor matrix}

\description{Correlation matrices may be shown graphically by using the image function to emphasize structure.  This is a particularly useful tool for showing the structure of  correlation matrices with a clear structure.  Partially meant for the pedagogical value of the graphic for teaching or discussing factor analysis and other multivariate techniques.
}
\usage{
corPlot(r,numbers=FALSE,colors=TRUE,n=51,main=NULL,zlim=c(-1,1),
  show.legend=TRUE, labels=NULL,n.legend=10,keep.par=TRUE,select=NULL,
   pval=NULL,cuts=c(.001,.01),cex,MAR,upper=TRUE,diag=TRUE,...)
cor.plot(r,numbers=FALSE,colors=TRUE,n=51,main=NULL,zlim=c(-1,1),
  show.legend=TRUE, labels=NULL,n.legend=10,keep.par=TRUE,select=NULL,
   pval=NULL,cuts=c(.001,.01),cex,MAR,upper=TRUE,diag=TRUE,...)
   
          
cor.plot.upperLowerCi(R,numbers=TRUE,cuts=c(.001,.01,.05),select=NULL,
      main="Upper and lower confidence intervals of correlations",...)
}
\arguments{
  \item{r}{A correlation matrix or the output of  \code{\link{fa}}, \code{\link{principal}} or \code{\link{omega}}.  }
  \item{R}{The object returned from  \code{\link{cor.ci}}  }
  \item{numbers}{Display the numeric value of the correlations.  Defaults to FALSE.}
  \item{colors}{Defaults to TRUE and colors use colors from the colorRampPalette from red through white to blue, but colors=FALSE will use a grey scale}
  \item{n}{The number of levels of shading to use.  Defaults to 51}
  \item{main}{A title. Defaults to "correlation plot"}
  \item{zlim}{The range of values to color -- defaults to -1 to 1}
  \item{show.legend}{A legend (key) to the colors is shown on the right hand side}
  \item{labels}{if NULL, use column and row names, otherwise use labels}
  \item{n.legend}{How many categories should be labelled in the legend?}
  \item{keep.par}{restore the graphic parameters when exiting}
  \item{pval}{scale the numbers by their pvals, categorizing them based upon the values of cuts}
  \item{cuts}{Scale the numbers by the categories defined by pval < cuts}
  \item{select}{Select the subset of variables to plot}
  \item{cex}{Character size.  Should be reduced a bit for large numbers of variables.}
  \item{MAR}{Allows for adjustment of the margins if using really long labels or big fonts}
  \item{upper}{Should the upper off diagonal matrix be drawn, or left blank?}
  \item{diag}{Should we show the diagonal?}
  \item{...}{Other parameters for axis (e.g., cex.axis to change the font size, srt to rotate the numbers in the plot)}
}
\details{When summarizing the correlations of large data bases or when teaching about factor analysis or cluster analysis, it is useful to graphically display the structure of correlation matrices.  This is a simple graphical display using the image function. 

The difference between mat.plot with a regular image plot is that the primary diagonal goes from the top left to the lower right. 
zlim defines how to treat the range of possible values. -1 to 1 and the color choice is more reasonable.  Setting it as c(0,1) will lead to negative correlations  treated as zero.  This is advantageous when showing general factor structures, because it makes the 0 white.  

The default shows a legend for the color coding on the right hand side of the figure.

Inspired, in part, by a paper by S. Dray (2008)  on the number of components problem. 

Modified following suggestions by David Condon and Josh Wilt to use a more meaningful color choice ranging from dark red (-1) through white (0) to dark blue (1).  Further modified to include the numerical value of the correlation.  (Inspired by the corrplot package).  These values may be scaled according the the probability values found in \code{\link{cor.ci}} or \code{\link{corr.test}}.

Unless specified, the font size is dynamically scaled to have a cex =  10/max(nrow(r),ncol(r).  This can produce fairly small fonts for large problems.  
The font size of the labels may be adjusted using cex.axis which defaults to one.  

By default \code{\link{cor.ci}} calls cor.plot.upperLowerCi and scales the correlations based upon "significance" values.  The correlations plotted are the upper and lower confidence boundaries.  To show the correlations themselves, call corPlot directly.

If using the output of \code{\link{corr.test}}, the upper off diagonal will be scaled by the corrected probability, the lower off diagonal the scaling is the uncorrected probabilities.

If using the output of \code{\link{corr.test}} or \code{\link{cor.ci}} as input to \code{\link{cor.plot.upperLowerCi}}, the upper off diagonal will be the upper bounds and the lower off diagonal the lower bounds of the confidence intervals. 


To do multiple \code{\link{corPlot}} on the same plot, specify that show.legend=FALSE and keep.par=FALSE.  

}
\author{William Revelle
}

\references{
Dray, Stephane (2008) On the number of principal components: A test of dimensionality based on measurements of similarity between matrices. Computational Statistics \& Data Analysis. 52, 4, 2228-2237.
}

\seealso{ \code{\link{fa}}, \code{\link{mat.sort}}, \code{\link{cor.ci}}, \code{\link{corr.test}}. 
}
\examples{
cor.plot(Thurstone,main="9 cognitive variables from Thurstone") 
#just blue implies positive manifold
#select just some variables to plot
cor.plot(Thurstone, zlim=c(0,1),main="9 cognitive variables from Thurstone",select=1:4)

#now red means less than .5
cor.plot(mat.sort(Thurstone),TRUE,zlim=c(0,1), 
       main="9 cognitive variables from Thurstone (sorted by factor loading) ")
simp <- sim.circ(24)
cor.plot(cor(simp),main="24 variables in a circumplex")

#scale by raw and adjusted probabilities
rs <- corr.test(sat.act[1:200,] ) #find the probabilities of the correlations
cor.plot(r=rs$r,numbers=TRUE,pval=rs$p,main="Correlations scaled by probability values") 
 #Show the upper and lower confidence intervals
cor.plot.upperLowerCi(R=rs,numbers=TRUE) 

#do multiple plots
op <- par(mfrow=c(2,2))
corPlot(ability,show.legend=FALSE,keep.par=FALSE,upper=FALSE)
f4 <- fa(ability,4)
corPlot(f4,show.legend=FALSE,keep.par=FALSE,numbers=TRUE)
om <- omega(ability,4)
corPlot(om,show.legend=FALSE,keep.par=FALSE,numbers=TRUE)
par(op)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line

