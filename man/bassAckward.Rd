\name{bassAckward}
\alias{bassAckward}
\alias{bassAckward.diagram}

\title{The Bass-Ackward factoring algorithm discussed by Goldberg
}
\description{
 Goldberg (2006) described a hierarchical factor structure organization from the ``top down".  The original idea was to do successive factor analyses from 1 to nf factors organized by factor score correlations from  one level to the next.  Waller (2007) discussed a simple way of doing this for components without finding the scores.  Using the factor correlations (from Gorsuch) to organize factors hierarchically results may be organized at many different levels. The algorithm may be applied to principal components (pca) or to true factor analysis.
}
\usage{
bassAckward(r, nfactors = 1, fm = "minres", rotate = "oblimin", scores = "tenBerge",
   adjust=TRUE, plot=TRUE,cut=.3, use = "pairwise", cor = "cor", weight = NULL,
   correct = 0.5,...)
bassAckward.diagram(x,digits=2,cut = .3,labels=NULL,marg=c(1.5,.5,1.0,.5),
   main="BassAckward",items=TRUE,sort=TRUE,lr=TRUE,curves=FALSE,organize=TRUE,
   values=FALSE,...) 
}
\arguments{
  \item{r}{A correlation matrix or a data matrix suitable for factoring}
  \item{nfactors}{Factors from 1 to nfactors will be extracted. If nfactors is a a vector, then just the number of factors specified in the vector will be extracted. (See examples).
}
  \item{fm}{Factor method.  The default is 'minres' factoring.  Although to be consistent with the original Goldberg article, we can also do principal components (fm ="pca").
}
  \item{rotate}{What type of rotation to apply.  The default for factors is oblimin.  Unlike the normal call to pca where the default is varimax, in bassAckward the default for pca is oblimin.}
 \item{scores}{What factor scoring algorithm should be used. The default is "tenBerge", other possibilities include "regression", or "bartlett"}
 \item{adjust}{If using any other scoring proceure that "tenBerge" should we adjust the correlations for the lack of factor score fit?}
  \item{plot}{By default draw a bassAckward diagram}
  \item{use}{How to treat missing data.  Use='pairwise" finds pairwise complete correlations.
}
  \item{cor}{What kind of correlation to find.  The default is Pearson.
}
  \item{weight}{Should cases be weighted?  Default, no.
}
  \item{correct}{If finding tetrachoric or polychoric correlations, what correction should be applied to empty cells (defaults to .5)
}
 
\item{x}{The object returned by bassAckward}
\item{digits}{Number of digits to display on each path}
\item{cut}{Values greater than the abs(cut) will be displayed in a path diagram.}
\item{labels}{Labels may be taken from the output of the bassAckward function or can be specified as a list.}
\item{marg}{Margins are set to be slightly bigger than normal to allow for a cleaner diagram}
\item{main}{The main title for the figure}
\item{items}{if TRUE, show the items associated with the factors}
\item{sort}{if TRUE, sort the items by factor loadings}
\item{lr}{Should the graphic be drawn left to right or top to bottom}
\item{curves}{Should we show the correlations between factors at the same level}
\item{organize}{Rename and sort the factors at two lowest levels for a more pleasing figure}
\item{values}{If TRUE, then show the percent variance accounted for by this factor.}
\item{...}{Other graphic parameters (e.g., cex)}
}
\details{
This is essentially a wrapper to the \code{\link{fa}} and \code{\link{pca}} combined with the \code{\link{faCor}} functions.  They are called repeatedly and then the weights from the resulting solutions are used to find the factor/component correlations. 

Although the default is do all factor solutions from 1 to the nfactors, this can be simplified by specifying just some of the factor solutions. Thus, for the 135 items of the spi, it is more reasonable to ask for 3,5, and 27 item solutions. 

The function \code{\link{bassAckward.diagram}} may be called using the \code{\link{diagram}} function or may be called directly.

The output from \code{\link{bassAckward.diagram}} is the sorted factor structure suitable for using \code{\link{fa.lookup}}.

Although not particularly pretty, it is possible to do Schmid-Leiman rotations at each level.  Specify the rotation as rotate="schmid".

}
\value{
\item{Call}{Echo the call}
\item{fm}{Echos the factor method used}
\item{fa}{A list of the factor loadings at each level}
\item{bass.ack}{A list of the factor correlations at each level}
\item{summary}{The factors at each level}
\item{sumnames}{Summary of the factor names}
\item{labels}{Factor labels including items for each level}
\item{r}{The correlation matrix analyzed}
\item{Phi}{The factor correlations at each level}
\item{fa}{The factor analysis loadings at each level, now includes Phi values}
\item{fa.vac}{The variance accounted for by each factor at each level}
}
\references{
Goldberg, L.R. (2006) Doing it all Bass-Ackwards: The development of hierarchical factor structures from the top down. Journal of Research in Personality, 40, 4, 347-358.

Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 


Waller, N. (2007), A general method for computing hierarchical component structures by Goldberg's Bass-Ackwards method, Journal of Research in Personality, 41, 4, 745-752,
 DOI: 10.1016/j.jrp.2006.08.005
}
\author{William Revelle}
\note{Goldberg calculated factor/component scores and then correlated these.  Waller suggests just looking at the unrotated components and then examining the correlations when rotating  different numbers of components.  I do not follow the Waller procedure, but rather find successive factors and then  find factor/component correlations following Gorsuch. 

It is important to note that the BassAckward solution is not a hierarchical solution in the standard meaning. The factors are not factors of factors as is found in a hierarchical model (e.g. \code{\link{sim.hierarchical}} or \code{\link{omega}}, but is merely way of organising solutions with a different number of factors. In each case, unlike \code{\link{omega}} the factors are of the original variables, not the lower level factors.  Thus, detailed statistics for any level of the hierarchy may be found by doing a factoring with that particular number of factors.  

To find basic statistics for the multiple factorings, examine the fa object.  For more detail just do the factoring at that level. 

To see the items associated with the lowest level factors, use \code{\link{fa.lookup}}. For the items associated with other levels, use \code{\link{fa.lookup}} specifying the level.  (See examples.)
}

\seealso{ \code{\link{fa}}, \code{\link{pca}},  \code{\link{omega}} and \code{\link{iclust}}  for alternative hierarchical solutions. 
\code{link{fa.lookup}} to show the items in the lowest level of the solution.

}
\examples{
bassAckward(Thurstone,4,main="Thurstone data set")
f.labels <- list(level1=cs(Approach,Avoid),level2=cs(PosAffect,NegAffect,Constraint), 
level3 = cs(Extraversion,Agreeableness,Neuroticism,Conscientiousness,Openness))

ba <- bassAckward(bfi[1:25],c(2,3,5),labels=f.labels,
        main="bfi data set from psychTools", values=TRUE)
print(ba,short=FALSE)
#show the items associated with the 5 level solution
fa.lookup(ba,dictionary=bfi.dictionary)
#now show the items associated with the 3 level solution
fa.lookup(ba$fa[[2]],dictionary=bfi.dictionary)
# compare the 3 factor solution to what get by extracting 3 factors directly
f3 <- fa(bfi[1:25],3)
f3$loadings - ba$fa[[2]]$loadings   # these are the same
 #do pca instead of factors  just summarize, don't plot
summary(bassAckward(bfi[1:25],c(1,3,5,7),fm="pca",
       main="Components",plot=FALSE))
##not run, but useful example
\donttest{  
f.labels <- list(level1 = cs(Neg,Pos,Constrain),level2 = cs(Extra,Neuro,Cons,Open,Agree),
level3 = cs(attnseeking,sociability,impulsivity,
   charisma,sensationseek,emotexpr,humor,anxiety,
   emotstab,irritability,wellbeing,industry,order,author,honesty,perfect,easygoing,
   selfcontrol,conservatism,creativity,introspect,art,
   intellect,conform,adaptability,compassion,trust))
if(require("psychTools") ) {
sp5 <- bassAckward(psychTools::spi[11:145], c(3,5,27),labels=f.labels,
           main="spi data set from psychTools")}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{ models}