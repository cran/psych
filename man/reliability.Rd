\name{reliability}
\alias{reliability}
\alias{plot.reliability}


\title{Reports 7 different estimates of scale reliabity including alpha, omega, split half}

\description{

Revelle and Condon, (2019) reviewed the problem of reliability in a tutorial meant to useful to the theoretician as well as the practitioner.  Although there are a number of functions in psych for estimating reliability of single scales, (e.g. \code{\link{alpha}} and \code{\link{omega}}), for split half reliability \code{\link{splitHalf}} or for finding test-retest reliability \code{\link{testRetest}} or multilevel reliability \code{\link{mlr}}, the \code{\link{reliability}} function  combines several of these functions to report these recommended measures for multiple scales. 

To quote from Revelle and Condon (2019) ``Reliability is a fundamental problem for measurement in all of science for `(a)ll measurement is befuddled by error' (p 294 McNemar, 1946). Perhaps because psychological measures are more befuddled than those of the other natural sciences, psychologists have long studied the problem of reliability.   

``Issues of reliability are fundamental to understanding how correlations between observed variables are (attenuated) underestimates of the relationships between the underlying constructs, how observed estimates of a person's score are biased estimates of their latent score, and how to estimate the confidence intervals around any particular measurement. Understanding the many ways to estimate reliability as well as the ways to use these estimates allows one to better assess individuals and to evaluate selection and prediction techniques. This is not just a problem for measurement specialists but for all who want to make theoretical inferences from observed data.

`` It is no longer  acceptable to report one  coefficient that is only correct if all items are exactly equally good measures of a construct.  Researchers are encouraged to report at least two coefficients (e.g., omega_h and omega_t) and then discuss why each is appropriate for the inference that is being made.  They are discouraged from reporting just alpha unless they can justify the assumptions implicit in using it (i.e., tau equivalence and unidimensionality)."  Here we make it easy to do so.
 
Although the \code{\link{alpha}} and \code{\link{omega}} functions will find reliability estimates for a single scale, and \code{\link{scoreItems}} and \code{\link{scoreOverlap}} will find alpha for multiple scales, it sometimes is convenient to call \code{\link{omega}} and \code{\link{splitHalf}} for multiple scales.  \code{\link{reliability}} takes a keys list (suitable for \code{\link{scoreItems}} ) and then finds hierarchical and total omega as well as split half reliabilities for each separate scale. 

\code{\link{plot.reliability}} takes the output of  \code{\link{reliability}} and displays it as a dot chart showing the values of both omegas as well as alpha and the distributions of split half reliabilities.
}
\usage{
reliability(keys=NULL, items, nfactors = 2, split = TRUE, raw=TRUE, plot=FALSE,hist=FALSE,
  n.sample=10000)
\method{plot}{reliability}(x,omega=TRUE,alpha=TRUE,split=TRUE,uni=TRUE,add=FALSE,
xlim=NULL, main=NULL,...)
}

\arguments{
  \item{keys}{A list of items to be scored (may be taken from a keys.list for \code{\link{scoreItems}}. This list may contain one or more keys. ) If keys are not specified, then all the items are used.

}
  \item{items}{The matrix or data.frame of items to be scored.  Can be substantially greater than the items included in keys. For just those items in each key are scored.
}
  \item{nfactors}{Omega is not well defined for two factors, but for small sets of items, two is the better choice. For larger number of items per scale, 3 is probably preferable.
}
  \item{split}{By default, find splitHalf reliabilities as well as the omega statistics. When plotting, split implies that raw was called in reliability.}
  
  \item{plot}{By default, suppress the omega plots for each scale.}
  
  
  
  \item{raw}{If TRUE, return a list of all the possible splits (up to n.samples).  Useful for graphic display. }
  
  \item{hist}{If TRUE then  split and raw are forced to TRUE and the histograms of the split half values are drawn.  (Otherwise, just return the values for later plotting)}
  
  \item{n.sample}{Normally defaults to 10,000.  This means that for up to 16 item tests, all possible splits are found.  choose(n,n/2)/2 explodes above that, eg. for all splits of the epi E scale requires 1,352,078 splits or 23.4 seconds on a MacBook Pro with a 2.4GHZ 8 core Intel Core I9.  Can be done, but do you want to do so?  }
  
  \item{x}{The object returned from reliability}
  \item{omega}{Add in the values of omega_h and omega_t}
  \item{uni}{Show the unidimensionality value from \code{\link{unidim}}.  }
  \item{alpha}{Add the value of alpha}
  \item{add}{Allows us to merge this figure with other ones}
  \item{main}{Defaults to "Split half distributions + omega, alpha"}
  \item{xlim}{The xlim of the plot}
  \item{...}{Other graphical parameters}
      
  }
\details{\code{\link{reliability}} is basically just a wrapper for \code{\link{omegah}}, \code{\link{unidim}} and \code{\link{splitHalf}}. Revelle and Condon (2019) recommended reporting at least three reliability statistics for any scale, here we make it easy to do. 

If the hist option is set to true, histgrams and density plots of the split half values for each test are also shown. The output from \code{\link{reliability}} can be passed to \code{\link{error.dots}} to show the reliability statistics for multiple scales graphically.  It is however more useful to just call the \code{\link{plot.reliability}} function to show the basic information.

For detailed analysis of any one scale, it is recommended to do a complete \code{\link{omega}} analysis, perhaps combined with a \code{\link{splitHalf}} analysis. The \code{\link{reliability}} function is just meant for the case where the user has multiple scales (perhaps scored using \code{\link{scoreItems}}) and then wants to get more complete reliability information for all of the scales.  

Following a suggestion, the ability to not bother with keys and just do omega and split half and draw the results has been added.  Either specify that keys=NULL, or just specify the items to use. (See the first example.)

\code{\link{plot.reliability}} provides a dot chart summary of the distributions of the split half values, as well as the estimates of omega and alpha and unidimensionality.  It can also be called by just issuing a plot command.
}


\value{

 \item{omega_h}{Omega_h is the (model based) hierarchical estimate of the general factor saturation of a scale.}
  \item{alpha}{The conventional alpha statistic (which is not model based)}
  \item{omega.tot}{A model based estimate of the total reliability of a scale}
  \item{Uni}{An experimental estimate of unidimensionality (from unidim)}
  \item{r.fit}{How well does the average r of the correlations reproduce the matrix?}
  \item{f.fit}{How well does a single factor reproduce the correlation matrix}
  \item{max.split}{The greatest split half reliability of the scale.  Found by finding all possible splits (if this is < 10,000) or sampled from 10,000 possible splits.}
  \item{min.split}{The lowest split half reliability of the scale.  An estimate of beta (see \code{\link{iclust}}). }
  \item{mean.r}{The average correlation of the items in the scale}
  \item{med.r}{The median correlation of the items.  If this differs from the mean, that is a sign of poor scale.}
  \item{splits}{A list of the split half values for all possible splits.}

}
\references{

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{https://personality-project.org/r/book/} 

Revelle, W. and Condon, D.M. (2019) Reliability from alpha to omega: A tutorial.  Psychological Assessment, 31, 12, 1395-1411. https://doi.org/10.1037/pas0000754.  \url{https://psyarxiv.com/2y3w9/} Preprint available from PsyArxiv 

}
\author{
William Revelle
}
\note{For much more information on reliability, see the help pages for \code{\link{omega}}  as well as the Revelle and Condon (2019) tutorial or the Revelle (in prep) chapter on reliability. 
}


\seealso{See Also \code{\link{omega}} to find more complete output for the various omega analyses,\code{\link{splitHalf}} to show more detail on split half estimates, \code{\link{scoreItems}} to find scores on multiple scales using unit weights, \code{\link{testRetest}} to find test retest reliabilities, \code{\link{mlr}} to find multilevel reliabilities.  

\code{\link{predicted.validity}} will call \code{\link{reliability}} and \code{\link{item.validity}} to use the average r information to find the asymptotic validity of a set of scales for a set of criteria.

}
\examples{
reliability(psychTools::ability) #an example of finding reliability for all items
rel <- reliability(psychTools::ability.keys,psychTools::ability)  #use keys to select scales
R <- cor(psychTools::ability,use="pairwise") #find the correlations to test 
rel.R <- reliability(psychTools::ability.keys,R) #this should be the same as rel
plot(rel.R) #versus all and subsets
all.equal(rel$result.df,rel.R$result.df ) #should be TRUE
reliability(psychTools::bfi.keys,psychTools::bfi)  #reliability when items are keyed negative

\dontrun{  
#this takes a few seconds but shows nice graphic displays

spi.rel <- reliability(psychTools::spi.keys,psychTools::spi,hist=TRUE) #graph them
spi.rel #show them

#plot them using plot.reliability
plot(spi.rel)  #draw the density distrbutions

plot(spi.rel,split=FALSE) #don't draw the split half density distribution
plot(spi.rel,omega=FALSE) # don't add omega values to the diagram
#or do this without the densities

#plot the first three values in a dot chart
error.dots(spi.rel$result.df[,1],sort=FALSE, xlim=c(.3,1),head=16,tail=16,
          main = expression(paste(omega[h], ~~~~ alpha,~~~~ omega[t])))
  #plot the omega_h values
error.dots(spi.rel$result.df[,2],sort=FALSE,pch=2,xlim=c(.3,1),head=16,tail=16,
       main="",labels="",add=TRUE)#add the alpha values
error.dots(spi.rel$result.df[,3],sort=FALSE, xlim=c(.3,1),head=16,tail=16, 
       pch=3,labels="", main="",add=TRUE) #and the omega_t values
       
#or, show the smallest and greatest split half, as well as alpha
error.dots(spi.rel$result.df[,4],sort=FALSE, xlim=c(.3,1),head=16,tail=16,
          main = expression(paste(beta, ~~~~ alpha,~~~~  glb)))
error.dots(spi.rel$result.df[,5],sort=FALSE,pch=5,xlim=c(.3,1),head=16,tail=16,
       main="",labels="",add=TRUE)#add the GLB values
error.dots(spi.rel$result.df[,2],sort=FALSE,pch=2,xlim=c(.3,1),head=16,tail=16,
       main="",labels="",add=TRUE)#add the alpha values
       

} 
}
\keyword{ multivariate }
\keyword{ models }

