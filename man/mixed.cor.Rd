\name{mixedCor}
\alias{mixedCor}
\alias{mixed.cor}
\title{Find correlations for mixtures of continuous, polytomous, and dichotomous variables}
\description{For data sets with continuous, polytomous and dichotmous variables, the absolute Pearson correlation is downward biased from the underlying latent correlation.  mixedCor finds Pearson correlations for the continous variables, \code{\link{polychoric}}s for the polytomous items, \code{\link{tetrachoric}}s for the dichotomous items, and the \code{\link{polyserial}} or \code{\link{biserial}} correlations for the various mixed variables. Results include the complete correlation matrix, as well as the separate correlation matrices and difficulties for the polychoric and tetrachoric correlations.
}
\usage{
mixedCor(data=NULL,c=NULL,p=NULL,d=NULL,smooth=TRUE,correct=.5,global=TRUE,ncat=8,
             use="pairwise",method="pearson",weight=NULL)
             
mixed.cor(x = NULL, p = NULL, d=NULL,smooth=TRUE, correct=.5,global=TRUE, 
        ncat=8,use="pairwise",method="pearson",weight=NULL)  #deprecated
}
\arguments{
  \item{data}{The data set to be analyzed (either a matrix or dataframe)}
  \item{c}{The names (or locations) of the continuous variables) (may be missing)}
  \item{x}{A set of continuous variables (may be missing) or, if p and d are missing, the variables to be analyzed.}
  \item{p}{A set of polytomous items (may be missing)}
  \item{d}{A set of dichotomous items (may be missing)}
  \item{smooth}{If TRUE, then smooth the correlation matix if it is non-positive definite}
  \item{correct}{When finding tetrachoric correlations, what value should be used to correct for continuity?}
  \item{global}{For polychorics, should the global values of the tau parameters be used, or should the pairwise values be used.  Set to local if errors are occurring.}
  \item{ncat}{The number of categories beyond which a variable is considered "continuous".}
  \item{use}{The various options to the \code{\link{cor}} function include "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". The default here is "pairwise"}
  \item{method}{The correlation method to use for the continuous variables. "pearson" (default), "kendall", or "spearman"}
  \item{weight}{If specified, this is a vector of weights (one per participant) to differentially weight participants. The NULL case is equivalent of weights of 1 for all cases. }
}
\details{
This function is particularly useful as part of the Synthetic Apeture Personality Assessment (SAPA) (\url{https://www.sapa-project.org/}) data sets where continuous variables (age, SAT V, SAT Q, etc) and mixed with polytomous personality items taken from the International Personality Item Pool (IPIP) and the dichotomous experimental IQ items that have been developed as part of SAPA (see, e.g., Revelle, Wilt and Rosenthal, 2010 or Revelle, Dworak and Condon, 2020.).  

This is a very computationally intensive function which can be speeded up considerably by using multiple cores and using the parallel package. (See the note for timing comparisons.) This adjusts the  number of cores to use when doing polychoric or tetrachoric. The greatest step in speed is going from 1 core to 2.  This is about a 50\% savings.  Going to 4 cores seems to have about at 66\% savings, and 8 a 75\% savings.  The number of parallel processes defaults to 2 but can be modified by using the \code{\link{options}} command:  options("mc.cores"=4) will set the number of cores to 4.

Item response analyses using \code{\link{irt.fa}} may be done separately on the polytomous and dichotomous items  in order to develop internally consistent scales. These scale may, in turn, be correlated with each other using the complete correlation matrix found by mixed.cor and using the \code{\link{score.items}} function.

This function is not quite as flexible as the hetcor function in John Fox's polychor package. 

Note that the variables may be organized by type of data:  continuous, polytomous, and dichotomous. This is done by simply specifying c, p, and d. This is advantageous in the case of some continuous variables having a limited number of categories because of subsetting.  

\code{\link{mixedCor}} is essentially a wrapper for \code{\link{cor}}, \code{\link{polychoric}}, \code{\link{tetrachoric}}, \code{\link{polydi}} and \code{\link{polyserial}}. It first identifies the types of variables, organizes them by type (continuous, polytomous, dichotomous), calls the appropriate correlation function, and then binds the resulting matrices together.   
}
\value{
 \item{rho}{The complete matrix}
 \item{rx}{The Pearson correlation matrix for the continuous items}
 \item{poly}{the polychoric correlation (poly$rho) and the item difficulties (poly$tau)}
 \item{tetra}{the tetrachoric correlation (tetra$rho) and the item difficulties (tetra$tau)}
}

\author{William Revelle}
\note{mixedCor was designed for the SAPA project (\url{https://www.sapa-project.org/}) with large data sets with a mixture of continuous, dichotomous, and polytomous data.  For smaller data sets, it is sometimes the case that the global estimate of the tau parameter will lead to unstable solutions.  This may be corrected by setting the global parameter = FALSE.  

\code{\link{mixedCor}} was added in April, 2017 to be slightly more user friendly. \code{\link{mixed.cor}} was deprecated in February, 2018.

When finding correlations between dummy coded SAPA data (e.g., of occupations), the real correlations are all slightly less than zero because of the ipsatized nature of the data.  This leads to a non-positive definite correlation matrix because the matrix is no longer of full rank. Smoothing will correct this, even though this might not be desired.  Turn off smoothing in this case.

 Note that the variables no longer need to be organized by type of data: first continuous, then polytomous, then dichotomous.  However, this automatic detection will lead to problems if the variables such as age are limited to less than 8 categories but those category values differ from the polytomous items.  The fall back is to specify x, p, and d. 
 
 If the number of alternatives in the polychoric data differ and there are some dicthotomous data, it is advisable to set correct=0.
 
Timing:  For large data sets, \code{\link{mixedCor}} takes a while.  Progress messages \code{\link{progressBar}}  report what is happening but do not actually report the rate of progress. ( The steps are initializing, Pearson r, polychorics, tetrachoric, polydi).  It is recommended to use the multicore option although the benefit between 2, 4 and 8 cores seems fairly small: For large data sets (e.g., 255K subjects, 950 variables), with 4 cores running  in parallel  (options("mc.cores=4") on a MacBook Pro with 2.8 Ghz Intel Core I7, it took 2,873/2,887 seconds elapsed time, 8,152/7,718 secs of user time, and 1,762/1,489 of system time (with and without smoothing). This is noticeabably better than the 4,842 elapsed time (7,313 user, 1,459 system) for 2 cores but not much worse than running 8 virtual cores, with an elapsed time of 2,629, user time of 13,460, and system time of 2,679. On a Macbook Pro with 2 physical cores and a 3.3 GHz Intel Cor I7, running 4 multicores took 4,423 seconds elapsed time, 12,781 seconds of user time, and 2,605 system time.  Running with 2 multicores, took slightly longer: 6,193 seconds elapsed time, 10,099 user time and 2,413 system time.  
 
 }

\references{
W.Revelle, J.Wilt, and A.Rosenthal. Personality and cognition: The personality-cognition link. In A.Gruszka, G. Matthews, and B. Szymura, editors,  Handbook of Individual Differences in Cognition: Attention, Memory and Executive
  Control, chapter 2, pages 27-49. Springer, 2010.

W Revelle, D. M. Condon,  J. Wilt,  J.A. French,  A. Brown, and L G. Elleman(2016) Web and phone based data collection using planned missing designs in Nigel G. Fielding and Raymond M. Lee and Grant Blank (eds) SAGE Handbook of Online Research Methods, Sage Publications, Inc.

W. Revelle, E.M. Dworak and D.M. Condon (2020) Exploring the persome: The power of the item in understanding personality structure. Personality and Individual Differences, \code{/10.1016/j.paid.2020.109905}

}

\seealso{ \code{\link{polychoric}}, \code{\link{tetrachoric}}, \code{\link{scoreItems}},  \code{\link{scoreOverlap}} \code{\link{scoreIrt}}
}
\examples{
data(bfi) 
r <- mixedCor(data=bfi[,c(1:5,26,28)])
r
#this is the same as
r <- mixedCor(data=bfi,p=1:5,c=28,d=26)
r #note how the variable order reflects the original order in the data
#compare to raw Pearson
#note that the biserials and polychorics are not attenuated
rp <- cor(bfi[c(1:5,26,28)],use="pairwise")
lowerMat(rp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{models}