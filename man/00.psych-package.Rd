\name{psych}
\alias{psych}
\alias{psych-package}
\alias{00.psych-package}
\docType{package}
\title{A package for personality, psychometric, and psychological research
}
\description{Overview of the psych package.

The psych package has been developed at Northwestern University to include functions most useful for personality and psychological research.  Some of the functions (e.g., \code{\link{read.clipboard}}, \code{\link{describe}},  \code{\link{pairs.panels}}, \code{\link{error.bars}} ) are useful for basic data entry and descriptive analyses.

Psychometric applications include routines for Very Simple Structure (\code{\link{VSS}}), Minimum Average Partial correlation (\code{\link{MAP}}), Item Cluster Analysis (\code{\link{ICLUST}}) and principal axes factor analysis (\code{\link{factor.pa}}), as well as functions to do Schmid Leiman transformations (\code{\link{schmid}}) to transform a hierarchical factor structure into a bifactor solution, and to calculate reliability coefficients alpha (\code{\link{score.items}}, \code{\link{score.multiple.choice}}), beta (\code{\link{ICLUST}}) and McDonald's omega (\code{\link{omega}} and  \code{\link{omega.graph}}). Guttman's six estimates of internal consistency reliability (\code{\link{guttman}}) and the six measures of Intraclass correlation coefficients (\code{\link{ICC}}) discussed by Shrout and Fleiss are also available.

The \code{\link{score.items}}, and \code{\link{score.multiple.choice}} functions may be used to form single or multiple scales from sets of dichotomous, multilevel, or multiple choice items by specifying scoring keys.

Additional functions make for more convenient descriptions of item characteristics.  Functions under development include 1 and 2 parameter Item Response measures.

 A number of procedures have been developed as part of the Synthetic Aperture Personality Assessment (SAPA) project.  These routines facilitate forming and analyzing composite scales equivalent to using the raw data but doing so by adding within and between cluster/scale item correlations. These functions include extracting clusters from factor loading matrices (\code{\link{factor2cluster}}), synthetically forming clusters from correlation matrices (\code{\link{cluster.cor}}), and finding multiple  ((\code{\link{mat.regress}}) and partial ((\code{\link{partial.r}}) correlations from correlation matrices.
 
 Functions to generate simulated data with particular structures include \code{\link{sim.circ}} (for circumplex structures), \code{\link{sim.item}} (for general structures) and \code{\link{sim.congeneric}} (for a specific demonstration of congeneric measurement).  The functions \code{\link{sim.congeneric}}  and \code{\link{sim.hierarchical}} can be used to create data sets with particular structural properties. A more general form for all of these is \code{\link{sim.structural}} for generating general structural models.

 Functions to apply various standard statistical tests include \code{\link{p.rep}} and its variants for testing the probability of replication, \code{\link{r.con}} for the confidence intervals of a correlation, \code{\link{r.test}} to test single, paired, or sets of correlations. 
 
 In order to study diurnal or circadian variations in mood, it is helpful to use circular statistics.  Functions to find the circular mean (\code{\link{circadian.mean}}), circular (phasic) correlations (\code{\link{circadian.cor}}) and the correlation between linear variables and circular variables (\code{\link{circadian.linear.cor}}) supplement a function to find the best fitting phase angle (\code{\link{cosinor}}) for measures taken with a fixed period (e.g., 24 hours).

 
 The most recent development version of the package is always available for download as a \emph{source} file from the repository at  \url{http://personality-project.org/r/src/contrib/}.  
 
}
\details{ The psych package was originally a combination of multiple source files maintained at the \url{http://personality-project.org/r} repository: ``useful.r", VSS.r., ICLUST.r, omega.r, etc.``useful.r" is a set of routines for easy data entry (\code{\link{read.clipboard}}),
simple descriptive statistics (\code{\link{describe}}), and splom plots combined with correlations (\code{\link{pairs.panels}}, adapted from the help files of pairs).  It is now a single package.

The \code{\link{VSS}} routines allow for testing the number of factors (\code{\link{VSS}}), showing plots (\code{\link{VSS.plot}}) of goodness of fit, and basic routines for estimating the number of factors/components to extract by using the \code{\link{MAP}}'s procedure, the examining the scree plot (\code{\link{VSS.scree}}) or comparing with the scree of an equivalent matrix of random numbers (\code{\link{VSS.parallel}}). 

In addition, there are routines for hierarchical factor analysis using Schmid Leiman tranformations (\code{\link{omega}}, \code{\link{omega.graph}}) as well as Item Cluster analysis (\code{\link{ICLUST}}, \code{\link{ICLUST.graph}}).

The more important functions in the package are for the analysis of multivariate data, with an emphasis upon those functions useful in scale construction of item composites. 

When given a set of items from a personality inventory, one goal is to combine these into higher level item composites. This leads to several questions:

1) What are the basic properties of the data?  \code{\link{describe}} reports basic summary statistics (mean, sd, median, mad, range,  minimum, maximum, skew, kurtosis, standard error) for vectors, columns of matrices, or data.frames. \code{\link{describe.by}} provides descriptive statistics, organized by one or more grouping variables.  \code{\link{pairs.panels}} shows scatter plot matrices (SPLOMs) as well as histograms and the Pearson correlation for scales or items. \code{\link{error.bars}} will plot variable means with associated confidence intervals. \code{\link{error.bars}} will plot confidence intervals for both the x and y coordinates.

2) What is the most appropriate number of item composites to form? After finding  either standard Pearson correlations, or finding tetrachoric or polychoric correlations using a wrapper (\code{\link{poly.mat}}) for John Fox's hetcor function, the dimensionality of the correlation matrix may be examined. The number of factors/components problem is a standard question of factor analysis, cluster analysis, or principal components analysis. Unfortunately, there is no agreed upon answer. The Very Simple Structure (\code{\link{VSS}}) set of procedures has been proposed as on answer to the question of the optimal number of factors.  Other procedures (\code{\link{VSS.scree}},  \code{\link{VSS.parallel}},  \code{\link{fa.parallel}}, and \code{\link{MAP}})  also address this question.  

3) What are the best composites to form?  Although this may be answered using principal components  (\code{\link{principal})} or principal axis factor analysis (\code{\link{factor.pa})} and to show the results graphically (\code{\link{fa.graph})}, it is sometimes more useful to address this question using cluster analytic techniques. (Some would argue that better yet is to use maximum likelihood factor analysis using \code{\link{factanal}} from the stats package.) Previous versions of \code{\link{ICLUST}} (e.g., Revelle, 1979)  have been shown to be particularly successful at doing this.  Graphical output from \code{\link{ICLUST.graph}} uses the Graphviz dot language and allows one to write files suitable for Graphviz.  

Graphical organizations of cluster and factor analysis output can be done using  \code{\link{cluster.plot}} which plots items by cluster/factor loadings and assigns items to that dimension with the highest loading.  

4) How well does a particular item composite reflect a single construct?  This is a question of reliability and general factor saturation.  Multiple solutions for this problem result in (Cronbach's) alpha (\code{\link{score.items}}), (Revelle's) Beta (\code{\link{ICLUST}}), and (McDonald's) \code{\link{omega}} (both omega hierarchical and omega total). Additional reliability estimates may be found in the \code{\link{guttman}} function.

5) For some applications, data matrices are synthetically combined from sampling different items for different people.  So called Synthetic Aperture Personality Assessement (SAPA) techniques allow the formation of large correlation or covariance matrices even though no one person has taken all of the items. To analyze such data sets, it is easy to form item composites based upon the covariance matrix of the items, rather than original data set.  These matrices may then be analyzed using a number of functions (e.g., \code{\link{cluster.cor}},   \code{\link{factor.pa}}, \code{\link{ICLUST}}, \code{\link{principal} }, \code{\link{mat.regress}}, and \code{\link{factor2cluster}}.

6) More typically, one has a raw data set to analyze. \code{\link{score.items}} will score data sets on multiple scales, reporting the scale scores, item-scale and scale-scale correlations, as well as coefficient alpha and alpha-1. Using a `keys' matrix (created by \code{\link{make.keys}} or by hand), scales can have overlapping or independent items. \code{\link{score.multiple.choice}} scores multiple choice items or converts multiple choice items to dichtomous (0/1) format for other functions. 

An additional set of functions generate simulated data to meet certain structural properties.  \code{\link{item.sim}} creates simple structure data,  \code{\link{sim.circ}} will produce circumplex structured data,  \code{\link{item.dichot}} produces circumplex or simple structured data for dichotomous items.  These item structures are useful for understanding the effects of skew, differential item endorsement on factor and cluster analytic soutions.

When examining personality items, some people like to discuss them as representing items in a two dimensional space with a circumplex structure.  Tests of circumplex fit \code{\link{circ.tests}} have been developed.  When representing items in a circumplex, it is convenient to view them in \code{\link{polar}} coordinates. 

Additional functions for testing the difference between two independent or dependent correlation \code{\link{r.test}}, to find the \code{\link{phi}} or \code{\link{Yule}} coefficients from a two by table, or to find the confidence interval of a correlation coefficient.    

Ten data sets are included: \code{\link{bfi}} represents 25 personality items thought to represent five factors of personality, \code{\link{iqitems}} has 14 multiple choice iq items. \code{\link{sat.act}} has data on self reported test scores by age and gender. \code{\link{galton} }  Galton's data set of the heights of parents and their children. \code{\link{peas}} recreates the original Galton data set of the genetics of sweet peas.  \code{\link{heights}} and
\code{\link{cubits}} provide even more Galton data, \code{\link{vegetables}} provides the Guilford preference matrix of vegetables.  \code{\link{cities}} provides airline miles between 11 US cities (demo data for multidimensional scaling).

\tabular{ll}{
Package: \tab psych\cr
Type: \tab Package\cr
Version: \tab 1.0-61 \cr
Date: \tab 2009-1-03\cr
License: \tab GPL version 2 or newer\cr
}
Index:


\link{psych}        A package for personality, psychometric, and psychological research.\cr
  \cr
Useful data entry and descriptive statistics\cr

\link{describe}                 Basic descriptive statistics useful for psychometrics\cr
\link{describe.by}              Find summary statistics by groups\cr
\link{headtail}                 combines the head and tail functions for showing data sets\cr
\link{read.clipboard}           shortcut for reading from the clipboard\cr
\link{read.clipboard.csv}       shortcut for reading comma delimited files from clipboard \cr
\link{read.clipboard.lower}            shortcut for reading lower triangular matrices from the clipboard\cr
\link{read.clipboard.upper}            shortcut for reading upper triangular matrices from the clipboard\cr
\link{pairs.panels}            SPLOM and correlations for a data matrix\cr
\link{multi.hist}             Histograms and densities of multiple variables arranged in matrix form\cr
\link{skew}                   Calculate skew for a vector, each column of a matrix, or data.frame\cr
\link{kurtosi}                Calculate kurtosis for a vector, each column of a matrix or dataframe\cr
\link{geometric.mean}        Find the geometric mean of a vector or columns of a data.frame \cr
\link{harmonic.mean}         Find the harmonic mean of a vector or columns of a data.frame \cr
\link{error.bars}             Plot means and error bars \cr
\link{error.bars.by}            Plot means and error bars for separate groups\cr
\link{error.crosses}          Two way error bars \cr
\link{interp.median}           Find the interpolated median, quartiles, or general quantiles. \cr
\link{rescale}                 Rescale data to specified mean and standard deviation \cr
\link{table2df}                Convert a two dimensional table of counts to a matrix or data frame \cr
\
\cr
\cr
Data reduction through cluster and factor analysis\cr
\cr
\link{factor.pa}              Do a principal Axis factor analysis\cr
\link{fa.graph}                Show the results of a factor analysis or principal components analysis graphically\cr
\link{principal}               Do an eigen value decomposition to find the principal components of a matrix\cr
\link{fa.parallel}            Scree test and Parallel analysis \cr
\link{guttman}                 8 different measures of reliability (6 from Guttman (1945)
\code{\link{ICLUST}}           Apply the ICLUST algorithm\cr
\link{ICLUST.graph}            Graph the output from ICLUST using the dot language\cr
\link{ICLUST.rgraph}           Graph the output from ICLUST using rgraphviz \cr
\link{poly.mat}                Find the polychoric correlations for items (uses J. Fox's hetcor \cr
\link{omega}                   Calculate the omega estimate of factor saturation (requires the GPArotation package)\cr
\link{omega.graph}             Draw a hierarchical or Schmid Leiman orthogonalized solution \cr 
\link{schmid}                  Apply the Schmid Leiman transformation to a correlation matrix\cr
\link{score.items}             Combine items into multiple scales and find alpha\cr
\link{score.multiple.choice}   Combine items into multiple scales and find alpha and basic scale statistics\cr
\link{smc}                     Find the Squared Multiple Correlation (used for initial communality estimates)\cr
\link{VSS}                     Apply the Very Simple Structure criterion to determine the appropriate number of factors.\cr
\link{VSS.parallel}            Do a parallel analysis to determine the number of factors for a random matrix\cr
\link{VSS.plot}                Plot VSS output\cr
\link{VSS.scree}               Show the scree plot of the factor/principal components\cr

\link{make.hierarchical}	  Generate simulated correlation matrices with hierarchical structure\cr
\link{MAP}                    Apply the Velicer Minimum Absolute Partial criterion for number of factors \cr\cr

Procedures particularly useful for Synthetic Aperture Personality Assessment\cr
\cr
\link{alpha.scale}             Find coefficient alpha for a scale (see also \link{score.items})\cr
\link{make.keys}             Create the keys file for score.items or cluster.cor            \cr
\link{correct.cor}            Correct a correlation matrix for unreliability\cr 
\link{count.pairwise}        Count the number of complete cases when doing pair wise correlations\cr
\link{cluster.cor}             find correlations of composite variables from larger matrix\cr
\link{cluster.loadings}        find correlations of items with  composite variables from a larger matrix\cr
\link{eigen.loadings}          Find the loadings when doing an eigen value decomposition\cr
\link{factor.pa}              Do a Principal Axis factor analysis and estimate factor scores\cr
\link{factor2cluster}         extract cluster definitions from factor loadings\cr
\link{factor.congruence}      Factor congruence coefficient\cr
\link{factor.fit}             How well does a factor model fit a correlation matrix\cr
\link{factor.model}            Reproduce a correlation matrix based upon the factor model\cr
\link{factor.residuals}       Fit = data - model\cr
\link{factor.rotate}          ``hand rotate" factors\cr
\link{guttman}                8 different measures of reliability
\link{mat.regress}            standardized multiple regression from raw or correlation matrix input\cr
\link{principal}               Do an eigen value decomposition to find the principal components of a matrix\cr
\cr
Functions for generating simulated data sets \cr
\link{sim.circ}               Generate a two dimensional circumplex item structure \cr
\link{sim.item}               Generate a two dimensional simple structrue with particular item characteristics \cr
\link{sim.congeneric}         Generate a one factor congeneric reliability structure \cr
\link{sim.structural}         Generate a multifactorial structural model  \cr
\link{sim.VSS}                  Generate simulated data for the factor model\cr
\link{phi.demo}               Create artificial data matrices for teaching purposes\cr
\link{sim.hierarchical}	  Generate simulated correlation matrices with hierarchical or any structure\cr 
\cr

Circular statistics (for circadian data analysis) \cr
\link{circadian.cor}  \cr
\link{circadian.linear.cor}  \cr
\link{circadian.mean}  \cr
\link{cosinor} \cr
 \cr
Miscellaneous functions\cr
\cr
\link{comorbidity}              Convert base rate and comorbity to phi, Yule and tetrachoric
\link{fisherz}                 Apply the Fisher r to z transform\cr
\link{fisherz2r}                Apply the Fisher z to r transform\cr
\link{ICC}                     Intraclass correlation coefficients  \cr
\link{cortest.mat}            Test for equality of two matrices (see also cortest.normal, cortest.jennrich ) \cr
\link{cortest.bartlett}        Test whether a matrix is an identity matrix \cr
\link{paired.r}               Test for the difference of two paired or two independent correlations\cr
\link{r.con}                    Confidence intervals for correlation coefficients \cr
\link{r.test}                  Test of significance of r, differences between rs. \cr
\link{p.rep}                   The probability of replication given a p, r, t, or F \cr
\link{phi}                     Find the phi coefficient of correlation from a 2 x 2 table \cr
\link{phi.demo}                Demonstrate the problem of phi coefficients with varying cut points \cr
\link{phi2poly}                Given a phi coefficient, what is the polychoric correlation\cr
\link{phi2poly.matrix}         Given a phi coefficient, what is the polychoric correlation (works on matrices)\cr
\link{polar}                   Convert 2 dimensional factor loadings to polar coordinates.\cr
\link{poly.mat}                Use John Fox's hetcor to create a matrix of correlations from a data.frame or matrix of integer values\cr
\link{polychor.matrix}         Use John Fox's polycor to create a matrix of polychoric  correlations from  a matrix of Yule correlations\cr
\link{scaling.fits}            Compares alternative scaling solutions and gives goodness of fits \cr
\link{thurstone}               Thurstone Case V scaling \cr
\link{tr}                      Find the trace of a square matrix \cr
\link{wkappa}                  weighted and unweighted versions of Cohen's kappa \cr

\link{Yule}                     Find the Yule Q coefficient of correlation \cr
\link{Yule.inv}                 What is the two by two table that produces a Yule Q with set marginals? \cr
\link{Yule2phi}                 What is the phi coefficient corresponding to a Yule Q with set marginals? \cr
\link{Yule2phi.matrix}           Convert a matrix of Yule coefficients to a matrix of phi coefficients. \cr
\link{Yule2phi.matrix}           Convert a matrix of Yule coefficients to a matrix of polychoric  coefficients. \cr
\cr
Functions that are under development and not recommended for casual use \cr
\link{irt.item.diff.rasch}    IRT estimate of item difficulty with assumption that theta = 0\cr
\link{irt.person.rasch}	      Item Response Theory estimates of theta (ability) using a Rasch like model\cr\cr
\cr
Data sets included in the psych package \cr
\link{bfi}                   represents 25 personality items thought to represent five factors of personality \cr
\link{cities}                The airline distances between 11 cities (used to demonstrate MDS) \cr
\link{epi.bfi}               13 personality scales \cr
\link{iqitems}               14 multiple choice iq items \cr
\link{sat.act}                Self reported ACT and SAT Verbal and Quantitative scores by age and gender\cr
\link{galton}                  Galton's data set of the heights of parents and their children \cr
\link{heights}                 Galton's data set of the relationship between height and forearm (cubit) length \cr
\link{cubits}                   Galton's data table of height and forearm length \cr
\link{peas}                    Galton`s data set of the diameters of 700 parent and offspring sweet peas \cr
\link{vegetables}              Guilford`s preference matrix of vegetables (used for thurstone) \cr

\link{test.psych}             Run a test of the major functions on 5 different data sets.  Primarily for development purposes. Although the output can be used as a demo of the various functions.
}
\note{Development versions (source code) of this package are maintained at the  repository \url{http://personality-project.org/r} along with further documentation.   Specify that you are downloading a source package.  
\cr
Some functions require other packages. Specifically, omega and schmid require the GPArotation package, and poly.mat, phi2poly and polychor.matrix requires John Fox's polychor package. ICLUST.rgraph and fa.graph require Rgraphviz.  i.e.:
\cr

\tabular{ll}{
function  \tab     requires\cr
\link{omega}    \tab     GPArotation \cr
\link{schmid}    \tab     GPArotation\cr
\link{poly.mat}   \tab    polychor\cr
\link{phi2poly}   \tab    polychor\cr
\link{polychor.matrix}\tab polychor \cr
\link{ICLUST.rgraph} \tab   Rgraphviz \cr
\link{fa.graph}  \tab      Rgraphviz \cr
}

}


\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illiniois \cr
        \url{http://personality-project.org/revelle.html}\cr


Maintainer: William Revelle <revelle@northwestern.edu>

}
\references{A general guide to personality theory and research may be found at the personality-project \url{http://personality-project.org}. See also the short guide to R at \url{http://personality-project.org/r}.
In addition, see 

Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{http://personality-project.org/r/book/}    }

\keyword{package}% __ONLY ONE__ keyword per line
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line


\examples{
#See the separate man pages 
test.psych()
}
