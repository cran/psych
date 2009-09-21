\name{00.psych}
\alias{psych}
\alias{psych-package}
\alias{00.psych-package}
\docType{package}
\title{A package for personality, psychometric, and psychological research
}
\description{Overview of the psych package.

The psych package has been developed at Northwestern University to include functions most useful for personality and psychological research.  Some of the functions (e.g., \code{\link{read.clipboard}}, \code{\link{describe}},  \code{\link{pairs.panels}}, \code{\link{error.bars}} ) are useful for basic data entry and descriptive analyses. Use help(package="psych") for a list of all functions.

Psychometric applications include routines (\code{\link{fa}} for principal axes (\code{\link{factor.pa}}), minimum residual (minres: \code{\link{factor.minres}}), and weighted least squares (\code{link{factor.wls}}  factor analysis  as well as functions to do Schmid Leiman transformations (\code{\link{schmid}}) to transform a hierarchical factor structure into a bifactor solution. Factor or components transformations to a target matrix include the standard Promax transformation (\code{\link{Promax}}), a transformation to a cluster target, or to any simple target matrix (\code{\link{target.rot}}) as well as the ability to call many of the GPArotation functions. Functions for determining the number of factors in a data matrix include Very Simple Structure (\code{\link{VSS}}) and Minimum Average Partial correlation (\code{\link{MAP}}). An alternative approach to factor analysis is Item Cluster Analysis (\code{\link{ICLUST}}).  Reliability coefficients alpha (\code{\link{score.items}}, \code{\link{score.multiple.choice}}), beta (\code{\link{ICLUST}}) and McDonald's omega (\code{\link{omega}} and  \code{\link{omega.graph}}) as well as Guttman's six estimates of internal consistency reliability (\code{\link{guttman}}) and the six measures of Intraclass correlation coefficients (\code{\link{ICC}}) discussed by Shrout and Fleiss are also available.

The \code{\link{score.items}}, and \code{\link{score.multiple.choice}} functions may be used to form single or multiple scales from sets of dichotomous, multilevel, or multiple choice items by specifying scoring keys.

Additional functions make for more convenient descriptions of item characteristics.  Functions under development include 1 and 2 parameter Item Response measures.

 A number of procedures have been developed as part of the Synthetic Aperture Personality Assessment (SAPA) project.  These routines facilitate forming and analyzing composite scales equivalent to using the raw data but doing so by adding within and between cluster/scale item correlations. These functions include extracting clusters from factor loading matrices (\code{\link{factor2cluster}}), synthetically forming clusters from correlation matrices (\code{\link{cluster.cor}}), and finding multiple  ((\code{\link{mat.regress}}) and partial ((\code{\link{partial.r}}) correlations from correlation matrices.
 
 Functions to generate simulated data with particular structures include \code{\link{sim.circ}} (for circumplex structures), \code{\link{sim.item}} (for general structures) and \code{\link{sim.congeneric}} (for a specific demonstration of congeneric measurement).  The functions \code{\link{sim.congeneric}}  and \code{\link{sim.hierarchical}} can be used to create data sets with particular structural properties. A more general form for all of these is \code{\link{sim.structural}} for generating general structural models.  These are discussed in more detail in the vignette (psych_for_sem).

 Functions to apply various standard statistical tests include \code{\link{p.rep}} and its variants for testing the probability of replication, \code{\link{r.con}} for the confidence intervals of a correlation, and \code{\link{r.test}} to test single, paired, or sets of correlations. 
 
 In order to study diurnal or circadian variations in mood, it is helpful to use circular statistics.  Functions to find the circular mean (\code{\link{circadian.mean}}), circular (phasic) correlations (\code{\link{circadian.cor}}) and the correlation between linear variables and circular variables (\code{\link{circadian.linear.cor}}) supplement a function to find the best fitting phase angle (\code{\link{cosinor}}) for measures taken with a fixed period (e.g., 24 hours).

 
 The most recent development version of the package is always available for download as a \emph{source} file from the repository at  \url{http://personality-project.org/r/src/contrib/}.  
 
}
\details{Two vignettes (overview.pdf) and psych_for_sem.pdf) are useful introductions to the package. They may be found as vignettes in R or may be downloaded from \url{http://personality-project.org/r/book/overview.pdf} and \url{http://personality-project.org/r/book/psych_for_sem.pdf}. 

The psych package was originally a combination of multiple source files maintained at the \url{http://personality-project.org/r} repository: ``useful.r", VSS.r., ICLUST.r, omega.r, etc.``useful.r" is a set of routines for easy data entry (\code{\link{read.clipboard}}),
simple descriptive statistics (\code{\link{describe}}), and splom plots combined with correlations (\code{\link{pairs.panels}}, adapted from the help files of pairs).  It is now a single package.

The \code{\link{VSS}} routines allow for testing the number of factors (\code{\link{VSS}}), showing plots (\code{\link{VSS.plot}}) of goodness of fit, and basic routines for estimating the number of factors/components to extract by using the \code{\link{MAP}}'s procedure, the examining the scree plot (\code{\link{VSS.scree}}) or comparing with the scree of an equivalent matrix of random numbers (\code{\link{VSS.parallel}}). 

In addition, there are routines for hierarchical factor analysis using Schmid Leiman tranformations (\code{\link{omega}}, \code{\link{omega.graph}}) as well as Item Cluster analysis (\code{\link{ICLUST}}, \code{\link{ICLUST.graph}}).

The more important functions in the package are for the analysis of multivariate data, with an emphasis upon those functions useful in scale construction of item composites. 

When given a set of items from a personality inventory, one goal is to combine these into higher level item composites. This leads to several questions:

1) What are the basic properties of the data?  \code{\link{describe}} reports basic summary statistics (mean, sd, median, mad, range,  minimum, maximum, skew, kurtosis, standard error) for vectors, columns of matrices, or data.frames. \code{\link{describe.by}} provides descriptive statistics, organized by one or more grouping variables.  \code{\link{pairs.panels}} shows scatter plot matrices (SPLOMs) as well as histograms and the Pearson correlation for scales or items. \code{\link{error.bars}} will plot variable means with associated confidence intervals. \code{\link{error.bars}} will plot confidence intervals for both the x and y coordinates.  \code{\link{corr.test}} will find the significance values for a matrix of correlations.

2) What is the most appropriate number of item composites to form? After finding  either standard Pearson correlations, or finding tetrachoric or polychoric correlations using a wrapper (\code{\link{poly.mat}}) for John Fox's hetcor function, the dimensionality of the correlation matrix may be examined. The number of factors/components problem is a standard question of factor analysis, cluster analysis, or principal components analysis. Unfortunately, there is no agreed upon answer. The Very Simple Structure (\code{\link{VSS}}) set of procedures has been proposed as on answer to the question of the optimal number of factors.  Other procedures (\code{\link{VSS.scree}},  \code{\link{VSS.parallel}},  \code{\link{fa.parallel}}, and \code{\link{MAP}})  also address this question.  

3) What are the best composites to form?  Although this may be answered using principal components  (\code{\link{principal}}), principal axis (\code{\link{factor.pa}}) or minimum residual (\code{\link{factor.minres}})  factor analysis (all part of the \code{\link{fa}} function) and to show the results graphically (\code{\link{fa.graph})}, it is sometimes more useful to address this question using cluster analytic techniques. (Some would argue that better yet is to use maximum likelihood factor analysis using \code{\link{factanal}} from the stats package.) Previous versions of \code{\link{ICLUST}} (e.g., Revelle, 1979)  have been shown to be particularly successful at forming maximally consistent and independent item composites.  Graphical output from \code{\link{ICLUST.graph}} uses the Graphviz dot language and allows one to write files suitable for Graphviz.  If Rgraphviz is available, these graphs can be done in R.

Graphical organizations of cluster and factor analysis output can be done using  \code{\link{cluster.plot}} which plots items by cluster/factor loadings and assigns items to that dimension with the highest loading.  

4) How well does a particular item composite reflect a single construct?  This is a question of reliability and general factor saturation.  Multiple solutions for this problem result in (Cronbach's) alpha (\code{\link{alpha}}, \code{\link{score.items}}), (Revelle's) Beta (\code{\link{ICLUST}}), and (McDonald's) \code{\link{omega}} (both omega hierarchical and omega total). Additional reliability estimates may be found in the \code{\link{guttman}} function.

5) For some applications, data matrices are synthetically combined from sampling different items for different people.  So called Synthetic Aperture Personality Assessement (SAPA) techniques allow the formation of large correlation or covariance matrices even though no one person has taken all of the items. To analyze such data sets, it is easy to form item composites based upon the covariance matrix of the items, rather than original data set.  These matrices may then be analyzed using a number of functions (e.g., \code{\link{cluster.cor}},   \code{\link{factor.pa}}, \code{\link{ICLUST}}, \code{\link{principal}}, \code{\link{mat.regress}}, and \code{\link{factor2cluster}}.

6) More typically, one has a raw data set to analyze. \code{\link{alpha}} will report several reliablity estimates as well as item-whole correlations for items forming a single scale, \code{\link{score.items}} will score data sets on multiple scales, reporting the scale scores, item-scale and scale-scale correlations, as well as coefficient alpha,  alpha-1 and G6+. Using a `keys' matrix (created by \code{\link{make.keys}} or by hand), scales can have overlapping or independent items. \code{\link{score.multiple.choice}} scores multiple choice items or converts multiple choice items to dichtomous (0/1) format for other functions. 

An additional set of functions generate simulated data to meet certain structural properties. \code{\link{sim.anova}} produces data simulating a 3 way analysis of variance (ANOVA) or linear model with or with out repeated measures. \code{\link{sim.item}} creates simple structure data,  \code{\link{sim.circ}} will produce circumplex structured data,  \code{\link{sim.dichot}} produces circumplex or simple structured data for dichotomous items.  These item structures are useful for understanding the effects of skew, differential item endorsement on factor and cluster analytic soutions.  \code{\link{sim.structural}} will produce correlation matrices and data matrices to match general structural models. (See the vignette).


When examining personality items, some people like to discuss them as representing items in a two dimensional space with a circumplex structure.  Tests of circumplex fit \code{\link{circ.tests}} have been developed.  When representing items in a circumplex, it is convenient to view them in \code{\link{polar}} coordinates. 

Additional functions for testing the difference between two independent or dependent correlation \code{\link{r.test}}, to find the \code{\link{phi}} or \code{\link{Yule}} coefficients from a two by table, or to find the confidence interval of a correlation coefficient.    

Ten data sets are included: \code{\link{bfi}} represents 25 personality items thought to represent five factors of personality, \code{\link{iqitems}} has 14 multiple choice iq items. \code{\link{sat.act}} has data on self reported test scores by age and gender. \code{\link{galton} }  Galton's data set of the heights of parents and their children. \code{\link{peas}} recreates the original Galton data set of the genetics of sweet peas.  \code{\link{heights}} and
\code{\link{cubits}} provide even more Galton data, \code{\link{vegetables}} provides the Guilford preference matrix of vegetables.  \code{\link{cities}} provides airline miles between 11 US cities (demo data for multidimensional scaling).

\tabular{ll}{
Package: \tab psych\cr
Type: \tab Package\cr
Version: \tab 1.0-80 \cr
Date: \tab 2009-9-30\cr
License: \tab GPL version 2 or newer\cr
}
Index:



\link{psych}      A package for personality, psychometric, and psychological research.\cr
 

Useful data entry and descriptive statistics\cr
\tabular{ll}{
\link{read.clipboard}       \tab      shortcut for reading from the clipboard\cr
\link{read.clipboard.csv}   \tab      shortcut for reading comma delimited files from clipboard \cr
\link{read.clipboard.lower}  \tab            shortcut for reading lower triangular matrices from the clipboard\cr
\link{read.clipboard.upper}   \tab           shortcut for reading upper triangular matrices from the clipboard\cr
\link{describe}        \tab           Basic descriptive statistics useful for psychometrics\cr
\link{describe.by}    \tab            Find summary statistics by groups\cr
\link{headtail}        \tab           combines the head and tail functions for showing data sets\cr

\link{pairs.panels}    \tab          SPLOM and correlations for a data matrix\cr
\link{corr.test}    \tab             Correlations, sample sizes, and p values  for a data matrix\cr
\link{cor.plot}    \tab             graphically show the size of correlations in a correlation matrix\cr
\link{multi.hist}      \tab         Histograms and densities of multiple variables arranged in matrix form\cr
\link{skew}           \tab          Calculate skew for a vector, each column of a matrix, or data.frame\cr
\link{kurtosi}        \tab          Calculate kurtosis for a vector, each column of a matrix or dataframe\cr
\link{geometric.mean}  \tab        Find the geometric mean of a vector or columns of a data.frame \cr
\link{harmonic.mean}  \tab         Find the harmonic mean of a vector or columns of a data.frame \cr
\link{error.bars}      \tab         Plot means and error bars \cr
\link{error.bars.by}   \tab           Plot means and error bars for separate groups\cr
\link{error.crosses}   \tab         Two way error bars \cr
\link{interp.median}  \tab           Find the interpolated median, quartiles, or general quantiles. \cr
\link{rescale}        \tab           Rescale data to specified mean and standard deviation \cr
\link{table2df}      \tab            Convert a two dimensional table of counts to a matrix or data frame \cr
}


Data reduction through cluster and factor analysis\cr

\tabular{ll}{
\link{fa}              \tab         Combined function for principal axis, minimum residual, or weighted least squares factor analysis\cr
\link{factor.pa}        \tab        Do a principal Axis factor analysis\cr
\link{factor.minres}    \tab        Do a minimum residual factor analysis\cr
\link{factor.wls}    \tab        Do a weighted least squares factor analysis\cr
\link{fa.graph}        \tab         Show the results of a factor analysis or principal components analysis graphically\cr
\link{fa.diagram}     \tab          Show the results of a factor analysis without using Rgraphviz \cr
\link{fa.sort}        \tab          Sort a factor or principal components output \cr
\link{principal}       \tab         Do an eigen value decomposition to find the principal components of a matrix\cr
\link{fa.parallel}     \tab         Scree test and Parallel analysis \cr
\link{factor.scores}     \tab           Estimate factor scores given a data matrix and factor loadings \cr
\link{guttman}         \tab          8 different measures of reliability (6 from Guttman (1945) \cr
\code{\link{ICLUST}}   \tab          Apply the ICLUST algorithm\cr
\link{ICLUST.graph}    \tab          Graph the output from ICLUST using the dot language\cr
\link{ICLUST.rgraph}   \tab          Graph the output from ICLUST using rgraphviz \cr
\link{poly.mat}        \tab          Find the polychoric correlations for items (uses J. Fox's hetcor) \cr
\link{omega}           \tab          Calculate the omega estimate of factor saturation (requires the GPArotation package)\cr
\link{omega.graph}      \tab         Draw a hierarchical or Schmid Leiman orthogonalized solution (uses Rgraphviz) \cr 
\link{schmid}            \tab        Apply the Schmid Leiman transformation to a correlation matrix\cr
\link{score.items}       \tab        Combine items into multiple scales and find alpha\cr
\link{score.multiple.choice}  \tab   Combine items into multiple scales and find alpha and basic scale statistics\cr
\link{smc}             \tab          Find the Squared Multiple Correlation (used for initial communality estimates)\cr
\link{VSS}            \tab           Apply the Very Simple Structure criterion to determine the appropriate number of factors.\cr
\link{VSS.parallel}    \tab          Do a parallel analysis to determine the number of factors for a random matrix\cr
\link{VSS.plot}       \tab           Plot VSS output\cr
\link{VSS.scree}      \tab           Show the scree plot of the factor/principal components\cr
\link{MAP}             \tab         Apply the Velicer Minimum Absolute Partial criterion for number of factors \cr
}

Functions for reliability analysis (some are listed above as well).
\tabular{ll}{
\link{alpha}      \tab         Find coefficient alpha and Guttman Lambda 6 for a scale (see also \link{score.items})\cr
\link{guttman}         \tab          8 different measures of reliability (6 from Guttman (1945) \cr
\link{omega}           \tab          Calculate the omega estimates of reliability (requires the GPArotation package)\cr
\link{ICC}          \tab             Intraclass correlation coefficients  \cr
\link{score.items}       \tab        Combine items into multiple scales and find alpha\cr
}

Procedures particularly useful for Synthetic Aperture Personality Assessment\cr


\tabular{ll}{
\link{alpha}      \tab         Find coefficient alpha and Guttman Lambda 6 for a scale (see also \link{score.items})\cr
\link{make.keys}        \tab         Create the keys file for score.items or cluster.cor            \cr
\link{correct.cor}      \tab        Correct a correlation matrix for unreliability\cr 
\link{count.pairwise}   \tab        Count the number of complete cases when doing pair wise correlations\cr
\link{cluster.cor}       \tab        find correlations of composite variables from larger matrix\cr
\link{cluster.loadings}  \tab        find correlations of items with  composite variables from a larger matrix\cr
\link{eigen.loadings}    \tab        Find the loadings when doing an eigen value decomposition\cr
\link{factor.pa}        \tab        Do a Principal Axis factor analysis and estimate factor scores\cr
\link{factor2cluster}    \tab       extract cluster definitions from factor loadings\cr
\link{factor.congruence}  \tab      Factor congruence coefficient\cr
\link{factor.fit}       \tab        How well does a factor model fit a correlation matrix\cr
\link{factor.model}    \tab          Reproduce a correlation matrix based upon the factor model\cr
\link{factor.residuals} \tab        Fit = data - model\cr
\link{factor.rotate}   \tab         ``hand rotate" factors\cr
\link{guttman}         \tab         8 different measures of reliability\cr
\link{mat.regress}     \tab         standardized multiple regression from raw or correlation matrix input\cr
}

Functions for generating simulated data sets \cr
\tabular{ll}{
\link{sim.anova}         \tab        Generate 3 independent variables and 1 or more dependent variables for demonstrating ANOVA and lm designs \cr
\link{sim.circ}         \tab        Generate a two dimensional circumplex item structure \cr
\link{sim.item}         \tab        Generate a two dimensional simple structrue with particular item characteristics \cr
\link{sim.congeneric}   \tab        Generate a one factor congeneric reliability structure \cr
\link{sim.structural}   \tab        Generate a multifactorial structural model  \cr
\link{sim.VSS}          \tab          Generate simulated data for the factor model\cr
\link{phi.demo}          \tab       Create artificial data matrices for teaching purposes\cr
\link{sim.hierarchical}	 \tab   Generate simulated correlation matrices with hierarchical or any structure\cr 
}


Graphical functions (require Rgraphviz) \cr
\tabular{ll}{
\link{structure.graph}  \tab           Draw a sem or regression graph \cr
\link{fa.graph}        \tab           Draw the factor structure from a factor or principal components analysis \cr
\link{omega.graph}      \tab           Draw the factor structure from an omega analysis (either with or without the Schmid Leiman transformation) \cr 
\link{ICLUST.graph}     \tab           Draw the tree diagram from ICLUST  \cr
}

Circular statistics (for circadian data analysis) \cr
\tabular{ll}{
\link{circadian.cor} \tab      Find the correlation with e.g., mood and time of day \cr
\link{circadian.linear.cor} \tab  Correlate a circular value with a linear value \cr
\link{circadian.mean} \tab      Find the circular mean of each column of a a data set     \cr
\link{cosinor} \tab            Find the best fitting phase angle for a circular data set \cr
}


 
Miscellaneous functions\cr
\cr
\tabular{ll}{
\link{comorbidity}   \tab             Convert base rate and comorbity to phi, Yule and tetrachoric\cr
\link{fisherz}      \tab             Apply the Fisher r to z transform\cr
\link{fisherz2r}    \tab             Apply the Fisher z to r transform\cr
\link{ICC}          \tab             Intraclass correlation coefficients  \cr
\link{cortest.mat}     \tab         Test for equality of two matrices (see also cortest.normal, cortest.jennrich ) \cr
\link{cortest.bartlett} \tab         Test whether a matrix is an identity matrix \cr
\link{paired.r}       \tab          Test for the difference of two paired or two independent correlations\cr
\link{r.con}          \tab            Confidence intervals for correlation coefficients \cr
\link{r.test}         \tab           Test of significance of r, differences between rs. \cr
\link{p.rep}          \tab           The probability of replication given a p, r, t, or F \cr
\link{phi}             \tab          Find the phi coefficient of correlation from a 2 x 2 table \cr
\link{phi.demo}        \tab          Demonstrate the problem of phi coefficients with varying cut points \cr
\link{phi2poly}         \tab         Given a phi coefficient, what is the polychoric correlation\cr
\link{phi2poly.matrix}  \tab         Given a phi coefficient, what is the polychoric correlation (works on matrices)\cr
\link{polar}           \tab          Convert 2 dimensional factor loadings to polar coordinates.\cr
\link{poly.mat}         \tab         Use John Fox's hetcor to create a matrix of correlations from a data.frame or matrix of integer values\cr
\link{polychor.matrix}   \tab        Use John Fox's polycor to create a matrix of polychoric  correlations from  a matrix of Yule correlations\cr
\link{scaling.fits}     \tab         Compares alternative scaling solutions and gives goodness of fits \cr
\link{thurstone}       \tab          Thurstone Case V scaling \cr
\link{tr}               \tab         Find the trace of a square matrix \cr
\link{wkappa}          \tab          weighted and unweighted versions of Cohen's kappa \cr
\link{Yule}            \tab           Find the Yule Q coefficient of correlation \cr
\link{Yule.inv}         \tab          What is the two by two table that produces a Yule Q with set marginals? \cr
\link{Yule2phi}         \tab          What is the phi coefficient corresponding to a Yule Q with set marginals? \cr
\link{Yule2phi.matrix}   \tab          Convert a matrix of Yule coefficients to a matrix of phi coefficients. \cr
\link{Yule2phi.matrix}   \tab          Convert a matrix of Yule coefficients to a matrix of polychoric  coefficients. \cr
}


Functions that are under development and not recommended for casual use \cr
\tabular{ll}{
\link{irt.item.diff.rasch}  \tab    IRT estimate of item difficulty with assumption that theta = 0\cr
\link{irt.person.rasch}	   \tab     Item Response Theory estimates of theta (ability) using a Rasch like model\cr\cr
}


Data sets included in the psych package \cr
\tabular{ll}{
\link{bfi}          \tab           represents 25 personality items thought to represent five factors of personality \cr
\link{bifactor}     \tab           8 different data sets with a bifactor structure  \cr
\link{cities}        \tab          The airline distances between 11 cities (used to demonstrate MDS) \cr
\link{epi.bfi}        \tab         13 personality scales \cr
\link{iqitems}        \tab         14 multiple choice iq items \cr
\link{msq}            \tab         75 mood items  \cr
\link{sat.act}        \tab          Self reported ACT and SAT Verbal and Quantitative scores by age and gender\cr

\link{galton}        \tab            Galton's data set of the heights of parents and their children \cr
\link{heights}       \tab            Galton's data set of the relationship between height and forearm (cubit) length \cr
\link{cubits}        \tab             Galton's data table of height and forearm length \cr
\link{peas}          \tab            Galton`s data set of the diameters of 700 parent and offspring sweet peas \cr
\link{vegetables}     \tab           Guilford`s preference matrix of vegetables (used for thurstone) \cr
}

A debugging function that may also be used as a demonstration of psych.
\tabular{ll}{
\link{test.psych}   \tab          Run a test of the major functions on 5 different data sets.  Primarily for development purposes. Although the output can be used as a demo of the various functions.
}

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
\link{structure.graph} \tab Rgraphviz \cr
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
