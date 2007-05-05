\name{00psych-package}
\alias{psych-package}
\alias{psych}
\docType{package}
\title{A package for personality, psychometric, and psychological research
}
\description{The psych package has been developed at Northwestern University to include functions most useful for personality and psychological research.  Some of the functions (e.g., \code{\link{describe}} and  \code{\link{pairs.panels}}  ) are useful for basic descriptive analyses.

Psychometric applications include routines for Very Simple Structure \code{\link{VSS}}, Item Cluster Analysis \code{\link{ICLUST}} as well as functions to do Schmid Leiman transformations \code{\link{schmid}} and to calculate reliability coefficients alpha (\code{\link{score.items}}), beta (\code{\link{ICLUST}}) and  \code{\link{omega}}.

Additional functions make for more convenient descriptions of item characteristics.  Functions under development include 1 and 2 parameter Item Response measures.

 A number of procedures have been developed as part of the Synthetic Aperture Personality Assessment (SAPA) project.  These routines facilitate forming and analyzing composite scales equivalent to using the raw data but doing so by adding within and between cluster/scale item correlations. These functions include extracting clusters from factor loading matrices (\code{\link{factor2cluster}}), synthetically forming clusters from correlation matrices (\code{\link{cluster.cor}}), and multiple correlation from matrices (\code{\link{mat.regress}}).
 
}
\details{ The psych package is a combination of multiple source files maintained at the \url{http://personality-project.org/r}: "useful.r", VSS.r., ICLUST.r, omega.r, etc.  The "useful.r" is a set of routines for easy data entry (\code{\link{read.clipboard}}),
simple descriptive statistics (\code{\link{describe}}), and splom plots combined with correlations \code{\link{pairs.panels}} (taken from the help files of cor).

The (\code{\link{VSS}} routines allow for testing the number of factors (\code{\link{VSS}}), showing plots (\code{\link{VSS.plot}}) of goodness of fit, and basic routines for psychometrics.

In addition, there are routines for hierarchical factor analysis using Schmid Leiman tranformations (\code{\link{omega}}) as well as Item Cluster analysis (\code{\link{ICLUST}}).

The more important functions in the package are for the analysis of multivariate data, with an emphasis upon those functions useful in scale construction of item composites.

When given a set of items from a personality inventory, one goal is to combine these into higher level item composites. This leads to several questions:

1) What is the most appropriate number of item composites to form?  The Very Simple Structure (\code{\link{VSS}}) set of procedures is geared to this question.  The number of factors/components problem is a standard question of factor analysis, cluster analysis, or principal components analysis. 

2) What are the best composites to form?  Although this may be answered using principal components  (\link{principal}) or factor analysis (\link{factor.pa}), it is sometimes more useful to address this question using cluster analytic techniques.  Previous versions of \code{\link{ICLUST}} (e.g., Revelle, 1979)  have been shown to be particularly successful at doing this.  Graphical output from \code{\link{ICLUST.graph}} uses the Graphviz dot language and requires Graphviz.  (As of April, 2007, R-graphviz did not work on the Intel-Mac.  Soon one should be able to use R-graphviz.)

3) How well does a particular item composite reflect a single construct?  This is a question of reliability and general factor saturation.  Multiple solutions for this problem result in (Cronbach's) alpha (\code{\link{score.items}}), (Revelle's) Beta (\code{\link{ICLUST}}), and (McDonald's) \code{\link{omega}}. Functions to estimate all three of these are included in psych.  

4) For some applications, data matrices are synthetically combined from sampling different items for different people.  So called Synthetic Aperture Personality Assessement (SAPA) techniques allow the formation of large correlation or covariance matrices even though no one person has taken all of the items. To analyze such data sets, it is easy to form item composites based upon the covariance matrix of the items, rather than original data set.  These matrices may then be analyzed using a number of functions (e.g., \code{\link{cluster.cor}},   \code{\link{factor.pa}}, \code{\link{ICLUST}}, \code{\link{principal} }, \code{\link{mat.regress}}, \code{\link{factor2cluster}}.

5) More typically, one has a raw data set to analyze. \code{\link{score.items}} will score data sets on multiple scales, reporting the scale scores, item-scale and scale-scale correlations, as well as coefficient alpha and alpha-1. Using a 'keys' matrix, scales can have overlapping or independent items.

\tabular{ll}{
Package: \tab psych\cr
Type: \tab Package\cr
Version: \tab 1.0-18\cr
Date: \tab 2007-5-8,\cr
License: \tab GPL version 2 or newer\cr
}
Index:

\link{psych}       A package for personality, psychometric, and psychological research.\cr
\cr
Useful data entry and descriptive statistics\cr

\link{describe}                Basic descriptive statistics useful for psychometrics\cr
\link{describe.by}             Find summary statistics by groups\cr
\link{read.clipboard}          shortcut for reading from the clipboard\cr
\link{read.clipboard.csv}      shortcut for reading comma delimited files from clipboard\cr
\link{pairs.panels}            SPLOM and correlations for a data matrix\cr
\link{panel.cor}                 part of pairs.panels\cr
\link{panel.cor.scale}           part of pairs.panels\cr
\link{panel.hist}                part of pairs.panels\cr
\link{multi.hist}             Histograms of multiple variables arranged in matrix form\cr
\link{skew}                   Calculate skew for a vector, each column of matrix, or data.frame\cr
\link{kurtosi}                Calculate kurtosis for a vector, each column of a matrix or dataframe\cr

\link{error.crosses}          Two way error bars \cr
\link{geometric.mean}        Find the geometric mean of a vector or columns of a data.frame \cr
\link{harmonic.mean}         Find the harmonic mean of a vector or columns of a data.frame \cr
\cr
Data reduction through cluster and factor analysis\cr
\cr
\code{\link{ICLUST}}           Apply the ICLUST algorithm\cr
\link{ICLUST.graph}            Graph the output from ICLUST using the dot language\cr

\link{omega}                   Calculate the omega estimate of factor saturation (requires the GPArotation  package\cr
\link{schmid}                  Apply the Schmid Leiman transformation to a correlation matrix\cr\cr
\link{score.items}    Combine items into multiple scales and find alpha\cr
\link{VSS}                     Apply the Very Simple Structure criterion to determine the appropriate number of factors.\cr
\link{VSS.parallel}            Do a parallel analysis to determine the number of factors for a random matrix\cr
\link{VSS.plot}                Plot VSS output\cr
\link{VSS.scree}               Show the scree plot of the factor/principal components\cr
\link{VSS.simulate}           Generate simulated data for the factor model\cr
\link{make.hierarchical}	  Generate simulated correlation matrices with hierarchical structure\cr\cr

Procedures particularly useful for Synthetic Aperture Personality Assessment\cr
\cr
\link{alpha.scale}             Find coefficient alpha for a scale (see also \link{score.items})\cr
\link{correct.cor}            Correct a correlation matrix for unreliability\cr 
\link{count.pairwise}        Count the number of complete cases when doing pair wise correlations\cr
\link{cluster.cor}             find correlations of composite variables from larger matrix\cr
\link{cluster.loadings}        find correlations of items with  composite variables from a larger matrix\cr
\link{eigen.loadings}          Find the loadings when doing an eigen value decomposition\cr
\link{factor.pa}              Do a principal Axis factor analysis\cr
\link{factor2cluster}         extract cluster definitions from factor loadings\cr
\link{factor.congruence}      Factor congruence coefficient\cr
\link{factor.fit}             How well does a factor model fit a correlation matrix\cr
\link{factor.model}            Reproduce a correlation matrix based upon the factor model\cr
\link{factor.residuals}       Fit = data - model\cr
\link{factor.rotate}          "hand rotate" factors\cr
\link{mat.regress}            multiple regression from matrix input\cr
\link{principal}               Do an eigen value decomposition to find the principal components of a matrix\cr
\cr
Miscellaneous functions\cr
\cr
\link{fisherz}                 Apply the Fisher r to z transform\cr
\link{paired.r}               Test for the difference of two paired correlations\cr
\link{phi2poly}                Given a phi coefficient, what is the polychoric\cr
\link{polychor.matrix}         Use John Fox's polycor to create a matrix of correlations\cr
\link{psycho.demo}             Create artificial data matrices for teaching purposes\cr
\cr
Functions that are under development and not recommended for casual use \cr
\link{irt.item.diff.rasch}    IRT estimate of item difficulty with assumption that theta = 0\cr
\link{irt.person.rasch}	      Item Response Theory estimates of theta (ability) using a Rasch like model\cr\cr

}
\note{Development versions of this package are maintained at the local repository \url{http://personality-project.org/r} along with further documentation.  Specify that you are downloading a source package.  


Some functions require other packages. Specifically, omega and schmid require the GPArotation  package, and phi2poly and polychor.matrix requires John Fox's polychor package. 
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
In addition, see An Introduction to Psychometric Theory with applications in R (Revelle, in preparation) at \url{http://personality-project.org/r/book/}}
}
\keyword{package}% __ONLY ONE__ keyword per line
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line


\examples{
#See the separate man pages 
}
