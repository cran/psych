\name{00.psych}
\alias{psych}
\alias{psych-package}
\alias{00.psych-package}
\docType{package}
\title{A package for personality, psychometric, and psychological research}
\description{Overview of the psych package.

The psych package has been developed at Northwestern University to include functions most useful for personality and psychological research.  Some of the functions (e.g., \code{\link[psychTools]{read.file}}, \code{\link[psychTools]{read.clipboard}}, \code{\link{describe}},  \code{\link{pairs.panels}}, \code{\link{error.bars}} and \code{\link{error.dots}}) are useful for basic data entry and descriptive analyses. Use help(package="psych") or objects("package:psych") for a list of all functions  Three vignettes are included as part of the package.  The intro vignette tells how to install psych and overview vignette provides examples of using psych in many applications.  In addition, there are a growing set of tutorials available on the \url{https://personality-project.org/r/} webpages. 

A companion package \code{\link[psychTools]{psychTools}} includes larger data set examples and one more vignette. 

Psychometric applications include routines (\code{\link{fa}} for  maximum likelihood (fm="mle"), minimum residual (fm="minres"), minimum rank (fm=minrank)  principal axes (fm="pa") and weighted least squares (fm="wls")  factor analysis  as well as functions to do Schmid Leiman transformations (\code{\link{schmid}}) to transform a hierarchical factor structure into a bifactor solution. Principal Components Analysis (\code{\link{pca}}) is also available.  Rotations  may be done using factor or components transformations to a target matrix include the standard Promax transformation (\code{\link{Promax}}), a transformation to a cluster target, or to any simple target matrix (\code{\link{target.rot}}) as well as the ability to call many of the GPArotation functions (e.g., oblimin, quartimin, varimax, geomin, ...). Functions for determining the number of factors in a data matrix include Very Simple Structure (\code{\link{VSS}}) and Minimum Average Partial correlation (\code{\link{MAP}}). 

An alternative approach to factor analysis is Item Cluster Analysis (\code{\link{ICLUST}}). This function is particularly appropriate for exploratory scale construction.

There are a number of functions for finding various reliability coefficients (see Revelle and Condon, 2019). These include the traditional  \code{\link{alpha}}  (found for multiple scales and with more useful output by \code{\link{scoreItems}}, \code{\link{score.multiple.choice}}), beta (\code{\link{ICLUST}}) and  both of McDonald's omega coefficients (\code{\link{omega}}, \code{\link{omegaSem}} and  \code{\link{omega.diagram}}) as well as Guttman's six estimates of internal consistency reliability (\code{\link{guttman}}) and the six measures of Intraclass correlation coefficients (\code{\link{ICC}}) discussed by Shrout and Fleiss are also available.  

Multilevel analyses may be done by \code{\link{statsBy}} and  \code{\link{multilevel.reliability}}.

The \code{\link{scoreItems}}, and \code{\link{score.multiple.choice}} functions may be used to form single or multiple scales from sets of dichotomous, multilevel, or multiple choice items by specifying scoring keys.  \code{\link{scoreOverlap}} correct interscale correlations for overlapping items, so that it is possible to examine hierarchical or nested structures.

Scales can be formed that best predict (after cross validation) particular criteria using \code{\link{bestScales}} using unit weighted or correlation weights. This procedure, also called  the \code{\link{BISCUIT}} algorithm (Best Items Scales that are Cross validated, Unit weighted, Informative, and Transparent) is a simple alternative to more complicated machine learning algorithms.


Additional functions make for more convenient descriptions of item characteristics  include 1 and 2 parameter Item Response measures.  The \code{\link{tetrachoric}}, \code{\link{polychoric}} and \code{\link{irt.fa}} functions are used to find 2 parameter descriptions of item functioning. \code{\link{scoreIrt}}, \code{\link{scoreIrt.1pl}} and \code{\link{scoreIrt.2pl}} do basic IRT based scoring.

 A number of procedures have been developed as part of the Synthetic Aperture Personality Assessment (SAPA \url{https://www.sapa-project.org/}) project.  These routines facilitate forming and analyzing composite scales equivalent to using the raw data but doing so by adding within and between cluster/scale item correlations. These functions include extracting clusters from factor loading matrices (\code{\link{factor2cluster}}), synthetically forming clusters from correlation matrices (\code{\link{cluster.cor}}), and finding multiple  ((\code{\link{setCor}}) and partial ((\code{\link{partial.r}}) correlations from correlation matrices.
 
\code{\link{setCor}} and \code{\link{mediate}} meet the desire to do regressions and mediation analysis from either raw data or from correlation matrices.  If raw data are provided, these functions can also do  moderation analyses.
 
 Functions to generate simulated data with particular structures include \code{\link{sim.circ}} (for circumplex structures), \code{\link{sim.item}} (for general structures) and \code{\link{sim.congeneric}} (for a specific demonstration of congeneric measurement).  The functions \code{\link{sim.congeneric}}  and \code{\link{sim.hierarchical}} can be used to create data sets with particular structural properties. A more general form for all of these is \code{\link{sim.structural}} for generating general structural models.  These are discussed in more detail in the vignette (psych_for_sem).

 Functions to apply various standard statistical tests include \code{\link{p.rep}} and its variants for testing the probability of replication, \code{\link{r.con}} for the confidence intervals of a correlation, and \code{\link{r.test}} to test single, paired, or sets of correlations. 
 

 
 In order to study diurnal or circadian variations in mood, it is helpful to use circular statistics.  Functions to find the circular mean (\code{\link{circadian.mean}}), circular (phasic) correlations (\code{\link{circadian.cor}}) and the correlation between linear variables and circular variables (\code{\link{circadian.linear.cor}}) supplement a function to find the best fitting phase angle (\code{\link{cosinor}}) for measures taken with a fixed period (e.g., 24 hours).
 
 A dynamic model of personality and motivation (the Cues-Tendency-Actions model) is include as (\code{\link{cta}}.
 
 A number of useful helper functions allow for data input (\code{\link[psychTools]{read.file}}), and data manipulation \code{\link{cs}} and \code{\link[psychTools]{dfOrder}},

 
 The most recent development version of the package is always available for download as a \emph{source} file from the repository at the PMC lab:
 
   install.packages("psych", repos = "https://personality-project.org/r/", type="source").
   
    This will provide the most recent version for PCs and Macs. 
 
}
\details{Three vignettes (overview.pdf and psych_for_sem.pdf) are useful introductions to the package. They may be found as vignettes in R or may be downloaded from 
\url{https://personality-project.org/r/psych/intro.pdf} \url{https://personality-project.org/r/psych/overview.pdf} and \url{https://personality-project.org/r/psych/psych_for_sem.pdf}.  In addition, there are a number of "HowTo"s available at \url{https://personality-project.org/r/} 

The more important functions in the package are for the analysis of multivariate data, with an emphasis upon those functions useful in scale construction of item composites. However, there are a number of very useful functions for basic data manipulation including  
\code{\link[psychTools]{read.file}}, \code{\link[psychTools]{read.clipboard}}, \code{\link{describe}},  \code{\link{pairs.panels}}, \code{\link{error.bars}} and \code{\link{error.dots}}) which are useful for basic data entry and descriptive analyses.



When given a set of items from a personality inventory, one goal is to combine these into higher level item composites. This leads to several questions:

1) What are the basic properties of the data?  \code{\link{describe}} reports basic summary statistics (mean, sd, median, mad, range,  minimum, maximum, skew, kurtosis, standard error) for vectors, columns of matrices, or data.frames. \code{\link{describeBy}} provides descriptive statistics, organized by one or more grouping variables. \code{\link{statsBy}} provides even more detail for data structured by groups including within and between correlation matrices, ICCs for group differences, as well as basic descriptive statistics organized by group.

 \code{\link{pairs.panels}} shows scatter plot matrices (SPLOMs) as well as histograms and the Pearson correlation for scales or items. \code{\link{error.bars}} will plot variable means with associated confidence intervals. \code{\link{errorCircles}} will plot confidence intervals for both the x and y coordinates.  \code{\link{corr.test}} will find the significance values for a matrix of correlations. \code{\link{error.dots}} creates a dot chart with confidence intervals.

2) What is the most appropriate number of item composites to form? After finding  either standard Pearson correlations, or finding tetrachoric or polychoric correlations,  the dimensionality of the correlation matrix may be examined. The number of factors/components problem is a standard question of factor analysis, cluster analysis, or principal components analysis. Unfortunately, there is no agreed upon answer. The Very Simple Structure (\code{\link{VSS}}) set of procedures has been proposed as on answer to the question of the optimal number of factors.  Other procedures (\code{\link{VSS.scree}},  \code{\link{VSS.parallel}},  \code{\link{fa.parallel}}, and \code{\link{MAP}})  also address this question.  \code{\link{nfactors}} combine several of these approaches into one convenient function. Unfortunately, there is no best answer to the problem.

3) What are the best composites to form?  Although this may be answered using principal components  (\code{\link{principal}}), principal axis (\code{\link{factor.pa}}) or minimum residual (\code{\link{factor.minres}})  factor analysis (all part of the \code{\link{fa}} function) and to show the results graphically (\code{\link{fa.diagram})}, it is sometimes more useful to address this question using cluster analytic techniques.  Previous versions of \code{\link{ICLUST}} (e.g., Revelle, 1979)  have been shown to be particularly successful at forming maximally consistent and independent item composites.  Graphical output from \code{\link{ICLUST.graph}} uses the Graphviz dot language and allows one to write files suitable for Graphviz.  If Rgraphviz is available, these graphs can be done in R.

Graphical organizations of cluster and factor analysis output can be done using  \code{\link{cluster.plot}} which plots items by cluster/factor loadings and assigns items to that dimension with the highest loading.  

4) How well does a particular item composite reflect a single construct?  This is a question of reliability and general factor saturation.  Multiple solutions for this problem result in (Cronbach's) alpha (\code{\link{alpha}}, \code{\link{scoreItems}}), (Revelle's) Beta (\code{\link{ICLUST}}), and (McDonald's) \code{\link{omega}} (both omega hierarchical and omega total). Additional reliability estimates may be found in the \code{\link{guttman}} function.

This can also be examined by applying \code{\link{irt.fa}} Item Response Theory techniques using factor analysis of the \code{\link{tetrachoric}} or \code{\link{polychoric}} correlation matrices and converting the results into the standard two parameter parameterization of item difficulty and item discrimination.  Information functions for the items suggest where they are most effective.

5) For some applications, data matrices are synthetically combined from sampling different items for different people.  So called Synthetic Aperture Personality Assessement (SAPA) techniques allow the formation of large correlation or covariance matrices even though no one person has taken all of the items. To analyze such data sets, it is easy to form item composites based upon the covariance matrix of the items, rather than original data set.  These matrices may then be analyzed using a number of functions (e.g., \code{\link{cluster.cor}},   \code{\link{fa}}, \code{\link{ICLUST}}, \code{\link{principal}}, \code{\link{mat.regress}}, and \code{\link{factor2cluster}}.

6) More typically, one has a raw data set to analyze. \code{\link{alpha}} will report several reliablity estimates as well as item-whole correlations for items forming a single scale, \code{\link{score.items}} will score data sets on multiple scales, reporting the scale scores, item-scale and scale-scale correlations, as well as coefficient alpha,  alpha-1 and G6+. Using a `keys' matrix (created by \code{\link{make.keys}} or by hand), scales can have overlapping or independent items. \code{\link{score.multiple.choice}} scores multiple choice items or converts multiple choice items to dichtomous (0/1) format for other functions. 

7) In addition to classical test theory (CTT) based scores of either totals or averages, 1 and 2 parameter IRT based scores may be found with \code{\link{scoreIrt.1pl}}, \code{\link{scoreIrt.2pl}} or more generally \code{\link{scoreIrt}}. Although highly correlated with CTT estimates, these scores take advantage of different item difficulties and are particularly appropriate for the problem of missing data. 

8) If the data has a multilevel structure (e.g, items nested within time nested within subjects) the \code{\link{multilevel.reliability}} aka \code{\link{mlr}} function will estimate generalizability coefficients for data over subjects, subjects over time, etc. \code{\link{mlPlot}} will provide plots for each subject of items over time. \code{\link{mlArrange}} takes the conventional wide output format and converts it to the long format necessary for some multilevel functions. Other functions useful for multilevel data include \code{\link{statsBy}} and \code{\link{faBy}}.    

An additional set of functions generate simulated data to meet certain structural properties. \code{\link{sim.anova}} produces data simulating a 3 way analysis of variance (ANOVA) or linear model with or with out repeated measures. \code{\link{sim.item}} creates simple structure data,  \code{\link{sim.circ}} will produce circumplex structured data,  \code{\link{sim.dichot}} produces circumplex or simple structured data for dichotomous items.  These item structures are useful for understanding the effects of skew, differential item endorsement on factor and cluster analytic soutions.  \code{\link{sim.structural}} will produce correlation matrices and data matrices to match general structural models. (See the vignette).


When examining personality items, some people like to discuss them as representing items in a two dimensional space with a circumplex structure.  Tests of circumplex fit \code{\link{circ.tests}} have been developed.  When representing items in a circumplex, it is convenient to view them in \code{\link{polar}} coordinates. 

Additional functions for testing the difference between two independent or dependent correlation \code{\link{r.test}}, to find the \code{\link{phi}} or \code{\link{Yule}} coefficients from a two by table, or to find the confidence interval of a correlation coefficient.    

Many data sets are included: \code{\link[psychTools]{bfi}} represents 25 personality items thought to represent five factors of personality, \code{\link[psychTools]{ability}} has 14 multiple choice iq items. \code{\link{sat.act}} has data on self reported test scores by age and gender. \code{\link[psychTools]{galton} }  Galton's data set of the heights of parents and their children. \code{\link[psychTools]{peas}} recreates the original Galton data set of the genetics of sweet peas.  \code{\link[psychTools]{heights}} and
\code{\link[psychTools]{cubits}} provide even more Galton data, \code{\link[psychTools]{vegetables}} provides the Guilford preference matrix of vegetables.  \code{\link[psychTools]{cities}} provides airline miles between 11 US cities (demo data for multidimensional scaling).

\tabular{ll}{
Package: \tab psych\cr
Type: \tab Package\cr
Version: \tab 2.2.5\cr
Date: \tab 2022--May--10\cr
License: \tab GPL version 2 or newer\cr
}
Partial Index:



\link{psych}      A package for personality, psychometric, and psychological research.\cr
 

Useful data entry and descriptive statistics\cr
\tabular{ll}{
\link[psychTools]{read.file}       \tab      search for, find, and read from file\cr
\link[psychTools]{read.clipboard}       \tab      shortcut for reading from the clipboard\cr
\link[psychTools]{read.clipboard.csv}   \tab      shortcut for reading comma delimited files from clipboard \cr
\link[psychTools]{read.clipboard.lower}  \tab            shortcut for reading lower triangular matrices from the clipboard\cr
\link[psychTools]{read.clipboard.upper}   \tab           shortcut for reading upper triangular matrices from the clipboard\cr
\link{describe}        \tab           Basic descriptive statistics useful for psychometrics\cr
\link{describe.by}    \tab            Find summary statistics by groups\cr
\link{statsBy}       \tab            Find summary statistics by a grouping variable,  including within and between correlation matrices. \cr 
\link{mlArrange}        \tab          Change multilevel data from wide to long format\cr 
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
\link{fa}              \tab         Combined function for principal axis, minimum 
residual,  weighted least squares, \cr
\tab and maximum likelihood factor analysis\cr
\link{factor.pa}        \tab        Do a principal Axis factor analysis   (deprecated)\cr
\link{factor.minres}    \tab        Do a minimum residual factor analysis (deprecated)\cr
\link{factor.wls}    \tab        Do a weighted least squares factor analysis (deprecated)\cr
\link{fa.graph}        \tab         Show the results of a factor analysis or principal components analysis graphically\cr
\link{fa.diagram}     \tab          Show the results of a factor analysis without using Rgraphviz \cr
\link{fa.sort}        \tab          Sort a factor or principal components output \cr
\link{fa.extension}   \tab          Apply the Dwyer extension for factor loadingss \cr
\link{principal}       \tab         Do an eigen value decomposition to find the principal components of a matrix\cr
\link{fa.parallel}     \tab         Scree test and Parallel analysis \cr
\link{fa.parallel.poly}     \tab     Scree test and Parallel analysis for polychoric matrices \cr
\link{factor.scores}     \tab           Estimate factor scores given a data matrix and factor loadings \cr
\link{guttman}         \tab          8 different measures of reliability (6 from Guttman (1945) \cr
\code{\link{irt.fa}}   \tab          Apply factor analysis to dichotomous items to get IRT parameters \cr
\code{\link{iclust}}   \tab          Apply the ICLUST algorithm\cr
\link{ICLUST.graph}    \tab          Graph the output from ICLUST using the dot language\cr
\link{ICLUST.rgraph}   \tab          Graph the output from ICLUST using rgraphviz \cr
\link{kaiser}     \tab           Apply kaiser normalization before rotating \cr
\link{polychoric}        \tab          Find the polychoric correlations for items  and find item thresholds\cr
\link{poly.mat}        \tab          Find the polychoric correlations for items (uses J. Fox's hetcor) \cr
\link{omega}           \tab          Calculate the omega estimate of factor saturation (requires the GPArotation package)\cr
\link{omega.graph}      \tab         Draw a hierarchical or Schmid Leiman orthogonalized solution (uses Rgraphviz) \cr 
\link{partial.r}        \tab         Partial variables from a correlation matrix \cr
\link{predict}          \tab         Predict factor/component scores for new data \cr
\link{schmid}            \tab        Apply the Schmid Leiman transformation to a correlation matrix\cr
\link{score.items}       \tab        Combine items into multiple scales and find alpha\cr
\link{score.multiple.choice}  \tab   Combine items into multiple scales and find alpha and basic scale statistics\cr
\link{set.cor}             \tab       Find Cohen's set correlation between two sets of variables \cr
\link{smc}             \tab          Find the Squared Multiple Correlation (used for initial communality estimates)\cr
\link{tetrachoric}            \tab    Find tetrachoric correlations and item thresholds \cr
\link{polyserial}            \tab    Find polyserial and biserial correlations for item validity studies \cr
\link{mixed.cor}          \tab       Form a correlation matrix from continuous, polytomous, and dichotomous items \cr
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
\link{omegaSem}           \tab       Calculate the omega estimates of reliability  using a Confirmatory model (requires the sem package)\cr
\link{ICC}          \tab             Intraclass correlation coefficients  \cr
\link{score.items}       \tab        Combine items into multiple scales and find alpha\cr
\link{glb.algebraic} \tab  The greates lower bound found by an algebraic solution (requires Rcsdp).  Written by  Andreas Moeltner \cr
}

Procedures particularly useful for Synthetic Aperture Personality Assessment\cr


\tabular{ll}{
\link{alpha}      \tab         Find coefficient alpha and Guttman Lambda 6 for a scale (see also \link{score.items})\cr
\link{bestScales} \tab A bootstrap aggregation function for choosing most predictive unit weighted items \cr
\link{make.keys}        \tab         Create the keys file for score.items or cluster.cor            \cr
\link{correct.cor}      \tab        Correct a correlation matrix for unreliability\cr 
\link{count.pairwise}   \tab        Count the number of complete cases when doing pair wise correlations\cr
\link{cluster.cor}       \tab        find correlations of composite variables from larger matrix\cr
\link{cluster.loadings}  \tab        find correlations of items with  composite variables from a larger matrix\cr
\link{eigen.loadings}    \tab        Find the loadings when doing an eigen value decomposition\cr
\link{fa}        \tab        Do a minimal residual or principal axis factor analysis and estimate factor scores\cr
\link{fa.extension}        \tab      Extend a factor analysis to a set of new variables\cr
\link{factor.pa}        \tab        Do a Principal Axis factor analysis and estimate factor scores\cr
\link{factor2cluster}    \tab       extract cluster definitions from factor loadings\cr
\link{factor.congruence}  \tab      Factor congruence coefficient\cr
\link{factor.fit}       \tab        How well does a factor model fit a correlation matrix\cr
\link{factor.model}    \tab          Reproduce a correlation matrix based upon the factor model\cr
\link{factor.residuals} \tab        Fit = data - model\cr
\link{factor.rotate}   \tab         ``hand rotate" factors\cr
\link{guttman}         \tab         8 different measures of reliability\cr
\link{mat.regress}     \tab         standardized multiple regression from raw or correlation matrix input\cr
\link{polyserial}       \tab      polyserial and biserial correlations with massive missing data\cr
\link{tetrachoric}      \tab    Find tetrachoric correlations and item thresholds \cr
}

Functions for generating simulated data sets \cr
\tabular{ll}{
\link{sim}               \tab       The basic simulation functions \cr
\link{sim.anova}         \tab       Generate 3 independent variables and 1 or more 
dependent variables for demonstrating ANOVA \cr 
\tab  and lm designs \cr
\link{sim.circ}         \tab        Generate a two dimensional circumplex item structure \cr
\link{sim.item}         \tab        Generate a two dimensional simple structure with 
particular item characteristics \cr
\link{sim.congeneric}   \tab        Generate a one factor congeneric reliability structure \cr
\link{sim.minor}       \tab         Simulate nfact major and nvar/2 minor factors  \cr
\link{sim.structural}   \tab        Generate a multifactorial structural model  \cr
\link{sim.irt}          \tab        Generate data for a 1, 2, 3 or 4 parameter logistic model\cr
\link{sim.VSS}          \tab          Generate simulated data for the factor model\cr
\link{phi.demo}          \tab       Create artificial data matrices for teaching purposes\cr
\link{sim.hierarchical} \tab        Generate simulated correlation matrices with hierarchical or any structure\cr 
\link{sim.spherical}    \tab        Generate three dimensional spherical data (generalization of circumplex to 3 space)\cr
}


Graphical functions (require Rgraphviz) -- deprecated \cr
\tabular{ll}{
\link{structure.graph}  \tab           Draw a sem or regression graph \cr
\link{fa.graph}        \tab           Draw the factor structure from a factor or principal components analysis \cr
\link{omega.graph}      \tab           Draw the factor structure from an omega analysis(either with or without the Schmid Leiman transformation) \cr 
\link{ICLUST.graph}     \tab           Draw the tree diagram from ICLUST  \cr
}


Graphical functions that do not require Rgraphviz \cr
\tabular{ll}{
\link{diagram}           \tab           A general set of diagram functions. \cr 
\link{structure.diagram}  \tab           Draw a sem or regression graph \cr
\link{fa.diagram}        \tab           Draw the factor structure from a factor or principal components analysis \cr
\link{omega.diagram}      \tab           Draw the factor structure from an omega analysis(either with or without the Schmid Leiman transformation) \cr 
\link{ICLUST.diagram}     \tab           Draw the tree diagram from ICLUST  \cr
\link{plot.psych}        \tab           A call to plot various types of output (e.g. from irt.fa, fa, omega, iclust \cr
\link{cor.plot}        \tab            A heat map display of correlations \cr
\link{scatterHist}     \tab            Bivariate scatter plot and histograms \cr
\link{spider}          \tab            Spider and radar plots (circular displays of correlations)
}

Circular statistics (for circadian data analysis) \cr
\tabular{ll}{
\link{circadian.cor}     \tab      Find the correlation with e.g., mood and time of day \cr
\link{circadian.linear.cor} \tab  Correlate a circular value with a linear value \cr
\link{circadian.mean}     \tab      Find the circular mean of each column of a a data set  \cr
\link{cosinor}       \tab            Find the best fitting phase angle for a circular data set \cr
}


 
Miscellaneous functions\cr
\cr
\tabular{ll}{
\link{comorbidity}   \tab             Convert base rate and comorbity to phi, Yule and tetrachoric\cr
\link[psychTools]{df2latex}   \tab               Convert a data.frame or matrix to a LaTeX table \cr
\link{dummy.code}     \tab          Convert categorical data to dummy codes \cr
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
\link{scaling.fits}     \tab         Compares alternative scaling solutions and gives goodness of fits \cr
\link{scrub}   \tab            Basic data cleaning \cr
\link{tetrachor}   \tab            Finds tetrachoric correlations \cr
\link{thurstone}       \tab          Thurstone Case V scaling \cr
\link{tr}               \tab         Find the trace of a square matrix \cr
\link{wkappa}          \tab          weighted and unweighted versions of Cohen's kappa \cr
\link{Yule}            \tab           Find the Yule Q coefficient of correlation \cr
\link{Yule.inv}         \tab          What is the two by two table that produces a Yule Q with set marginals? \cr
\link{Yule2phi}         \tab          What is the phi coefficient corresponding to a Yule Q with set marginals? \cr
\link{Yule2tetra}        \tab          Convert one or a matrix of Yule coefficients to tetrachoric coefficients. \cr
}


Functions that are under development and not recommended for casual use \cr
\tabular{ll}{
\link{irt.item.diff.rasch}  \tab    IRT estimate of item difficulty with assumption that theta = 0\cr
\link{irt.person.rasch}	   \tab     Item Response Theory estimates of theta (ability) using a Rasch like model\cr\cr
}


Data sets included in the psych or psychTools package \cr
\tabular{ll}{
\link[psychTools]{bfi}          \tab           represents 25 personality items thought to represent five factors of personality \cr
\link{Thurstone}     \tab           8 different data sets with a bifactor structure  \cr
\link[psychTools]{cities}        \tab          The airline distances between 11 cities (used to demonstrate MDS) \cr
\link[psychTools]{epi.bfi}        \tab         13 personality scales \cr
\link[psychTools]{iqitems}        \tab         14 multiple choice iq items \cr
\link[psychTools]{msq}            \tab         75 mood items  \cr
\link{sat.act}        \tab          Self reported ACT and SAT Verbal and Quantitative scores by age and gender\cr
\link{Tucker}         \tab          Correlation matrix from Tucker \cr
\link[psychTools]{galton}        \tab            Galton's data set of the heights of parents and their children \cr
\link[psychTools]{heights}       \tab            Galton's data set of the relationship between height and forearm (cubit) length \cr
\link[psychTools]{cubits}        \tab             Galton's data table of height and forearm length \cr
\link[psychTools]{peas}          \tab            Galton`s data set of the diameters of 700 parent and offspring sweet peas \cr
\link[psychTools]{vegetables}     \tab           Guilford`s preference matrix of vegetables (used for thurstone) \cr
}

A debugging function that may also be used as a demonstration of psych.
\tabular{ll}{
\link{test.psych}   \tab          Run a test of the major functions on 5 different data 
sets.  Primarily for development purposes.\cr 
\tab  Although the output can be used as a demo of the various functions.
}

}
\note{Development versions (source code) of this package are maintained at the  repository \url{https://personality-project.org/r/} along with further documentation.   Specify that you are downloading a source package.  
\cr
Some functions require other packages. Specifically, omega and schmid require the GPArotation package,  ICLUST.rgraph and fa.graph require Rgraphviz but have alternatives using the diagram functions.  i.e.:
\cr

\tabular{ll}{
function  \tab     requires\cr
\link{omega}    \tab     GPArotation \cr
\link{schmid}    \tab     GPArotation\cr
\link{ICLUST.rgraph} \tab   Rgraphviz \cr
\link{fa.graph}  \tab      Rgraphviz \cr
\link{structure.graph} \tab Rgraphviz \cr
\link{glb.algebraic} \tab Rcsdp \cr
}

}

\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illiniois \cr
        \url{https://personality-project.org/revelle.html}\cr


Maintainer: William Revelle <revelle@northwestern.edu>

}
\references{A general guide to personality theory and research may be found at the personality-project \url{https://personality-project.org/}. See also the short guide to R at \url{https://personality-project.org/r/}.
In addition, see 

Revelle, W. (in preparation) An Introduction to Psychometric Theory with applications in R. Springer. at \url{https://personality-project.org/r/book/}   

Revelle, W. and Condon, D.M. (2019) Reliability from alpha to omega: A tutorial.  Psychological Assessment, 31, 12, 1395-1411. https://doi.org/10.1037/pas0000754.  \url{https://psyarxiv.com/2y3w9/} Preprint available from PsyArxiv  }

\keyword{package}% __ONLY ONE__ keyword per line
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line


\examples{
#See the separate man pages 
#to test most of the psych package run the following
#test.psych()   
}
