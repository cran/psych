\name{NEWS}
\title{News for Package 'psych'}
\section{Changes in psych version 1.2.9 (2012-08-31)}{
  \subsection{Introduction}{
    \itemize{
      \item Version 1.2.9 is the development release of the psych package. It is available at the repository at \url{http://personality-project.org/r}. The released version is 1.2.8.  The second digit reflects the year (i.e., 2012), the third set the month (i.e., 1.2.8 was released in August of 2012, the last two digits reflect the day of any modifications).
      \item The psych package includes functions and data sets to do classic and modern psychometrics and to analyze personality  and experimental psychological data sets. The psych package has been developed as a supplement to courses in research methods in psychology, personality research, and graduate level psychometric theory. The functions are a supplement to the text (in progress): An introduction to psychometric theory with applications in R. 
      \item Additional functions are added sporadically.
      \item This news file reports changes that have been made as the package has been developed.
    }
  }
  \subsection{To do}{
    \itemize{
      \item Suggestions are welcome, but the current list includes:
      \item Add confirmatory clustering to ICLUST
      \item Get cluster scores in ICLUST -- analogous to factor scores (requested by Ben Shalet)
      \item Add signal detetection function as an option for dichotmous data. (in progress)
      \item Add the ability to create multiple groups in sim.irt functions
      \item Make pairs.panels not change options for graphics window   (fixed, maybe?)
      \item Add labels to factor.plot for more than 2 factors?
      \item Add goodness of fit tests for fa.extension applied to the extension variables. 
      \item Find canonical loadings in set.cor
      
      \item Add omega factor extension figure option (requested by Sylia Wilson)
      \item Fix bug in describe.by for the case of NULL categories for the by variable. (Reported by Nipa Phojanamongkolkij).
      \item Add option to sim.congeneric (and others) to do multiple simulations with same theta values. 
     \item add option to do subject density and item density plot to IRT plot.
      \item note that the df in set.cor should reflect the number of cases actually used
     
      \item add the ability to scale radar plots of raw data from min to max, and add a scale to radar and spider plots (e.g. circular histograms)
      \item add the ability to add labels to lavaan.diagram    
     \item modify structure  diagram so that the error arrows are to the left for lr printing
     \item add the von Neumann test of successive differences
    
    }
  }
  } 
  
  \section{Changes in psych version 1.2.8 (2012-8-25)}{
   \subsection{Additions}{
\itemize{
 \item Updated the overview vignette
 \item Added the ability to do target rotations (TargetQ) in omega.
  \item Improved documentation of fa to explain how oblique SS are found.
  \item A better set of items is used for the iq example (iqitems).
  \item add CVE estimate to fa (see Reise 2010) (had already done this for omega).
  \item Changed ICLUST.sort so that it can sort the loadings of principal (requested by Gudmundur Arnkelsson and then print them with the number of digits requested in the print command. 
  \item Modified score.irt so that subjects who miss all items or pass all items are given an estimate based upon the (product) of the difficulty of the items they miss (pass) and then adjusted based upon half the quantile difference from 0 (if they miss all items) and 100 (if they pass all items).
  \item Modified sim.omega to allow for specifying a general factor.  This allows for tests of not just the bias in the case of no general factor, but also the ability to detect a general factor.  Also modified it to include calls to omegaSem.
   \item modified iclust.diagram so that cluster names can be specified rather than all start with C1 ... Cn  (requested by Michael Kubovy)
   \item Slightly improved the documentation for r.test so it is clear which correlation is which. 
   \item Modified iclust so that the fit statistic is based upon the off-diagonal elements unless otherwise specified.  Use diagonal =TRUE to get fits matching previous analyses. 
   \item Modified the print.psych.iclust function to print out the Root Mean Square Residual correlation.  (It was previously reported in iclust output, but not printed.)
    }
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
       \item  fixed a bug in score.items, such that if the number of subjects is equal to the number of items, no scores were returned.  Reported by Jeromy Anglim. 
    \item Modified factor.stats so that rare condition of an exact fit in omega does not lead to an error.
    \item Fixed plot.irt so that item labels are correct
    \item Modified omegaSem to pass n.obs to sem 
    \item Fixed df2latex so that it will on text variables as well as numeric variables   
   }
   } 
  }
  
  \section{Changes in psych version 1.2.7 (2012-7-31)}{
   \subsection{Additions}{
\itemize{
  \item Improved statsBy to find within and between group correlations.
  \item Added a data set, withinBetween, to graphically show the issue of multilevel correlations.
  \item Added a simulation function to generate multilevel data (sim.multilevel)
  \item Added a function (irt.responses) to plot responses as a function of the latent score for multiple choice alternatives.
  \item Modified progressBar so that it only shows dots when using the console, but not when using a text file for output (i.e., with Sweave).
    }
  
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
       \item fixed plot.poly and plot.irt so that multiple scales for one factor can be plotted.
   }
   } 
  }
  
\section{Changes in psych version 1.2.6 (2012-6-20)}{
   \subsection{Additions}{
\itemize{
   \item Slight modification to cor.plot to allow for better control of multiple plots on the same page. 
   \item Added the ability to put values into the cor.plot output (numbers=TRUE)
   \item Added a function (lowerUpper) to combine two symmetric matrices and output one as the above diagonal entries, and the other as the below diagonal entries of a square matrix. 
   \item Added a function (topBottom) to act like headTail but to do it without ellipses. 
   \item In a gradual switch to camelCase, changed the name of headtail to headTail, and describe.by to describeBy. Warning messages are issued about the change.
   \item Added a function, statsBy, to find summary statistics (means, sample sizes, standard deviations) by a grouping variable.  Similar to describeBy but somewhat easier to use. Also will report (as an option) the pooled within group correlations. 
   \item Added a function to find sample size weighted correlations, cor.wt.  Useful for SAPA types of analysis. 
   \item Modified superMatrix to take list input so that more than two matrices can be combined at once (suggested by David Condon).
   \item Added errorCircles, a greatly improved version of error.crosses. 
   \item Added the option to specify starting communalities in the fa function (requested by David Gosar). 
    }
  
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
\item cor.plot was actually switching row and colnames (but since it is mainly used on symmetric matrices, this was not a problem).  Fixed.
\item principal failed when handling missing data in data.frame. Fixed. (reported by Neil Stewart)

   }
   } 
  }
  
   \section{Changes in psych version 1.2.4 (2012-4-30)}{
   \subsection{Additions}{
\itemize{
 \item (Note that 1.2.3 was not released but was replaced with 1.2.4)
 \item fixed ICLUST so that it converts covar matrices to correlations. 
 \item Modified pairs.panels so that the method of correlation can be specified (suggested by Carsten Dormann)
 \item added error values to the structure.diagram/lavaan diagram as an option (suggested by José Luis Cañadas )
 \item Changed the oblique.scores option in fa to be FALSE (use the structure the matrix) following a request by Niels Waller.
   \item Fixed the size of the boxes in dia.rect so that omega.diagram and  fa.diagram provide reasonable size boxes for long names (requested by Erich Studerus)
    \item modified output of corr.test so that it announces whether or not it is adjusting the probability values for multiple tests.
  \item added the ability to export dot.plot commands without requiring Rgraphviz. (requested by Erich Studerus).  Although just implemented for fa.graph, omega.diagram functions can be exported using fa.graph applied to the schmid$sl object from omega.
  \item Added 3 options to the skew and kurtosis functions, and thus to the describe and describe.by functions as well.  (Suggested by Bruce Dudek).
  \item Added the ability to plot general factors (ala omega) in the fa.diagram function.
  \item Added the function lowerCor to find correlations (with pairwise deletion as a default) and print out the lowerMat of the results.
    }
  
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
\item A documentation error in r.test was fixed (reported by Nicholas S. Holtzman).  The example of Steiger case B was incorrectly ordered and thus led to incorrect results. Specifying by name (rather than location) is now recommended.
 \item In fa, the rotation was not returned, fixed.  (reported by Niels Waller)
 \item predict.psych would not work for a single case (reported by Jonathan Williams).  Fixed.  Also made the documentation a bit clearer to explain the way the factor scores are standardized.
 \item By adding a global-local parameter to mixed.cor it now seems to work on those data sets with badly distributed polytomous data. (Reported by Eric Smith and David Condon.)
   }
   } 
  }
  
  \section{Changes in psych version 1.2.1 (2012-1-31)}{
   \subsection{Additions}{
  \itemize{
    \item Modified scrub to be a bit easier to use (suggested by Josh Wilt)
   \item Added a check in mixed.cor so that if the data are incorrect, announce it first, rather than waiting until late in the process. 
   \item Added a progress indicator to particularly slow functions (e.g. tetrachoric, polychoric, biserial).  This will show up when using these functions in other functions (such as irt.fa or mixed.cor). 
   \item Modified cluster.cor so that in the case of correlations that are NA, it will not find smcs of the raw correlation matrix.  In addition, it will impute scale level correlations based upon the average between scale correlation (previously it was reporting an NA).  Requested by David Condon.
   \item Added the helper function lowerMat to print the lower submatrix.  Added this as an option to the print.psych function.
   \item Modified mixed.cor so that it will automatically (if desired) find the Pearson, polychoric, and tetrachoric correlation matrices and then reorganize the matrix back to the original form.
   \item Minor modifications to tetrachoric and polychoric to give slightly more helpful error messages in the case of missing data and corrections for continuity.
   \item Minor modifications to summary to give cleaner summary stats for fa and irt.fa analyses.
   \item Added a residual and plot.residuals function to plot either qq plots of the residuals, or a cor.plot of the residuals.
  }
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
   \item Fixed a few references that were working but doing partial matching in print.psych, test.psych, pairs.panels (e.g, in print.psych digit=digits was changed to digits=digits, test.psych nf was changed to nfactors, col was changed to col.smooth in pairs.panels). This had been flagged by linux checking as a note about "partial matching" for a long time. Finally figured out what it meant.
 \item corrected biplot.psych to properly take the xlim and ylim parameters (reported by Andreas M. Brandmaier)
 \item Correction to tetrachoric (and related functions) to return NA in the case of no cases for x or y.  (Reported by David Condon)
 \item Fixed fa so that the score.oblique option actually works. (Reported by Jessica Jaynes).
   }
   } 
  }
  
  \section{Changes in psych version 1.1.12 (2011-12-30)}{
   \subsection{Additions}{
  \itemize{ 
 
  \item Modified scrub to be more general.
  \item Added factor structure as an object in the fa function and modified the print.psych.fa function to note that it is the pattern matrix being reported. 
  \item Modified the ICLUST output so that the eigenvalues reflect the independent contribution of each cluster. 
  \item Added a comment and an option to ICLUST.sort about the order in which clusters are sorted (suggested by Gudmundur Arnkelsson).
  \item Added the rug option to pairs.panels.  This will, by default draw a rug underneath the histograms.
  \item Added a covar option to principal to find principal components from covariance matrices as well as correlation matrices.
  \item Cleaned up the objects returned from fa.parallel.poly to more closely match those from fa.parallel. 
  \item Added the ability to find and report residual correlations in set.cor. Minor formatting corrections to set.cor and mat.regress.
  
  }
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
  \item  Fixed a bug in fa such that the oblique.scores option did not work. 
  \item  fa.parallel.poly was ignoring the fm parameter and was not printing correctly.
   }
   } 
  }


\section{Changes in psych version 1.1.1122 (2011-11-22)}{
   \subsection{Additions}{
  \itemize{ 
  \item Added various rotation options to principal to match those in fa.
  \item Added "components" as an option to factor score to  score principal components using just the component loadings.  This is most appropriate when not rotating the components.
  \item Added the Harman.5 (socio-demographic) data set to allow comparisons with SAS. 
  }
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
  \item The way component scores were found had changed from "regression" to "tenBerge" without being documented.  Changed the documentation and added the method of scoring as an option.  Problem reported by Alexander Weiss.  
   }
   } 
  }
  
\section{Changes in psych version 1.1.1111 (2011-11-11)}{
  \subsection{Additions}{
  \itemize{ 
  \item This is basically a rerelease of version 1.1.10 trying to solve a problem with 32 bit machines and one test that fails. Originally call 1.11.01 but that one still failed the tests.
  \item Added ... as an option in multi.hist (following a suggestion by Rui Barradas)
  \item Added the all=FALSE option to test.psych.  This turns off one test to avoid a problem in the testing on Solaris 32 bit machines.
  }
  }
  
   \subsection{Bugs Fixed}{
   \itemize{
   \item factor.scores and factor.stats were failing in tests for Linux operating systems, but not Macs or PCs.  Added some protection against complex eigen values.  
   \item fixed problem in irt.scores reported by David Condon 
   \item fixed irt.fa so it will properly plot for polytomous items
   }
   }
  
  }
\section{Changes in psych version 1.1.10 (2011-10-15)}{
  \subsection{Additions}{
  \itemize{ 
  \item Started to modify all calls to sd and mean to use apply to respond to their deprecation for data.frames and matrices
  \item Changed the release numbering system to reflect year and month of release.  
  \item Added a parcels function to form item parcels of size 2 or size 3.
  \item Added a jitter option to the factor.plot function.
  \item Added 3 more factor scoring options to the factor.scores function and changed the default scoring option in fa to be the tenBerge
  \item Added the kaiser function to do kaiser normalization.
  \item Added the df2latex function to make LaTeX tables.
    \item Added plot.circular to draw radar and spider plots
  	\item Added a short function to create dummy codes (dummy.code)
    \item Added a score.irt function to find IRT based scores. 
    \item Added a table of information by attribute level to be (silently) produced  in the plot(irt.fa) function.
    \item Added bifactor and biquartimin as  rotation options to fa (based on the Jennrich and Bentler 2011 papers).
    \item Implemented an improvement to print.psych function to use switch.  (Following a suggestion by Joshua Wiley ). This does not change the functionality, but makes the code easier to debug and to change.
     \item Added an option to alpha  to automatically flip items if they seem to be negatively correlated with total score.  (Suggested by Jeremy Miles).
     \item Serious modification to cor.plot to allow for colorRamp data.  Makes much prettier correlation plots. (Suggested by David Condon and Joshua Wilt).
 	 \item Modified sim and sim.simplex to create State Trait Auto Regressive simplex structures (suggested by Deepika Anand).
     \item Added a correct option to irt.fa (passes correct to tetrachoric).
     \item Added the ability for omega to be rerun on the same correlation matrix found by previous omega runs.  This speeds up analyses of large matrices but particularly that of tetrachoric/polychoric matrices.
     \item Modified headtail so that it will not choke on mixed numeric and string data
     \item Added a smoothing option to the tetachoric correlation function to smooth the resulting matrix if it is not positive definite.
     \item Added cor.smooth to do a principal components based smoothing for correlations.
     \item Added global=TRUE option to polychoric (and tetrachoric) correlations. If global=FALSE, each correlation is found using the pairwise taus.  This seems to match John Fox's polycor function.  If global=TRUE, the taus are the datawise values. This will differ from polycor in the case of a great deal of missing data.  Somewhat faster and perhaps more appropriate given that we are trying to find a global solution. 
     \item Added a correction for multiple tests to the p values in the \code{\link{corr.test}} function and added a new function to just report probabilities \code{\link{corr.p}}
     }
     }
   \subsection{Bugs Fixed}{
   \itemize{
   \item Fixed a very subtle bug in score.items for the case of impute="none" where some subjects have no non-missing data. This affected the estimate of the average number of items (num.ob.item) per scale. Reported by Ben Schalet.
    \item Fixed bug in factor.stats that gave much too narrow confidence intervals for RMSEA (Thanks to Rick Zinbarg and Ashley Kendall for pointing out the problem).
    \item revised factor extension to properly handle oblique solutions
    \item fa now reports the call correctly -- since introduction of iterations it had reported the wrong call
    \item Fixed mixed.cor so that it can handle several special cases.
    \item Perhaps finally fixed pairs.panels so it leaves us with the original parameter settings
    \item Fixed a problem in print.psych such that mat.regress output would not print.  (Reported by Ben Schalet.)
    \item Fixed plot.poly so that it correctly plots the ICC results, particularly for negatively keyed items.
    \item Fixed fa so that the covar option works.  (It had been killed when iterations were introduced.)
    \item Modified fa, fa.irt, factor.stats, factor.scores so that they handle bad data more gracefully.  
   }
  }
}
\section{Changes in psych version 1.0.98 (2011-06-10)}{
  \subsection{Additions}{
  \itemize{
 	 \item Improvements to set.cor to report the canonical correlations and their chi squares.  Also added a plot function to set.cor so that it will plot the eigen values of the canonical correllations.
 	 \item Added the cushny data set to show the original drug data from student.
     }
  }
   \subsection{Bugs Fixed}{
   \itemize{
    \item Resolved why the denomiator df in set.cor does not precisely match the Cohen example.  This seems to be because Systat is doing an unbiased estimator and is thus multiplying by (u+v)/v.
   }
  }
  }  
\section{Changes in psych version 1.0.97 (2011-05-15)}{
  \subsection{Additions}{
  \itemize{
  \item Added a polychoric and tetrachoric option to fa.parallel
  \item Added a keys option to plot.psych  (plot.poly and plot.irt) to allow for plotting information functions for selected items while retaining the factor loadings from the complete set. (Requested by Ben Schalet.)
  \item Added a summary of the graphic output for irt.fa that includes the average information for each item as well as where the item is most informative. (Suggested by David Condon.)
  \item Integrated the interactions between fa, fa.poly, irt.fa, fa.diagram, plot, etc. to make commands more consistent.
  \item added a factor analysis for tetrachoric or polychoric correlations (fa.poly) that will work without doing an irt analysis using irt.fa.
  \item Added bootstrapped confidence intervals for loadings and correlations in factor  analysis (fa), for reliability estimates in omega (omega) and fa.poly.(Originally requested by Steve Miller and Renee Engeln-Maddox polychoric correlations, extended to the omega case.)
   \item Modified fa.diagram to include graphic output for fa.extension.
  \item Added Cohen's set correlation to the mat.regress function. (Requested by Amanda Uliaszek.)
  \item Changed the order of the parameters in the mat.regress function to be more consistent with the y~x notation of most regression functions.
  \item  Added an oblique.scores option to principal.  This means that if oblique rotations are done (not the default) then the scores from principal will be oblique as well.  
   }
  }
  \subsection{Bugs Fixed}{
   \itemize{
    \item None yet
   }
  }
  }
  \section{Changes in psych version 1.0.96 (2011-04-04)}{
  \subsection{Additions}{
  \itemize{
   \item Factor extension (fa.extension) now can extend an omega as well as a normal factor analysis.
    \item Modified fa so that the factor scores are found for the oblique factors rather than the orthogonal factors. This now produces factor scores that have roughly the same correlations as do the factors.  Factor score statistics remain the same. This leads to a different set of factor scores than factanal finds, for it reports the factor scores for the orthogonal factors. 
   
    \item Modified response.frequencies to incorporate a very nice suggestion by Joshua Wiley to identify unique item responses that are not necessarily integers.  
   }
  }
  \subsection{Bugs Fixed}{
   \itemize{
      \item Modified fa.extension to correctly find extended oblique loadings.
      \item check added to omega for bizarre case of 1 factor extracted through PC.  (This should not be done, for it makes no sense, but nonetheless, it was leading to omega_h > omega_t.) A warning for this condition is now issued and omega_h = omega_t.  Reported by Alison Lewis.
   }
  }
  }
\section{Changes in psych version 1.0.95 (2011-03-30)}{
  \subsection{Additions}{
 \itemize{
 \item Fixed News so that it works as a NEWS.rd file
 \item Added an example to msq
 \item Added a function to perform factor extension (fa.extension)
 \item Added the Dwyer and Gorsuch data sets as examples of fa.extension.
 \item Added lavaan.diagram to the set of diagram functions.  Modified to work with lavaan 0.4-7
 \item Added scatter.hist to draw a scatter plot and associated histograms and densities.  
\item Modified score.items so that it reports the number of missing responses for each person for each scale.  It had previously reported the total number of missing responses but did not break this down by scale.  Requested by Ryne Sherman.  Changed the documentation to reflect that missing is calculated.
\item Modified error.bars.by so that it can not draw lines (lines=FALSE) if desired. (Requested by Ryne Sherman.)
\item Added to  error.bars.Rd to describe how to supply means and standard errrors. 
\item Made various data files available as separate rather than collated. (e.g., Thurstone, Bechtoldt,Holzinger) in the bifactor data set (they remain there as well). (Requested by Michael Friendly.) 
\item Changed VSS.scree and scree to provide a few more options.
\item Cleaned up circ.tests to be more compatible with rest of psych.
\item Clean up output from schmid to more match other functions.
\item Used tools::compactPDF on the vignettes to make them smaller.  Also switched to using png instead of pdf for some of the graphics in the vignettes. This reduced the size drastically.
  }
 }
\subsection{Bugs Fixed}{ 
\itemize{
	\item fixed error.bars.by so that it draws the alpha level instead of 2 alpha level error bars.
	\item fixed fa for special case of 3 x 3 matrix with correlations of exactly .5 which would lead to a non-conforming array problem
 	\item fixed partial.r to not round output.  For rounded output, the print.psych function defaults to 2 digits (reported by David Freedman)     
 	\item Finally fixed the output of print.psych.fa so that it labels sorted factors independent of cut (Reported by Reinhold Hatzinger )
   }
}
}
\section{Changes in psych version 1.0.94 (2011-01-01)}{
  \subsection{Additions}{
    \itemize{
      \item Eliminated references to polychor in Yule2poly, 
    }
  }
  \subsection{Bugs Fixed }{
    \itemize{
      \item fixed read.clipboard.upper so that it will read labels
      \item fixed  mat.regress so that it can work with a single predictor
    }
  }
}
\section{Changes in psych version 1.0.93 (2010-12-22)}{
  \subsection{Additions}{
    \itemize{
      \item Added mixed.cor to do mixtures of Pearson, polychoric and tetrachoric correlations
      \item Added legends to error.bars.by  bar graph
      \item Modified score.items so it will work on correlation matrices as well as raw data.
      \item Improved the plot.irt and plot.poly functions to plot irt output.  Called by plot.psych or generic plot call.
      \item Added a trivial little function, bi.bars to do paired histograms
      \item Modified irt.fa to allow for multiple factors
      \item Modified scrub to be a general recoding function
      \item Modified error.bars.by to allow for ... parameters to be passed into axis
      \item Modified diagram so it will dispatch a number of diagram functions (fa, prinicpal, omega, iclust, lavaan).
      \item Added a biplot.psych function to do biplots of fa or principal results.
      \item Added an unbiased estimate of kurtosis in the kurtosi function.
      \item Added a predict function to fa and principal (predict.fa and predict.principal, both called from predict.psych) 
      \item Added lavaan.diagram to draw the output from a cfa or sem from the lavaan package.
      \item Added a function to calculate Mardia's test for multivariate normality
      \item Added lengends and improved the examples for pairs.panels.by
      \item Removed the "short" option in score.items
      \item Added the ability to score scales with greatly missing data (for SAPA). Modified the impute option to include = "none"
      \item Added the term KR20 and KR21 to documentation of alpha
      \item Modified pairs.panels so that with the cor option it can draw regressions and give correlations.
      \item Modified the documentation for error.crosses to give a better example.
      \item Added the ability to read rownumbers (variable names?) in first column of read.clipboard.lower and read.clipboard.upper
      \item Tweaked the zlim default for cor.plot to be c(-1,1) rather than c(0,1)
      \item Minor adjustments to Yule and phi.demo to call tetrachoric instead of polycor, thus making the package work with fewer dependencies.
      \item Cleaned up the fa summary output 
    }
  }
  \subsection{Bugs Fixed }{
    \itemize{
      \item Fixed dia.curved.arrow to draw one headed curved arrows (as documented)
      \item Fixed a problem with fa.print option sort  -- Although giving the correct communalities and uniquenesses in the unsorted option, when sorted, it was not sorting h2 or u2. 
    }
  }
}
\section{Changes in psych version 1.0.92 (2010-09-22)}{
  \subsection{Additions}{
    \itemize{
      \item Added omega total from sem to the omegaSem function.
    }
  }
  \subsection{Bugs Fixed }{
    \itemize{
      \item error.bars.by was ignoring the user supplied ylim.  (reported by Helena Chui)  Fixed.
      \item Was drawing incorrect  sl solution in omega.graph (reported by Mark Difford).fixed.
      \item omega.Rd and schmid.Rd incorrectly referred to fm=ml as fm=mle. Fixed.
    }
  }
}
\section{Changes in psych version 1.0.91 (2010-09-15)}{
  \subsection{Additions}{
    \itemize{
      \item Added iclust and vss as an alternative names to ICLUST and vss to be consistent with most naming conventions.
      \item Added the simple utility scrub to do data cleaning.
      \item Added biserial and polyserial to do biserial/polyserial correlations (not to be confused with point biserial)
      \item Modified score.multiple.choice so that it handles missing data without using imputation
      \item Changed the bfi data set to include more subjects and to correctly define the Openness scale.
      \item Modified score.multiple.choice and score.items to report response frequencies for each item (frequencies based upon number of non-missing responses).
      \item Added the burt data set of 11 emotional variables (from Burt, 1915) 
      \item Added draw.tetra to draw tetrachoric correlations for fixed tau values (teaching demo)
      \item Added plot.irt and plot.poly to the plot.psych function. These have options to draw the item characteristic function, the item information function, and the test information.
      \item Added irt.fa to do factor analysis of dichotomous or discrete items by first forming a tetrachoric or polychoric correlation matrix. Then convert the results to IRT form. 
      \item Added tetrachoric and polychoric functions to find tetrachoric/polychoric correlations as well as item discrimination parameters.
      \item Added omegaSem and omegaFromSem to calculate omega from a confirmatory factor model using John Fox's sem package.
      \item Added (with permission) Bond's Logical Operations Test (BLOT) data set.
    }
  }
  \subsection{Bugs Fixed}{
    \itemize{
      \item sim.irt was not properly handling vectors of item discrimination or guessing parameters.  Fixed.
      \item Minor problem in fa.stats that would lead to an error message if chisquare was NaN. 
      \item Minor bug in reverse.code for case of 1 variable (why one would want to use 1 the function for one variable is unclear.)
      \item tetrachoric would not apply the correction for continuity when given a table input but would when given normal data.  Fixed.
    }
  }
}
\section{Changes in psych version 1.0.90 (2010-07-07)}{
  \subsection{Additions}{
    \itemize{
      \item Added a few parameters to fa.parallel to a) set the y label and b) to show or not show a legend. Useful for multipanel displays.
      \item Minor improvements to fa.diagram to allow for input of factor loadings and intercorrelations from other functions.
      \item Added sim.omega to simulate multifactorial hierarchical structures as a test (and demonstration) of omega.
      \item Minor changes to the parameter call list for sim.minor to be compatible with sim.omega.
    }
  }
  \subsection{Bugs Fixed}{
    \itemize{
      \item Fixed a problem in ICLUST where it would sometimes not stop at the right number of clusters for n.clus = 1
      \item Fixed problem in omega.diagram where it would sometimes use the incorrect line type for hierarchical (non-Schmid-Leiman) diagrams.
      \item Fixed bug in schmid where it would not rotate exactly orthogonal factors and find exactly omega = 0.
      \item Added a warning message that oblique rotation does not work, and does a cluster based rotation instead in the case of exactly orthogonal factors.  This is just a problem for artificial data sets.
      \item History (Both of these are probably due to a problem in GPArotation which has since been fixed.)
    }
  }
}
\section{Changes in psych version 1.0.89 (2010-06-21)}{
  \subsection{Additions}{
    \itemize{
      \item Rearranged the order of parameters for fa so that it is easier to call.
      \item Added various sem stats to fa output (e.g. RMSEA, BIC, rms, crms).  These are also provided for omega output.
      \item Add fit statistics to omega for g solution only.  
      \item Added a measure of general factor adequacy in Omega.  This is just the percent of common variance that is general factor variance for each item. High values and low variances suggest a more meaningful estimate than lower average values or higher variances.
      \item Various clean up of output of ICLUST and iclust.diagram.  Added color option for positive and negative loadings.
      \item Changed ICC to work on complete data (i.e., using na.omit) as the default (following a question) by Ross Culloch
      \item Changed pairs.panels to allow for color choice of histogram (requested by Elaine Kuo) 
      \item Changed describe to handle numeric but annotated vectors (suggested by Soeren Vogel)
      \item Clarified documentation of cohen.kappa and allowed for  using explicitly non-numerical categories  (as suggested by Peter Ehlers  and Scot McNary)
      \item Increased the usefulness of error.bars to allow for input of x values, y values, and confidence interval. Also added the ability to use 1 standard deviation.   
      \item Minor change to print.psych.fa to round before using formatting. 
      \item Removed the cut <- min(cut,max(abs(load))/2)  following a request by  Reinhold Hatzinger to be consisent with factanal
      \item Added an option to omega to allow for analysis of data sets from just the factor loadings and correlation matrices. 
      \item Add a scree function to replace VSS.scree
    }
  }
  \subsection{Bugs Fixed}{
    \itemize{
      \item Fixed a bug in pairs.panels for the case of not smoothing but with data ellipses (reported by Hirayuki Sato)
      \item Fixed a serious bug in omega for the case of non-positive manifold items (e.g., personality items). The sign of the group factors was not correctly being set which would lead to incorrect estimates of the multiple R2.
      \item Fixed bug in ICLUST.diagram  where it would draw on top of other figures in rare cases.
      \item Fix bug in ICLUST for the path coefficients -- in the case of negative loadings, some of these were much too small.
      \item Correctly color the paths in omega.diagram. 
      \item Fixed print.psych.fa to correctly showing item numbers for sorted output. 
      \item Corrected tetrachor so it would not inappropriately claim that some cells had no data (this did not affect the correlations, just gave an inappropriate warning).
    }
  }
}
\section{Changes in psych version 1.0.88 (2010-04-24)}{
  \subsection{Additions}{
    \itemize{
      \item Added within option for error.bars and error.bars.by to allow for within subject error estimates.
      \item Reformatted this file (NEWS) so that the command  (e.g.,) news(query=Version > "1.0.85",package="psych")  provides nice information.
      \item Added more detail to the overview vignette to discuss how to find omega.
      \item Added the block.random function to block randomize for experimental designs.
    }
  }
  \subsection{Bugs Fixed}{
    \itemize{
      \item Corrected schmid (and thus, omega) so that it finds orthogonal factors before rotating them obliquely.  This corrected an error introduced in version 1.0.86 when the default rotation for fa was switched to be oblimin. 
      \item correct=FALSE option in tetrachor actually works
    }
  }
}
\section{Changes in psych version 1.0.87 (2010-04-04)}{
  \subsection{Bugs fixed}{
    \itemize{
      \item error.bars (and error.bars.by) would not properly draw < 3 variables
    }
  }
}
\section{Changes in psych version 1.0.86 (2010-03-26)}{
  \subsection{Additions}{
    \itemize{
      \item Added logit, logistic and logistic.grm functions to help teach about IRT.
      \item Added the ability to factor analyze covariance matrices (with resulting loadings of covariances rather than correlations) (suggested by Andreas Moeltner)
      \item Added two new functions to estimate the greatest lower bound (glb.algebraic, written by Andreas Moeltner) and glb.fa.  
      \item Added read.clipboard.fwf to allow reading of fixed width files from clipboard. (Originally this required adding a patch to the read.fwf function but that was subsequentally fixed.)
      \item Minor adjustment to ICLUST.cluster so that it will produce fewer clusters (if desired) than it thinks are appropriate.  That is, if n.clus is less than the appropriate solution, a warning is issued, but n.clus clusters are found.
      \item Changed the default rotation in fa to be oblimin.  (Subsequently discovered that this introduced a bug into schmid and therefore omega.)
      \item Added the glb.communality function to calculate the greatest lower bound by using factor communality estimates (based upon correspondece with Andreas Moeltner, although his algorithm (glb.algebraic, which is now included) which uses the Rcsdp provides slightly different estimates).
      \item added names to communalities in principal (requested by Gumundur Arnkelsson).
      \item Minor change to headtail
      \item Added the sim.rasch function to simulate 1 parameter IRT Rasch models.
      \item Added the data sets bock and income. bock is the source of the lsat data,income is US family income from 2008.
      \item Seriously modified the estimation of Cohen's Kappa in the wkappa function and added a new function: cohen.kappa to calculate kappa from raw data or from similarity matrices.  Also calculates the variances following Fleiss, Cohen, and Blashfield, and then find the confidence limits.
      \item Modified the sim.minor function to simulate correlation matrices with nfact major factors and nvar/2 minor factors with larger number of variables.   
    }
  }
  \subsection{Bugs fixed}{
    \itemize{
      \item Fixed print.psych.fa to give correct communalities for sorted factors (had been not sorting the communalities, although it had the factor loadings).
      \item Minor change to comorbidity to clean up output (using print.psych)
      \item Fixed guttman to correctly estimate lambda 4 (error pointed out by Andreas Moeltner). 
    }
  }
}
\section{Changes in psych version 1.0.85 (2009-12-20)}{
  \subsection{Additions}{
    \itemize{
      \item Modified the output from fa.parallel to print summary statistics.  Also added an option to use smcs as communality estimates.
      \item Added the sim.minor function to generate major and minor factor structures. 
      \item Modified fa to use first derivatives based upon the model rather than just empirically derived.  This results in an impressive speed improvement.
      \item Continued to improve the dia.x functions.  Labeled the entire set of functions diagram (for easier search).
    }
  }
}
\section{Changes in psych version 1.0.84 (2009-11-30)}{
  \subsection{Additions}{
    \itemize{
      \item Improved documentation for the dia functions to show the power of structure diagrams
      \item Improved structure.diagram to allow graphs to go left to right (default) or bottom to top (which means causal arrows go left to right)
      \item Added the affect data set which includes data from two experimental studies of affect.
      \item Added some features from pairs to allow group coloring in pairs.panels
      \item Modified factor.parallel to allow for other factoring models (including minres)
      \item Modified fa to give eigenvalues of the correlation matrix, and of the common factor space.
      \item Modified pairs.panels so it can plot even if the correlations are NA
      \item Added a legend and improved the output of cor.plot.
      \item Having learned how to use strwidth, have cleaned up the dia.rect, fa.diagram, structure.diagram, and iclust.diagram functions.
      \item Changed the way fa.parallel handles correlation matrices if the number of subjects is not specified. Rather than issuing a warning and proceeding as if n.obs = nvars, it now assumes n.obs =100.
      \item Improved documentation for fa to explain minimum correlation of factor score estimates. 
    }
  }
  \subsection{Bugs fixed}{
    \itemize{
      \item Fixed a problem with printing of the principal (components) results
      \item Fixed bug in the print.psych.vss function that returned the wrong value (although identifying the correct number) for the  Velicer MAP criterion. (Reported by Maike Luhmann).
    }
  }
}
\section{Changes in psych version 1.0.83 (2009-10-26)}{
  \subsection{Additions}{
    \itemize{
      \item Added Suggests Rgraphviz to DESCRIPTION file to get around CRAN test for 2.10 (as I should have done according to documentation for packages).  Removed all tests of functions which require RGraphviz for 2.9.2 because Rgraphviz does not work on SnowLeopard. 
    }
  }
}
\section{Changes in psych version 1.0.82 (2009-10-26)}{
  \subsection{Additions}{
    \itemize{
      \item Improved the output for mat.regress to allow print() and summary().
      \item Corrected bug to fa for case of oblique and negatively correlated factors (reported by Erich Studerus).
      \item Corrected problem in describe.by to properly identify the matrix output 
    }
  }
  \subsection{Bugs fixed}{
    \itemize{
      \item Minor improvements to the dia.x functions 
      \item Fixed the formatting of the news file to follow the appropriate format
      \item Removed the old supplementary  manual from the docs folder to make for a cleaner installation
    }
  }
}
\section{Changes in psych version 1.0.81 (2009-10-04)}{
  \subsection{Additions}{
    \itemize{
      \item Improved the documentation for the dia.X functions and made some of the calls simpler.  
    }
  }
}
\section{Changes in psych version 1.0.80 (2009-09-30)}{
  \subsection{Additions}{
    \itemize{
      \item Added the ability to draw "diagrams" to replace the functions that required Rgraphviz.  These are path diagrams with rectangles for observed variables and ellipses for latent variables.  A set of graphic functions for drawing diagrams includes dia.ellipse, dia.rect, dia.arrow, dia.curve and dia.self.  This allows the  function omega and ICLUST to draw structural diagrams without using Rgraphviz.
      \item The documentation for these new functions is a bit less than desirable and the dia primitives will be modified slightly for the next release.
    }
  }
}
\section{Changes in psych version 1.0.79 (2009-08-25)}{
  \subsection{Additions}{
    \itemize{
      \item added the Schmid data set which includes the original Schmid-Leiman correlation matrix as well as data from Chen, West, and Sousa (2006).
      \item modified geometric.mean and harmonic.mean to allow for missing data.
      \item Cleaned up describe.by so that it can handle matrix output of multiple grouping variables and then print correctly.
    }
  }
}
\section{Changes in psych version 1.0.78 (2009-07-28)}{
  \subsection{Additions}{
    \itemize{
      \item Modified factor.stats to allow for statistics even if the matrix is very singular (i.e., N var >> n.obs).
    }
  }
  \subsection{Bugs fixed}{
    \itemize{
      \item Corrected bug to ICLUST introduced in version 74 that can, in some conditions, produce negative betas. 
    }
  }
}
\section{Changes in psych version 1.0.77}{
  \subsection{Bugs fixed}{
    \itemize{
      \item Fixed an error in printing communalities for rotated factor solutions reported by Tamaki Hattori. Also resolved a problem of extra lines appearing in error.bars.  Minor revision to polar. 
    }
  }
}
\section{Changes in psych version 1.0.76}{
  \subsection{Bugs fixed}{
    \itemize{
      \item Finally resolved the problem of not passing the CRAN tests for Mac (even though developed on a Mac, the test package assumes no other packages).
    }
  }
}
\section{Changes in psych version 1.0.75 (2009-07-01)}{
  \subsection{Additions}{
    \itemize{
      \item Reinstituted reporting communalities and uniquenesses for principal components and for factor analysis.  This feature had been dropped sometime in the past and reported missing by Ista Zahn.
      \item Added two new data sets: Harman.Holzinger and Harman.Burt.  The former is used by Peter Bentler as an example of alternative methods for estimating   reliability, the latter is a nice example of what happens if the correlation matrix is singular.  smc estimates are inappropriate in this case and need to be checked.
    }
  }
}
\section{Changes in psych version 1.0.74 (2009-06-25)}{
  \subsection{Additions}{
    \itemize{
      \item Added a weighting option to ICLUST.  Weighted beta finds the average between cluster correlation and takes the ratio of the  n ^2 * average between cluster correlation to the total variance of the pooled clusters as a measure of the general factor.  That is, Unweighted beta calculates beta based upon the correlation between  two clusters, corrected for test length using the Spearman-Brown prophecy formala, while weighted beta finds the average interitem correlation between the items  within two clusters and then finds beta from this.  That is, for two clusters A and B of size N and M with between average correlation rb, weighted beta is (N+M)^2 rb/(Va +Vb + 2Cab).  Raw (unweighted) beta is 2rab/(1+rab) where rab = Cab/sqrt(VaVb).   Weighted beta seems a more appropriate estimate and is now the default.  Unweighted beta is still available for consistency with prior versions. (Unfortunately, in doing this, a bug was introduced that could produce negative betas.  Fixed in Version 78).
      \item Modified the fa function to include maximimum likelihood factor analysis.
    }
  }
}
\section{Changes in psych version 1.0.73 (2009-06-08)}{
  \subsection{Additions}{
    \itemize{
      \item Added gls fit option to fa
      \item Added mat.sort to sort correlation matrices by factor loadings.  
      \item Revised the overview vignette to make more readable.
      \item Fixed (I hope) the documentation error that has been preventing being loaded onto CRAN for Mac.
    }
  }
}
\section{Changes in psych version 1.0.72 (2009-06-02)}{
  \subsection{Additions}{
    \itemize{
      \item Added factor stats to omega
    }
  }
  \subsection{Bugs fixed}{
    \itemize{
      \item Minor fix to print.psych.omega to correct printing error introduced in .71
    }
  }
}
\section{Changes in psych version 1.0.71 (2009-06-01)}{
  \subsection{Additions}{
    \itemize{
      \item Combined principal axis, unweighted least squares (minres) and weighted least squares into one combined function (fa). Eventually the three other functions (factor.pa,factor.wls, and factor.minres) will be phased out. Added more rotation and transformation options (basically by making calls to the GPArotation package).
      \item factor.stats now reports the factor weights matrix for factor scores.
      \item A new factor analysis function factor.wls will give weighted least squares (wls).
      \item Omega will now give an estimate for  1  factor and a message about how there are alternatives for 2 factors.
    }
  }
  \subsection{Bugs fixed}{
    \itemize{
      \item Fixed a bug in print.psych.omega so that it will now report the item numbers correctly when sorting.
      \item Fixed call in phi.demo so that it checks if the polycor package is installed.  This was causing problems in the Mac version on CRAN.
    }
  }
}
\section{Changes in psych version 1.0.70 (2009-05-25)}{
  \subsection{Additions}{
    \itemize{
      \item Added factor.stats function to report a variety of statistics for factors or components.  In particular, added measures of factorial adequacy.  factor.stats is now called by the factor.pa, factor.minres, and principal, so it is easier to adjust all three of these functions.
      \item Added a new rotation option to target.rot -- if keys are not specified, it will rotate to a cluster solution.
      \item minor tweak to sim.hierarchical to make it consistent with the Jensen and Weng paper
      \item added to omega to calculate omege_infinity
      \item corrected bug in principal for case of singular matrices (reported by Gudmundur Arnkelsson) 
      \item corrected bug in print.psych so that it now sorts by factor loadings correctly 
    }
  }
}
\section{Changes in psych version 1.0.69 (2009-04-25)}{
  \subsection{Additions}{
    \itemize{
      \item added cor.plot to graphically display correlation structures
      \item modified cosinor to do regression as well as optimization for fits
      \item minor correction to ICLUST to catch non-invertible matrices
    }
  }
}
\section{Changes in psych version 1.0.68 (2009-04-15)}{
  \subsection{Additions}{
    \itemize{
      \item Added method="minres" to do minimal residual (ULS) factor analysis to the factor.pa function. 
      \item Added the target.rot function to do targeted rotation.
      \item Added the sim.anova function to simulate 3 way ANOVA or linear models.
      \item Minor cleanups of code in other functions
    }
  }
}
\section{Changes in psych version 1.0.67 (2009-03-15)}{
  \subsection{Additions}{
    \itemize{
      \item Have included the vignettes in the package
      \item Serious modifications to mat.regress to allow for raw data and covariances to be used.
      \item Added corrections for item overlap to use estimated communalities of items to alpha, cluster.cor, etc.
      \item Added estimates of G6+  (Modified Guttman 6) reliabililty to alpha, score.items
    }
  }
}
\section{Changes in psych version 1.0.66 (2009-02-15)}{
  \subsection{Additions}{
    \itemize{
      \item Introduced two package vignettes (overview and psych_for_sem)
      \item minor modification to ICLUST with respect to labeling
      \item addition to factor.congruence to allow for comparing omega solutions with factor solutions
      \item modified describe to allow for descriptions of categorical variables.
    }
  }
}
\section{Changes in psych version 1.0.64 (2009-02-01)}{
  \subsection{Additions}{
    \itemize{
      \item Added the function alpha to find Cronbach's alpha and Guttman's Lambda 6 for total scales as well as with single item deleted.
      \item Minor enhancements to the structure.graph and structure.sem functions. Further refinements to the psych_for_sem.pdf tutorial.  
    }
  }
}
\section{Changes in psych version 1.0.63 (2009-01-15)}{
  \subsection{Additions}{
    \itemize{
      \item Added structure.graph and structure.sem to create graphical output and pre-sem commands for sem.
      \item Modified omega.graph to return sem commands for running the sem function from John Fox.
      \item Added Promax to report factor intercorrelations for promax rotations.
      \item bug fix in VSS and omega to allow for promax rotations.
      \item correction to winsor so that it does not sort the vectors of the input data.frame.
      \item bug fix in principal (for problem of oblique rotations)
      \item continued improvements in the print.psych and summary.psych functions 
      \item cleaned up sim.congeneric to reduce the amount of output (now handled through print.psych)
      \item added super.matrix function (to combine two matrices -- useful for constructing keys and for doing examples in structure.graph)
      \item Improvements to structure graph.
    }
  }
}
\section{Changes in psych version 1.0.62 (2008-12-31)}{
  \subsection{Additions}{
    \itemize{
      \item added the structure.graph function
      \item added circadian functions 
      \item renamed some functions to make the index order easier to follow.
    }
  }
}
\section{Changes in psych version 1.0.58 (2009-09-15)}{
  \subsection{Additions}{
    \itemize{
      \item Introduced the guttman function to estimate the 6 Guttman (1945) estimates (as well as a few others).
    }
  }
}
\section{Changes in psych version 1.0.54 (2009-08-15)}{
  \subsection{Additions}{
    \itemize{
      \item Added print.psych and summary.psych generic functions to clean up the output of many of the functions. 
    }
  }
}
\section{Changes in psych version 1.0.50 (2007-11-15)}{
  \subsection{Additions}{
    \itemize{
      \item Introduced several new functions:
      \item headtail   (combines head and tail functions to produce the first N and last M rows of a matrix or dataframe)
      \item error.bars and 
      \item error.bars.by     plot means with error bars by a grouping variable
      \item omega revised to allow for negatively keyed items.  Detects items that should be reversed.
    }
  }
}

\section{Changes in psych version 1.0.42 (2007-03-24)}{
  \subsection{Additions}{
    \itemize{
      \item A relatively stable release. 
    }
  }
}

\section{Changes in psych version 1.0.40 (2007-08-15)}{
  \subsection{Additions}{
    \itemize{
      \item score.multiple.choice   provides scores and item statistics for multiple choice items.
      \item circ.sim modified to be more general item simulations
      \item circ.simulation
      \item item.sim
      \item polar converts two dimensional factor loadings to polar coordinates.
    }
  }
}

\section{Changes in psych version 1.0.17 (2007-05-06)}{
  \subsection{Additions}{
    \itemize{
      \item The first release to CRAN
    }
  }
}
