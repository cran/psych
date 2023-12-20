\name{corTest}
\alias{corTest}
\alias{corr.test}
\alias{corr.p}
\title{Find the correlations, sample sizes, and probability values between elements of a matrix or data.frame.   }
\description{Although the cor function finds the correlations for a matrix,  it does not report probability values.  cor.test does, but for only one pair of variables at a time.  corr.test uses cor to find the correlations for either complete or pairwise data and reports the sample sizes and probability values as well. For symmetric matrices, raw probabilites are reported below the diagonal and correlations adjusted for multiple comparisons above the diagonal. In the case of different x and ys, the default is to adjust the probabilities for multiple tests. Both corr.test and corr.p return raw and adjusted confidence intervals for each correlation. 
}
\usage{
corTest(x, y = NULL, use = "pairwise",method="pearson",adjust="holm", 
    alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corr.test(x, y = NULL, use = "pairwise",method="pearson",adjust="holm", 
    alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corr.p(r,n,adjust="holm",alpha=.05,minlength=5,ci=TRUE)
}

\arguments{
  \item{x}{A matrix or dataframe }
  \item{y}{A second matrix or dataframe with the same number of rows as x }
  \item{use}{use="pairwise" is the default value and will do pairwise deletion of cases. use="complete" will select just complete cases. }
  \item{method}{method="pearson" is the default value.  The alternatives to be passed to cor are "spearman" and "kendall".  These last two are much slower, particularly for big data sets. }
  \item{adjust}{What adjustment for multiple tests should be used? ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). See \code{\link{p.adjust}} for details about why to use "holm" rather than "bonferroni"). }
  \item{alpha}{alpha level of confidence intervals}
  \item{r}{A correlation matrix}
  \item{n}{Number of observations if using corr.p. May be either a matrix (as returned from corr.test, or a scaler. Set to n - np if finding the significance of partial correlations. (See below). }
  \item{ci}{By default, confidence intervals are found.  However, this leads to a noticable slowdown of speed, particularly for large problems.  So, for just the rs, ts and ps, set ci=FALSE}
  \item{minlength}{What is the minimum length for abbreviations.  Defaults to 5.}
  \item{normal}{By default, probabilities for method="spearman" and method="kendall" are found by normal theory.  If normal=="FALSE", then repetitive calls are made to cor.test.  This is much slower, but gives more accurate p values. exact is set to be FALSE which means that exact p values for small samples are not found given the problem of ties. }
  
  }
\details{corr.test uses the \code{\link{cor}} function to find the correlations, and then applies a t-test to the individual correlations using the formula
\deqn{t = \frac{r * \sqrt(n-2)}{\sqrt(1-r^2)}
}{t = r* sqrt(n-2)/sqrt(1-r^2) }

\deqn{se = \sqrt(\frac{1-r^2}{n-2}) }{se = sqrt((1-r^2)/(n-2))}

The t and Standard Errors are returned as objects in the result, but are not normally displayed. Confidence intervals are found and printed if using the print(short=FALSE) option.  These are found by using the fisher z transform of the correlation and then taking the range r +/- qnorm(alpha/2) *  se and the standard error of the z transforms is \deqn{se = \sqrt(\frac {1}{n-3}) }{se = sqrt(1/(n-3))}.   These values are then back transformed to be in correlation units. They are returned in the ci object. 

Note that in the case of method=="kendall" since these are the normal theory confidence intervals they are slightly too big.

The probability values may be adjusted using the Holm (or other) correction.  If the matrix is symmetric (no y data), then the original p values are reported below the diagonal and the adjusted above the diagonal.  Otherwise, all probabilities are adjusted (unless adjust="none").  This is made explicit in the output. Confidence intervals are shown for raw and adjusted probabilities in the ci object.

For those who like the conventional use of "magic asterisks" to show (stars) to represent conventional levels of significance, the object stars is returned (but not shown)).  See the examples.

\code{\link{corr.p}} may be applied to the results of \code{\link{partial.r}} if n is set to n - s (where s is the number of variables partialed out)  Fisher, 1924. 
}
\value{
  \item{r}{The matrix of correlations}
  \item{n}{Number of cases per correlation}
  \item{t}{value of t-test for each correlation}
  \item{p}{two tailed probability of t for each correlation.  For symmetric matrices, p values adjusted for multiple tests are reported above the diagonal. }
  \item{se}{standard error of the correlation}
  \item{ci}{the alpha/2 lower and upper values.}
  \item{ci2}{ci but with the adjusted pvalues as well.  This was added after tests showed we were breaking some packages that were calling the ci object without bothering to check for its dimensions.}
  \item{ci.adj}{These are the adjusted ((Holm or Bonferroni) confidence intervals. If asking to not adjust, the Holm adjustments for the confidence intervals are shown anyway, but the probability values are not adjusted and the appropriate confidence intervals are shown in the ci object.  }
  \item{stars}{For those people who like to show magic asterisks denoting ``statistical significance" the stars object flags those correlation values that are unlikely given normal theory. See the last example for  how to print these neatly.}
}

\note{For very large matrices (> 200 x 200), there is a noticeable speed improvement if confidence intervals are not found.

That adjusted confidence intervals are shown even when asking for no adjustment might be confusing.  If you don't want adjusted intervals, just use the ci object.  The adjusted values are given in the ci.adj object. }

\seealso{  \code{\link{cor.test}} for tests of a single correlation,  Hmisc::rcorr for an equivalant function, \code{\link{r.test}} to test the difference between correlations, and \code{\link{cortest.mat}} to test for equality of two correlation matrices. 

Also see \code{\link{cor.ci}} for bootstrapped confidence intervals of Pearson, Spearman, Kendall, tetrachoric or polychoric correlations.  In addition \code{\link{cor.ci}} will find bootstrapped estimates of composite scales based upon a set of correlations (ala \code{\link{cluster.cor}}). 

In particular, see \code{\link{p.adjust}} for a discussion of p values associated with multiple tests.

Other useful functions related to finding and displaying correlations include \code{link{corPlot}} to graphically display the correlation matrix,  and \code{\link{lowerCor}} for finding the correlations and then displaying the lower off diagonal using the \code{\link{lowerMat}} function. \code{\link{lowerUpper}} to compare two correlation matrices.  Also see \code{\link{pairs.panels}} to show the correlations and scatter plots.  


}
\examples{
ct  <- corTest(attitude)
#ct <- corr.test(attitude)  #find the correlations and give the probabilities
ct #show the results


cts <- corr.test(attitude[1:3],attitude[4:6]) #reports all values corrected for multiple tests

#corr.test(sat.act[1:3],sat.act[4:6],adjust="none")  #don't adjust the probabilities

#take correlations and show the probabilities as well as the confidence intervals
print(corr.p(cts$r,n=30),short=FALSE)  

#don't adjust the probabilities
print(corr.test(sat.act[1:3],sat.act[4:6],adjust="none"),short=FALSE)  

#print out the stars object without showing quotes
print(corr.test(attitude)$stars,quote=FALSE)  #note that the adjusted ps are given as well


kendall.r <- corr.test(bfi[1:40,4:6], method="kendall", normal=FALSE)
#compare with 
cor.test(x=bfi[1:40,4],y=bfi[1:40,6],method="kendall", exact=FALSE)
print(kendall.r,digits=6)
}

\keyword{multivariate }
\keyword{ models }
