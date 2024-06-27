\name{bigCor}
\alias{bigCor}
\title{Find large correlation matrices by stitching together smaller ones found more rapidly}
\description{When analyzing many subjects (ie. 100,000 or more) with many variables (i.e. 1000 or more) core R can take a long time and sometime exceed  memory limits (i.e. with 600K subjects and 6K variables).  bigCor runs (in parallel if multicores are available) by breaking the variables into subsets (of size=size), finding all subset correlations, and then stitches the resulting matrices into one large matrix.   Noticeable improvements in speed compared to cor.  
}
\usage{
bigCor(x, size = NULL, use = "pairwise",cor="pearson",correct=.5)
}

\arguments{
  \item{x}{A data set of numeric variables}
  \item{size}{What should the size of the subsets be? Defaults to NCOL (x)/20  
}
  \item{use}{The standard correlation option.  "pairwise" allows for missing data
}
\item{cor}{Defaults to Pearson correlations, alteratives are polychoric and spearman  }
\item{correct}{Correction for continuity for polychoric correlations. (see \code{\link{polychoric}})}
}
\details{
The data are divided into subsets of size=size.  Correlations are then found for each subset and pairs of subsets. 

Time is roughly linear with the number of cases and increases by the square of the number of variables.  The benefit of more cores is noticeable. It seems as if with 4 cores, we should use sizes to split it into 8 or 12 sets.  Otherwise we don't actually use all cores efficiently. 

There is some overhead in using multicores.  So for smaller problems (e.g. the 4,000 cases of the  145 items of the psychTools::spi data set, the timings are roughly .14 seconds for bigCor (default size) and .10 for normal cor. For small problems, this actually  gets worse as we use more cores.    The cross over point seems to be at roughly 5K subjects.  (updated these timings to recognize the M1 Max chip.  An increase of 4x in speed!  They had been .44 and .36.)


The basic loop loops over the subsets. When the size is a integer subset of the number of variables and is a multiple of the number of cores, the multiple cores will be used more.  Notice the benefit of 660/80 versus 660/100.  But this breaks down if we try 660/165.  Further notice the benefit when using a smaller subset (55) which led to the 4 cores being used more.  

The following timings are included to help users tinker with parameters:



 Timings (in seconds) for various problems with 645K subjects on an 8 core Mac Book Pro with a 2.4 GHZ Intell core i9.

options(mc.cores=4) (Because we have 8 we can work at the same time as we test this.)

First test it with 644,495 subjects and 1/10 of the number of possible variables.  Then test it for somewhat fewer variables.  
\tabular{lllll}{
Variables  \tab size         \tab 2 cores     \tab 4 cores    \tab compared to normal cor function \cr


660  \tab 100  \tab  430   \tab 434 \tab  430 \cr

660  \tab 80   \tab  600 \tab   348  \tab  notice the improvement with 8ths  \cr

660  \tab 165  \tab \tab    666  \tab (Stitching seems to have been very slow) \cr

660  \tab 55   \tab \tab    303  \tab  Even better if we break it into 12ths! \cr
\cr

 500  \tab 100  \tab  \tab 332   \tab 322 secs  \cr
\cr

480 \tab 120    \tab 408 \tab  365   \tab 315    Better to change the size \cr 

480 \tab 60    \tab  358 \tab   206  \tab This leads to 8 splits \cr

\cr 
}
We also test it with fewer subjects.  Time is roughly linear with number of subjects.

\tabular{lllll}{
Variables  \tab size      \tab 2 cores     \tab 4 cores    \tab compared to normal cor function 
Further comparisons with fewer subjects (100K) \cr
\cr 
480  \tab  60  \tab 57 \tab    31 \tab 47 with normal cor.  Note the effect of n subjects!  \cr
200 \tab 50  \tab  19.9   \tab  13.6 \tab 27.13 \cr
100  \tab 25    \tab 4.6  \tab 3.5   \tab 5.85  \cr
 \cr 
} 
One last comparison, 10,000 subjects, showing the effect of getting the proper size value. You can tune on these smaller sets of subjects before trying large problems.

\tabular{lllll}{
Variables  \tab size       \tab 2 cores     \tab 4 cores    \tab compared to normal cor function 
\cr
480  \tab 120 \tab 5.2    \tab   5.1  \tab  4.51 \cr
480  \tab 60 \tab 2.9    \tab   2.88  \tab  4.51 \cr
480  \tab   30  \tab 2.65  \tab  2.691  \cr
480  \tab 20   \tab 2.73 \tab    2.77 \tab  \cr
480   \tab   10   \tab 2.82  \tab    2.97 \tab too many splits? \cr

200  \tab  50   \tab 2.18  \tab 1.39 \tab 2.47 for normal cor (1.44 with 8 cores 2.99 with 1 core) \cr
200  \tab  25   \tab 1.2  \tab 1.17 \tab 2.47 for normal cor \cr (1.16 with 8 cores,  1.17 with 1 core) \cr
100  \tab 25    \tab .64  \tab .52   \tab .56 \cr
} 


Timings updated in 2/23 using a MacBook Pro with M1 max chip
  10,000 subjects 953 variables  suggests that a very small size (e.g. 20) is probably optimal  
\tabular{lllllll}{
Variables  \tab size     \tab 2 cores     \tab 4 cores  \tab 8 cores   \tab compared to normal cor function 
\cr
953 \tab 20 \tab 7.92    \tab 4.55  \tab 2.88 \tab 11.04 \cr
953 \tab 30 \tab 7.98    \tab 4.88   \tab 3.15 \tab 11.04 \cr
953 \tab 40 \tab 8.22 \tab 5.14 \tab 3.63 \tab 11.16  \cr
953 \tab 60 \tab 8.51 \tab 5.59 \tab 3.93 \tab 11.16  \cr
953 \tab 80  \tab 8.31 \tab 5.59 \tab 4.14 \tab 11.16  \cr
953 \tab 120  \tab 8.33 \tab 6.22  \tab 4.75 \tab 11.16  \cr
}

 

}

\note{
Does not seem to work with data.tables
}

\value{
The correlation matrix
}
\references{
Examples of large data sets with massively missing data are taken from the SAPA project. e.g.,

 William Revelle, Elizabeth M. Dworak, and David M. Condon (2021) Exploring the persome: The power of the item in understanding personality structure. Personality and Individual Differences, 169, \doi{10.1016/j.paid.2020.109905}
 
 David Condon (2018)The SAPA Personality Inventory: an empirically-derived, hierarchically-organized self-report personality assessment model.  PsyArXiv /sc4p9/ \doi{10.31234/osf.io/sc4p9}



}
\author{William Revelle }



\seealso{\code{\link{pairwiseCountBig}} which will do the same, but find the count of observations per cell.}


\examples{
R <- bigCor(bfi,10)
#compare the results with 
r.bfi <- cor(bfi,use="pairwise")
all.equal(R,r.bfi)
}

\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{multivariate }% __ONLY ONE__ keyword per line
% Use only one keyword per line.
