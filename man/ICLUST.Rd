\name{ICLUST}
\alias{ICLUST}
\alias{iclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ICLUST: Item Cluster Analysis   -- Hierarchical cluster analysis using psychometric principles }
\description{A common data reduction technique is to cluster cases (subjects). Less common, but particularly useful in psychological research, is to cluster items (variables). This may be thought of as an alternative to factor analysis, based upon a much simpler model. The cluster model is that the correlations between variables reflect that each item loads on at most one cluster, and that items that load on those clusters correlate as a function of their respective loadings on that cluster and items that define different clusters correlate as a function of their respective cluster loadings and the intercluster correlations. 
Essentially, the cluster model is a Very Simple Structure factor model of complexity one (see \code{\link{VSS}}).

This function applies the ICLUST algorithm to hierarchically cluster items to form composite scales. Clusters are combined if coefficients alpha and beta will increase in the new cluster.

Alpha, the mean split half correlation, and beta, the worst split half correlation, are estimates of the reliability and general factor saturation of the test.  (See also the \code{\link{omega}} function to estimate McDonald's coeffients   \eqn{\omega_h}{omega hierarchical} and   \eqn{\omega_t}{omega total})


}
\usage{
ICLUST(r.mat, nclusters=0, alpha=3, beta=1, beta.size=4, alpha.size=3,
correct=TRUE,correct.cluster=TRUE, reverse=TRUE, beta.min=.5, output=1, digits=2,labels=NULL,cut=0,
n.iterations = 0,title="ICLUST",plot=TRUE)

#ICLUST(r.mat)    #use all defaults
#ICLUST(r.mat,nclusters =3)    #use all defaults and if possible stop at 3 clusters
#ICLUST(r.mat, output =3)     #long output shows clustering history
#ICLUST(r.mat, n.iterations =3)  #clean up solution by item reassignment
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.mat}{ A correlation matrix or data matrix/data.frame. (If r.mat is not square i.e, a correlation matrix, the data are correlated using pairwise deletion. }
  \item{nclusters}{Extract clusters until nclusters remain (default will extract until the other criteria are met or 1 cluster, whichever happens first). See the discussion below for alternative techniques for specifying the number of clusters. }
  \item{alpha}{ Apply the increase in alpha criterion  (0) never or for (1) the smaller, 2) the average, or 3) the greater of the separate alphas. (default = 3) }
  \item{beta}{  Apply the increase in beta criterion (0) never or for (1) the smaller, 2) the average, or 3) the greater of the separate betas. (default =1) }
  \item{beta.size}{ Apply the beta criterion after clusters are of beta.size  (default = 4)}
  \item{alpha.size}{ Apply the alpha criterion after clusters are of size alpha.size (default =3) }
  \item{correct}{ Correct correlations for reliability (default = TRUE) }
  \item{correct.cluster}{Correct cluster -sub cluster correlations for reliability of the sub cluster , default is TRUE))} 
  \item{reverse}{Reverse negative keyed items (default = TRUE}
  \item{beta.min}{ Stop clustering if the beta is not greater than beta.min (default = .5) }
  \item{output}{ 1) short, 2) medium, 3 ) long output (default =1)}
  \item{labels}{vector of item content or labels}
  \item{cut}{sort cluster loadings > absolute(cut) (default = 0) }
  \item{n.iterations}{iterate the solution n.iterations times to "purify" the clusters (default = 0)}
  \item{digits}{ Precision of digits of output (default = 2) }
  \item{title}{ Title for this run }
  \item{plot}{Should ICLUST.rgraph be called automatically for plotting (requires Rgraphviz default=TRUE)}
}
\details{ Extensive documentation and justification of the algorithm is available in the original MBR 1979 \url{http://personality-project.org/revelle/publications/iclust.pdf} paper.  Further discussion of the algorithm and sample output is available on the personality-project.org web page: \url{http://personality-project.org/r/r.ICLUST.html} 

The results are best visualized using  \code{\link{ICLUST.graph}}, the results of which can be saved as a dot file for the Graphviz program.  \url{http://www.graphviz.org/}.  With the installation of Rgraphviz, ICLUST will automatically provide cluster graphs.

 A common problem in the social sciences is to construct scales or composites of items to measure constructs of theoretical interest and practical importance. This process frequently involves administering a battery of items from which those that meet certain criteria are selected. These criteria might be rational, empirical,or factorial. A similar problem is to analyze the adequacy of scales that already have been formed and to decide whether the putative constructs are measured properly. Both of these problems have been discussed in numerous texts,  as well as in myriad articles. Proponents of various methods have argued for the importance of face validity, discriminant validity, construct validity, factorial homogeneity, and theoretical importance. 
 
Revelle (1979) proposed that hierachical cluster analysis could be used to estimate a new coefficient (beta) that was an estimate of the  general factor saturation of a test.  More recently, Zinbarg, Revelle, Yovel and Li (2005) compared McDonald's Omega to Chronbach's alpha and Revelle's beta. They conclude that   \eqn{\omega_h}{omega} hierarchical is the best estimate.  An algorithm for estimating \code{\link{omega} } is available as part of this package. 

Revelle and Zinbarg (2009) discuss alpha, beta, and omega, as well as other estimates of reliability. 


ICLUST was completely rewritten for the psych package.  Please email me if you want help with this version of ICLUST or if you desire more features.

A requested feature (not yet available) is to specify certain items as forming a cluster.  That is, to do confirmatory cluster analysis.  

The program currently has three primary functions: cluster, loadings, and graphics.  

Although ICLUST will give what it thinks is the best solution in terms of the number of clusters to extract, the user will sometimes disagree.  To get more clusters than the default solution, just set the nclusters parameter to the number desired.  However, to get fewer than meet the alpha and beta criteria, it is sometimes necessary to set alpha=0 and beta=0 and then set the nclusters to the desired number.   

Clustering 24 tests of mental ability

 A sample output using the 24 variable problem by Harman can be represented both graphically and in terms of the cluster order.  Note that the graphic is created using GraphViz in the dot language.  \code{\link{ICLUST.graph}} produces the dot code for Graphviz.  Somewhat lower resolution graphs with fewer options are available in the \code{\link{ICLUST.rgraph}} function which requires Rgraphviz.  Dot code can be viewed directly in Graphviz or can be tweaked using commercial software packages (e.g.,OmniGraffle)

Note that for this problem, with these parameters, the data form one large cluster. (This is consistent with the Very Simple Structure (\code{\link{VSS}}) output as well, which shows a clear one factor solution for complexity 1 data.)  See below for an example with this same data set, but with more stringent parameter settings.  

At least for the Harman 24 mental ability measures, it is interesting to compare the cluster pattern matrix with the oblique rotation solution from a factor analysis.  The factor congruence of a four factor oblique pattern solution with the four cluster solution is > .99 for three of the four clusters and > .97 for the fourth cluster.  

To see the graphic output go to \url{http://personality-project.org/r/r.ICLUST.html} or use \code{\link{ICLUST.rgraph}} (requires Rgraphviz). 

}
\value{
 \item{title }{Name of this run}
 \item{results}{A list containing the step by step cluster history, including which pair was grouped, what were the alpha and betas of the two groups and of the combined group.  
  
  Note that the alpha values are ``standardized alphas'' based upon the correlation matrix, rather than the raw alphas that will come from \code{\link{score.items}}
  
  The print.psych and summary.psych functions will print out just the must important results.}
 \item{corrected}{The raw and corrected for alpha reliability cluster intercorrelations.}
 \item{clusters}{a matrix of -1,0, and 1 values to define cluster membership.}
 \item{purified}{A list of the cluster definitions and cluster loadings of the purified solution.
   
   To show just the most salient items, use the cutoff option in \code{\link{print.psych}} }
 \item{cluster.fit, structure.fit, pattern.fit}{There are a number of ways to evaluate how well any factor or cluster matrix reproduces the original matrix. Cluster fit considers how well the clusters fit if only correlations with clusters are considered.  Structure fit evaluates R = CC' while pattern fit evaluate R = C  inverse (phi) C'  where C is the cluster loading matrix, and phi is the intercluster correlation matrix.}
}

\references{Revelle, W. Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 1979, 14, 57-74.

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma. Psychometrika, 2009. 

\url{http://personality-project.org/revelle/publications/iclust.pdf} \cr
See also  more extensive documentation at 
\url{http://personality-project.org/r/r.ICLUST.html} and \cr
Revelle, W. (in prep) An introduction to psychometric theory with applications in R. To be published by Springer.  (working draft available at  \url{http://personality-project.org/r/book/} }
\author{William Revelle   }
        

\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{omega}} }

\examples{
test.data <- Harman74.cor$cov
ic.out <- ICLUST(test.data)
summary(ic.out)
ic.out <- ICLUST(test.data,nclusters =4)  #use all defaults and stop at 4 clusters
print(ic.out)
plot(ic.out)    #this shows the spatial representation



}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
