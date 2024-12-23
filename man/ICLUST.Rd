\name{iclust}
\alias{ICLUST}
\alias{iclust}
\title{iclust: Item Cluster Analysis   -- Hierarchical cluster analysis using psychometric principles }
\description{A common data reduction technique is to cluster cases (subjects). Less common, but particularly useful in psychological research, is to cluster items (variables). This may be thought of as an alternative to factor analysis, based upon a much simpler and more intuitive model. The cluster model is that the correlations between variables reflect that each item loads on at most one cluster, and that items that load on those clusters correlate as a function of their respective loadings on that cluster and items that define different clusters correlate as a function of their respective cluster loadings and the intercluster correlations. 
Essentially, the cluster model is a Very Simple Structure factor model of complexity one (see \code{\link{VSS}}).

\code{\link{iclust}} function applies the iclust algorithm (Revelle, 1979) to hierarchically cluster items to form composite scales. Clusters are combined if coefficients alpha and beta will increase in the new cluster.

 \eqn{\alpha}{alpha}, the mean split half correlation, and  \eqn{\beta}{beta}, an estimate of the general factor saturation based upon the correlation between the  most distinct partions of the correlation matrix, are estimates of the reliability and general factor saturation of the test.  (See also the \code{\link{omega}} function to estimate McDonald's coeffients   \eqn{\omega_h}{omega hierarchical} and   \eqn{\omega_t}{omega total}). Other reliability estimates are found in the \code{\link{reliability}} and \code{\link{splitHalf}}.


}
\usage{
iclust(r.mat, nclusters=0, alpha=3, beta=1, 
	beta.size=4, alpha.size=3,
	correct=TRUE,correct.cluster=TRUE,
	reverse=TRUE, beta.min=.5, output=1,
	digits=2,labels= NULL,cut=0, 
	n.iterations =0, 
	title="ICLUST",
	cor="cor", use="pairwise",
	plot=TRUE, 
	weighted=TRUE, 
	cor.gen=TRUE,SMC=TRUE,purify=TRUE,diagonal=FALSE,
	n.obs = NA, reliability=FALSE)

ICLUST(r.mat, nclusters=0, alpha=3, beta=1, 
	beta.size=4, alpha.size=3,
	correct=TRUE,correct.cluster=TRUE, 
	reverse=TRUE, beta.min=.5, output=1, 
	digits=2,labels=NULL,cut=0,
	n.iterations = 0,
	title="ICLUST",
	cor="cor", use="pairwise", plot=TRUE,
	weighted=TRUE,
	cor.gen=TRUE,SMC=TRUE,purify=TRUE,diagonal=FALSE,
	n.obs=NA,reliability=FALSE)


#iclust(r.mat)    #use all defaults
#iclust(r.mat,nclusters =3)    #use all defaults and if possible stop at 3 clusters
#ICLUST(r.mat, output =3)     #long output shows clustering history
#ICLUST(r.mat, n.iterations =3)  #clean up solution by item reassignment
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.mat}{ A correlation matrix or data matrix/data.frame. (If r.mat is not square i.e, is not a correlation matrix, the data are correlated using the options specified in cor.  The default is to  use  pairwise deletion for Pearson correlations. Alternatives include tetrachoric or polychoric correlations.  }
  \item{nclusters}{Extract clusters until nclusters remain (default will extract until the other criteria are met or 1 cluster, whichever happens first). See the discussion below for alternative techniques for specifying the number of clusters. }
  \item{alpha}{ Apply the increase in alpha criterion  (0) never or for (1) the smaller, 2) the average, or 3) the greater of the separate alphas. (default = 3).   }
  \item{beta}{  Apply the increase in beta criterion (0) never or for (1) the smaller, 2) the average, or 3) the greater of the separate betas. (default =1).  By setting a larger nvalue, the application of the beta criterion becomes more stringent.  See the examples. }
  \item{beta.size}{ Apply the beta criterion after clusters are of beta.size  (default = 4).  Setting this to a smaller values is more strigent.  }
  \item{alpha.size}{ Apply the alpha criterion after clusters are of size alpha.size (default =3). }
  \item{correct}{ Correct correlations for reliability (default = TRUE). }
  \item{correct.cluster}{Correct cluster -sub cluster correlations for reliability of the sub cluster , default is TRUE))} 
  \item{reverse}{Reverse negative keyed items (default = TRUE).  This is the principal that a reversed score item makes sense.  Probably not appropriate if clustering people. }
  \item{beta.min}{ Stop clustering if the beta is not greater than beta.min (default = .5) }
  \item{output}{ 1) short, 2) medium, 3) long, 4) very long output. (default =1). To see the cluster statistics at each level, use 3. To see the within cluster correlations at each level, use 4.  (This produces a great deal of output, but is useful to understand the algorithm.) }
  \item{labels}{Vector of item content or labels. If NULL, then the colnames are used.  If FALSE, then labels are V1 .. Vn.  See the examples for how to use a dictionary.}
  \item{cut}{Sort cluster loadings > absolute(cut) (default = 0) }
  \item{n.iterations}{Iterate the solution n.iterations times to "purify" the clusters (default = 0)}
  \item{digits}{ Precision of digits of output (default = 2) }
  \item{title}{ Title for this run. }
  \item{cor}{What kind of correlation should be applied to raw data, defaults to "pearson" with pairwise complete, options include "tet", "poly", "mixed","spearman".}
  \item{use}{See options for use in the cor function.  Defaults to "pairwise". }
  \item{plot}{Should ICLUST.diagram be called automatically for plotting (does not require Rgraphviz default=TRUE)}
  \item{weighted}{Weight the intercluster correlation by the size of the two clusters (TRUE) or do not weight them (FALSE). A third option ("min"), added 6/11/24 is to estimate the general saturation by the minimum of the within and between interitem cluster correlations.}
  \item{cor.gen}{When correlating clusters with subclusters, base the correlations on the general factor (default) or general + group (cor.gen=FALSE)}
  \item{SMC}{When estimating cluster-item correlations, use the smcs as the estimate of an item communality (SMC=TRUE) or use the maximum correlation (SMC=FALSE).}  
  \item{purify}{Should clusters be defined as the original groupings (purify = FALSE) or by the items with the highest loadings on those original clusters? (purify = TRUE)  }
  \item{diagonal}{Should the diagonal be included in the fit statistics.  The default is not to include it.  Prior to 1.2.8, the diagonal was included.} 
  \item{n.obs}{Number of observations (if using a correlation matrix).  Specify if using a correlation matrix in order to some  fit statistics that depend upon sample size (e.g., chi square, RMSEA, etc.)}
  \item{reliability}{Report various alternative estimates of reliability and return the reliability object. This will lead to problems if some clusters are just items.  For well formed clusters, set reliability =TRUE or just use the reliability function on the keys object.  }
  
}
\details{ Extensive documentation and justification of the algorithm is available in the original MBR 1979 \url{https://personality-project.org/revelle/publications/iclust.pdf} paper.  Further discussion of the algorithm and sample output is available on the personality-project.org web page: \url{https://personality-project.org/r/r.ICLUST.html} 

 A common problem in the social sciences is to construct scales or composites of items to measure constructs of theoretical interest and practical importance. This process frequently involves administering a battery of items from which those that meet certain criteria are selected. These criteria might be rational, empirical,or factorial. A similar problem is to analyze the adequacy of scales that already have been formed and to decide whether the putative constructs are measured properly. Both of these problems have been discussed in numerous texts,  as well as in myriad articles. Proponents of various methods have argued for the importance of face validity, discriminant validity, construct validity, factorial homogeneity, and theoretical importance. 
 
Revelle (1979) proposed that hierachical cluster analysis could be used to estimate a new coefficient (beta) that was an estimate of the  general factor saturation of a test.  More recently, Zinbarg, Revelle, Yovel and Li (2005) compared McDonald's Omega to Cronbach's alpha and Revelle's beta. They conclude that   \eqn{\omega_h}{omega} hierarchical is the best estimate.  An algorithm for estimating \code{\link{omega} } is available as part of the psych package. 

Revelle and Zinbarg (2009) discuss alpha, beta, and omega, as well as other estimates of reliability.

Reise and Haviland (2024) give a  nice discussion of alpha and beta for applied scale construction.


The original ICLUST program was written in FORTRAN in the early 1970s to run on CDC and IBM mainframes and was then modified to run in PC-DOS.  The R version of iclust is a completely new version  written for the psych package.  

A requested feature (not yet available) is to specify certain items as forming a cluster.  That is, to do confirmatory cluster analysis.  

The program currently has three primary functions: cluster, loadings, and graphics.  
\eqn{\beta}{beta}.  Conceptually \eqn{\beta}{beta} is similar to the worst split half reliability  (see e.g. \code{\link{splitHalf}})  but it is better conceived of as an estimate of the general factor saturation of a test. 

Consider the matrix M composed of four submatrices of size n x n and m x m

 \tabular{lll}{

    \tab   Rx  \tab Rxy \cr
 M = \tab  Rxy  \tab Ry  \cr
 } 
 with average within and between item correlations of 
 
   \tabular{lll}{

    \tab   rxx  \tab rxy \cr
 av.r  = \tab rxy  \tab ryy  \cr
 } 
 
 Then \eqn{Rx = \Sigma rxx } with n^2 elements and  \eqn{Rxy = \Sigma rxy } with n * m elements.  The total variance of the M (Vm) matrix is just the sum of the four submatrices.  If M has been partitioned such that rxy is minimized (which is the goal of clustering), then the correlations in Rxy reflect just the general factor in this test.  Thus 
 
 \eqn{\beta = \frac{n * m * rxy}{Vm}}{beta =\frac{n * m * rxy}{Vm}}.

That is, the general variance of the test is assumed to be  what the two most unrelated parts share.    


An important element of \code{\link{iclust}} is the ability to estimate coefficient \eqn{\beta}{beta}. 
In June, 2009, the option of weighted versus unweighted beta was introduced.  Unweighted beta calculates beta based upon the correlation between  two clusters, corrected for test length using the Spearman-Brown prophecy formala, while weighted beta finds the average interitem correlation between the items  within two clusters and then finds beta from this.  That is, for two clusters X and Y of size N and M with between average correlation rxy, weighted beta is (N+M)^2 rxy/(Vx +y + 2Cxy).  Raw (unweighted) beta is 2Rxy/(1+Rxy) where Rxy = Cxy/sqrt(VxVy).   Weighted beta seems a more appropriate estimate and is now the default.  Unweighted beta is still available for consistency with prior versions.

Added in June, 2024 is yet a third option , the minimum of rxx, ryy, and rxy where rxy is the average between cluster correlation.  This is logically the best option, but is still  being tested.  weighted=TRUE is still the default.

Also modified in June, 2009 was the way of correcting for item overlap when calculating the cluster-subcluster correlations for the graphic output.  This does not affect the final cluster solution, but does produce slightly different path values.  In addition, there are two ways to solve for the cluster - subcluster correlation.

Given the covariance between two clusters, Cab with average rab = Cab/(N*M), and cluster variances Va and Vb with Va = N  + N*(N-1)*ra  then the  correlation of cluster A with the combined cluster AB  is either

a) (((N+M)^2)rab + Cab)/sqrt(Va*Va)   (option cor.gen=TRUE) or \cr
b) (Va - N + Nra + Cab)/sqrt(Vab*Va)    (option cor.gen=FALSE)


The default is to use cor.gen=TRUE.


Although iclust will give what it thinks is the best solution in terms of the number of clusters to extract, the user will sometimes disagree.  To get more clusters than the default solution, just set the nclusters parameter to the number desired.  However, to get fewer than meet the alpha and beta criteria, it is sometimes necessary to set alpha=0 and beta=0 and then set the nclusters to the desired number.  

Hierarchical clustering is not guaranteed to give the "best" solution. Thus, following the original clustering, items are correlated with all clusters and reassigned based upon their highest correlation. This is shown in item by cluster Structure matrix.  The first column (O) is the original cluster definition, the second column (P) is the purified cluster solution.  This is seen in the example clustering of \code{\link{bfi}} data set.
The figure shows the original clustering and matches the O column of the Cluster Structure.

Thus, there are two sets of scoring keys reported: the purified keys, keys, and the original, unpurified keys, keys.org.  The reliability object reported is based upon the purified keys.  The alphas reported in the iclust.diagram are based upon the original keys.  To find the other statistics reported in the reliability object, run the \code{\link{reliability}} function on the keys.org  object.  (see the bfi cluster example below.)
  


Clustering 24 tests of mental ability

 A sample output using the 24 variable problem by Harman can be represented both graphically and in terms of the cluster order. The default is to produce graphics using the \code{\link{diagram}} functions.  An alternative is to use the Rgraphviz package (from BioConductor).  Because this package is sometimes hard to install, there is an alternative option (\code{\link{ICLUST.graph}} to write  dot language instructions for subsequent processing.  This will create a  graphic instructions suitable for any viewing program that uses the dot language.  \code{\link{ICLUST.rgraph}} produces the dot code for Graphviz.  Somewhat lower resolution graphs with fewer options are available in the \code{\link{ICLUST.diagram}} function which does not require Rgraphviz.  Dot code can be viewed directly in Graphviz or can be tweaked using commercial software packages (e.g., OmniGraffle)
 

Note that for the Harman 24 variable problem, with the default parameters, the data form one large cluster. (This is consistent with the Very Simple Structure (\code{\link{VSS}}) output as well, which shows a clear one factor solution for complexity 1 data.)  

An alternative solution is to ask for a somewhat more stringent set of criteria and require an increase in the size of beta for all clusters greater than 3 variables.  This produces a 4 cluster solution.

It is also possible to use the original parameter settings, but ask for a 4 cluster solution.

At least for the Harman 24 mental ability measures, it is interesting to compare the cluster pattern matrix with the oblique rotation solution from a factor analysis.  The factor congruence of a four factor oblique pattern solution with the four cluster solution is > .99 for three of the four clusters and > .97 for the fourth cluster.  The cluster pattern matrix is returned as an invisible object in the output. 

In September, 2012, the fit statistics (pattern fit and cluster fit) were slightly modified to (by default) not consider the diagonal (diagonal=FALSE).  Until then, the diagonal was included in the cluster fit statistics.  The pattern fit is analogous to factor analysis and is based upon the model = P x Structure where Structure is Pattern * Phi.  Then R* = R - model and fit is the ratio of sum(r*^2)/sum(r^2) for the off diagonal elements.  

The results are best visualized using  \code{\link{ICLUST.graph}}, the results of which can be saved as a dot file for the Graphviz program.  https://www.graphviz.org/. The \code{\link{iclust.diagram}} is called automatically to produce cluster diagrams.  The resulting diagram is not quite as pretty as what can be achieved in dot code but is quite adequate if you don't want to use an external graphics program. With the installation of Rgraphviz, ICLUST can also provide cluster graphs.

As of May, 2024, it is much easier to take the cluster results and if using raw data to find cluster scores.  ICLUST now returns two keys lists which can be passed to various scoring functions. See the discussion above about raw keys and purified keys.

As of June, 2024, following a discussion with Steven Reise, I have added a call to the \code{\link{reliability}} function which is passed the keys.list for each cluster.  This allows for alternative estimates of reliability for the clusters.  
}
\value{
 \item{title }{Name of this analysis}
 \item{results}{A list containing the step by step cluster history, including which pair was grouped, what were the alpha and betas of the two groups and of the combined group.  
  
  Note that the alpha values are ``standardized alphas'' based upon the correlation matrix, rather than the raw alphas that will come from \code{\link{scoreItems}}
  
  The print.psych and summary.psych functions will print out just the must important results.}
 \item{corrected}{The raw and corrected for alpha reliability cluster intercorrelations.}
 \item{clusters}{a matrix of -1, 0, and 1 values to define cluster membership.}
 \item{purified}{A list of the cluster definitions and cluster loadings of the purified solution.  These are sorted by importance (the eigenvalues of the clusters). The cluster membership from the original (O) and purified (P) clusters are indicated along with the cluster structure matrix.  These item loadings are the same as those found by the \code{\link{scoreItems}} function and are found by correcting the item-cluster correlation for item overlap by summing the item-cluster covariances with all except that item and then adding in the smc for that item. These resulting correlations are then corrected for scale reliability.  
   
   To show just the most salient items, use the cutoff option in \code{\link{print.psych}} }
 \item{cluster.fit, structure.fit, pattern.fit}{There are a number of ways to evaluate how well any factor or cluster matrix reproduces the original matrix. Cluster fit considers how well the clusters fit if only correlations with clusters are considered.  Structure fit evaluates R = CC' while pattern fit evaluate R = C  inverse (phi) C'  where C is the cluster loading matrix, and phi is the intercluster correlation matrix.}
 \item{pattern}{The pattern matrix loadings.   Pattern is just C  inverse (Phi).  The pattern matrix is conceptually equivalent to that of a factor analysis, in that the pattern coefficients are b weights of the cluster to the variables, while the normal cluster loadings are correlations of the items with the cluster.  The four cluster and four factor pattern matrices for the Harman problem are very similar.}
 
\item{order}{This is a vector of the variable names in the order of the cluster diagram. This is useful as a tool for sorting correlations matrices.  See the last example.}
 
 \item{keys}{A list of keys for scoring. These are the keys representing the "purified" clusters }
 \item{keys.org}{A list of scoring keys for the orginal (unpurified) solution.}
 \item{reliability}{The output of the \code{\link{reliability}} function applied to each cluster using the keys list.   This is a useful to examine splithalf reliablities -- which may or may not be the same as beta.   }
 \item{stats}{A list of various fit statistics.  Similar to those from \code{\link{fa}}.  }
}



\references{
Reise, S.P. and Haviland, M. G. (2024) Understanding Alpha and Beta and Sources of Common Variance: Theoretical Underpinnings and a Practical Example. Journal of Personality Assessment <doi.org:10.1080/00223891.2024.2420175>

Revelle, W. Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 1979, 14, 57-74.

Revelle, W. and Zinbarg, R. E. (2009) Coefficients alpha, beta, omega and the glb: comments on Sijtsma. Psychometrika, 2009. 

\url{https://personality-project.org/revelle/publications/iclust.pdf} \cr
See also  more extensive documentation at 
\url{https://personality-project.org/r/r.ICLUST.html} and \cr
Revelle, W. (in prep) An introduction to psychometric theory with applications in R. To be published by Springer.  (working draft available at  \url{https://personality-project.org/r/book/} }
\author{William Revelle   }
        
\note{iclust draws graphical displays with or without using Rgraphiviz. Because of difficulties installing Rgraphviz on many systems, the default it not even try using it.  With the introduction of the \code{\link{diagram}} functions, iclust now draws using \code{\link{diagram}} which is not as pretty as using Rgraphviz, but more stable. However, Rgraphviz can be used by using \code{\link{ICLUST.rgraph}} to  produces slightly better graphics. It is also possible to export dot code in the dot language for further massaging of the graphic. This may be done using \code{\link{ICLUST.graph}}.  This last option is probably preferred for nice graphics which can be massaged in any dot code program (e.g., graphviz (https://graphviz.org)  or a commercial program such as OmniGraffle. 


To view the cluster structure more closely, it is possible to save the graphic output as a pdf and then magnify this using a pdf viewer.  This is useful when clustering a large number of variables.

In order to sort the clusters by cluster loadings, use \code{\link{iclust.sort}}.

By default, the correlations used for the similarity matrix  are Pearson correlations. It is of course possible to use \code{\link{tetrachoric}} or \code{\link{polychoric}} to form the correlation matrix for later analysis.   This option is now included directly in the function and can be specified with the cor option(e.g., cor="tet").

iclust can also be used to organize complex correlation matrices.  Thus, by clustering the items in a correlation matrix, sorting that matrix by the cluster loadings using \code{\link{mat.sort}}, and then plotting with  \code{\link{corPlot}}.  See the penultimate example.  An alterative way of ordering the variables is to use the order object which is just a vector of the variable names sorted by the way they appear in the tree diagram.  (See the final example.)
}
\seealso{ \code{\link{iclust.sort}}, \code{\link{ICLUST.diagram}}, \code{\link{ICLUST.graph}}, \code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{omega}} }

\examples{
test.data <- Harman74.cor$cov
ic.out <- iclust(test.data,title="ICLUST of the Harman data")
summary(ic.out)

#use all defaults and stop at 4 clusters
ic.out4 <- iclust(test.data,nclusters =4,title="Force 4 clusters", reliability=FALSE)  
summary(ic.out4)
ic.out1 <- iclust(test.data,beta=3,beta.size=3)  #use more stringent criteria
ic.out  #more complete output 
plot(ic.out4)    #this shows the spatial representation
#use a dot graphics viewer on the out.file
#dot.graph <- ICLUST.graph(ic.out,out.file="test.ICLUST.graph.dot")  
#show the equivalent of a factor solution 
fa.diagram(ic.out4$pattern,Phi=ic.out4$Phi,main="Pattern taken from iclust") 


#Try running the above example with the Harman74.cor data set from datasets.

 
#demonstrating the use of labels using the bfi data set 

ic<- iclust(bfi[,1:25], labels=bfi.dictionary[1:25,2] ,
    title="ICLUST of bfi data set")
ic  #show the output
cluster.scores <- scoreItems(ic$keys, bfi[,1:25])  #find cluster scores
original.scores <- scoreItems(ic$keys.org,bfi[,1:25]) #find cluster scores

#The effect of a more stringent beta criterion
ic.beta3 <- iclust(bfi[1:25], beta=3, title="ICLUST with more stringent beta")
ic.beta3

#beta is not the worst split half, but rather a general factor estimate
#consider 9 items representing 3 factors
F <- matrix(c(rep(.6,3),rep(0,9),rep(.6,3),rep(0,9),rep(.6,3)), ncol=3)
R <- F \%*\% t(F)
diag(R) <- 1
ic <- iclust(R)
ic$beta     #0 because there is nothing in common between these three clusters
#but
splitHalf(R) # split half = .19  because the split is formed from V1.. V4 and V5 .. V9.
    
#organize a correlation matrix based upon cluster solution.
R <- cor(bfi, use="pairwise")
ic <- iclust(R,plot=FALSE, reliability = FALSE) #suppress the plot
R.s <- mat.sort(R, ic)
corPlot(R.s, main ="bfi sorted by iclust loadings")
#compare with 
corPlot(R[ic$order,ic$order] ,main="bfi sorted by iclust order")




}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
