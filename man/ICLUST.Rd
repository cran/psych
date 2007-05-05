\name{ICLUST}
\alias{ICLUST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ICLUST: Item Cluster Analysis   -- Hierarchical cluster analysis using psychometric principles }
\description{A common data reduction technique is to cluster cases (subjects). Less common, but particularly useful in psychological research, is to cluster items (variables). This may be thought of as an alternative to factor analysis, based upon a much simpler model. The cluster model is that the correlations between variables reflect that each item loads on at most one cluster, and that items that load on those clusters correlate as a function of their respective loadings on that cluster and items that define different clusters correlate as a function of their respective cluster loadings and the intercluster correlations. 
Essentially, the cluster model is a Very Simple Structure factor model of complexity one (see \code{\link{VSS}}).

This function applies the ICLUST algorithm to hierarchically cluster items to form composite scales. Clusters are combined if coefficients alpha and beta will increase in the new cluster.

Alpha, the mean split half correlation, and beta, the worst split half correlation, are estimates of the reliability and general factor saturation of the test.  (See also the \code{\link{omega}} function to estimate McDonald's coeffient omega.)


}
\usage{
ICLUST(r.mat, nclusters=1, alpha=3, beta=1, beta.size=4, alpha.size=3,
correct=TRUE, reverse=TRUE, beta.min=.5, output=1, digits=2,labels=NULL,cut=0,
n.iterations = 0,title="ICLUST")

#ICLUST(r.mat)    #use all defaults
#ICLUST(r.mat,nclusters =3)    #use all defaults and if possible stop at 3 clusters
#ICLUST(r.mat, output =3)     #long output shows clustering history
#ICLUST(r.mat, n.iterations =3)  #clean up solution by item reassignment
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.mat}{ A correlation matrix or data matrix/data.frame. (If r.mat is not square i.e, a correlation matrix, the data are correlated using pairwise deletion. }
  \item{nclusters}{Extract clusters until nclusters remain (default =1)  }
  \item{alpha}{ Apply the increase in alpha criterion  (0) never or for (1) the smaller, 2) the average, or 3) the greater of the separate alphas. (default = 3) }
  \item{beta}{  Apply the increase in beta criterion (0) never or for (1) the smaller, 2) the average, or 3) the greater of the separate betas. (default =1) }
  \item{beta.size}{ Apply the beta criterion after clusters are of beta.size  (default = 4)}
  \item{alpha.size}{ Apply the alpha criterion after clusters are of size alpha.size (default =3) }
  \item{correct}{ Correct correlations for reliability (default = TRUE) }
  \item{reverse}{Reverse negative keyed items (default = TRUE}
  \item{beta.min}{ Stop clustering if the beta is not greater than beta.min (default = .5) }
  \item{output}{ 1) short, 2) medium, 3 ) long output (default =1)}
  \item{labels}{vector of item content or labels}
  \item{cut}{sort cluster loadings > absolute(cut) (default = 0) }
  \item{n.iterations}  {iterate the solution n.iterations times to "purify" the clusters (default = 0)}
  \item{digits}{ Precision of digits of output (default = 2) }
  \item{title}{ Title for this run }
}
\details{ Extensive documentation and justification of the algorithm is available in the original MBR 1979 \url{http://personality-project.org/revelle/publications/iclust.pdf} paper.  Further discussion of the algorithm and sample output is available on the personality-project.org web page: \url{http://personality-project.org/r/r.ICLUST.html} 

The results are best visualized using  \code{\link{ICLUST.graph}}, the results of which can be saved as a dot file for the Graphviz program.  \url{http://www.graphviz.org/}

 A common problem in the social sciences is to construct scales or composites of items to measure constructs of theoretical interest and practical importance. This process frequently involves administering a battery of items from which those that meet certain criteria are selected. These criteria might be rational, empirical,or factorial. A similar problem is to analyze the adequacy of scales that already have been formed and to decide whether the putative constructs are measured properly. Both of these problems have been discussed in numerous texts,  as well as in myriad articles. Proponents of various methods have argued for the importance of face validity, discriminant validity, construct validity, factorial homogeneity, and theoretical importance. 
 
Revelle (1979) proposed that hierachical cluster analysis could be used to estimate a new coefficient (beta) that was an estimate of the  general factor saturation of a test.  More recently, Zinbarg, Revelle, Yovel and Li (2005) compared McDonald's Omega to Chronbach's alpha and Revelle's beta. They conclude that omega is the best estimate.  An algorithm for estimating \code{\link{omega} } is available as part of this package. 


This R version is a completely new version of ICLUST.  Although early testing suggests it is stable, let me know if you have problems. Please email me if you want help with this version of ICLUST or if you desire more features.

The program currently has three primary functions: cluster, loadings, and graphics.  

Clustering 24 tests of mental ability

 A sample output using the 24 variable problem by Harman can be represented both graphically and in terms of the cluster order.  Note that the graphic is created using GraphViz in the dot language.  \code{\link{ICLUST.graph}} produces the dot code for Graphviz.  I have also tried doing this in Rgraphviz with less than wonderful results.  Dot code can be viewed directly in Graphviz or can be tweaked using commercial software packages (e.g.,OmniGraffle)

Note that for this problem, with these parameters, the data formed one large cluster. (This is consistent with the Very Simple Structure (\code{\link{VSS}}) output as well, which shows a clear one factor solution for complexity 1 data.)  See below for an example with this same data set, but with more stringent parameter settings.  

To see the graphic output go to \url{http://personality-project.org/r/r.ICLUST.html} . 

}
\value{
  \item{title }{Name of this run}
  \item{results}{A list containing}
   \item{clusters}{a matrix of -1,0, and 1 values to define cluster membership.}
  \results {The step by step cluster history, including which pair was grouped, what were the alpha and betas of the two groups and of the combined group}
   \item{corrected}{The raw and corrected for alpha reliability cluster intercorrelations.}
   \item{purified}{A list of the cluster definitions and cluster loadings of the purified solution}
}

\references{Revelle, W. Hierarchical Cluster Analysis and the Internal Structure of Tests. Multivariate Behavioral Research, 1979, 14, 57-74. 
\url{http://personality-project.org/revelle/publications/iclust.pdf} \cr
See also  more extensive documentation at 
\url{http://personality-project.org/r/r.ICLUST.html} }
\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illinois \cr
        \email{  revelle@northwestern.edu           }\cr
 \url{http://personality-project.org/revelle.html}   }
        

\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{omega} }}

\examples{
\dontrun{
test.data <- Harman74.cor$cov
ic.out <- ICLUST(test.data)         #use all defaults
out.file <- file.choose(new=TRUE)   #create a new file to write the plot commands to 
ICLUST.graph(ic.out,out.file,title = "ICLUST of Harman's 24 mental variables" )   

ic.out <- ICLUST(test.data,nclusters =3)  #use all defaults and if possible stop at 3 clusters
ICLUST.graph(ic.out,out.file,title = "ICLUST of 24 mental variables with forced 3 cluster solution")

ICLUST(test.data, output =3)     #long output shows clustering history

ic.out <- ICLUST(test.data,,nclusters=4, n.iterations =3)  #clean up solution by item reassignment
ICLUST.graph(ic.out,out.file,title = "ICLUST of 24 mental variables with forced 4 cluster solution")
ic.out    #shows the output on the console
}

#produces this output
#ICLUST(Harman74.cor$cov)
#$title
#[1] "ICLUST"
#
#$clusters
#      VisualPerception                  Cubes         PaperFormBoard                  Flags     GeneralInformation 
#                     1                      1                      1                      1                      1 
# PargraphComprehension     SentenceCompletion     WordClassification            WordMeaning               Addition 
#                     1                      1                      1                      1                      1 
#                  Code           CountingDots StraightCurvedCapitals        WordRecognition      NumberRecognition 
#                     1                      1                      1                      1                      1 
#     FigureRecognition           ObjectNumber           NumberFigure             FigureWord              Deduction 
#                     1                      1                      1                      1                      1 
#      NumericalPuzzles       ProblemReasoning       SeriesCompletion     ArithmeticProblems 
#                     1                      1                      1                      1 
#
#$corrected
#     [,1]
#[1,]    1
#
#$loadings
#                       [,1]
#VisualPerception       0.57
#Cubes                  0.36
#PaperFormBoard         0.40
#Flags                  0.46
#GeneralInformation     0.62
#PargraphComprehension  0.62
#SentenceCompletion     0.60
#WordClassification     0.63
#WordMeaning            0.62
#Addition               0.43
#Code                   0.54
#CountingDots           0.44
#StraightCurvedCapitals 0.57
#WordRecognition        0.41
#NumberRecognition      0.38
#FigureRecognition      0.50
#ObjectNumber           0.45
#NumberFigure           0.51
#FigureWord             0.44
#Deduction              0.59
#NumericalPuzzles       0.58
#ProblemReasoning       0.58
#SeriesCompletion       0.66
#ArithmeticProblems     0.62
#
#$fit
#$fit$clusterfit
#[1] 0.78
#
#$fit$factorfit
#[1] 0.78
#
#
#$results
#    Item/Cluster Item/Cluster similarity correlation alpha1 alpha2 beta1 beta2 size1 size2 rbar1 rbar2   r1   r2 alpha
#C1           V23          V20       1.00        0.51   0.51   0.51  0.51  0.51     1     1  0.51  0.51 0.59 0.59  0.68
#C2            V9           V5       1.00        0.72   0.72   0.72  0.72  0.72     1     1  0.72  0.72 0.78 0.78  0.84
#C3            V7           V6       1.00        0.72   0.72   0.72  0.72  0.72     1     1  0.72  0.72 0.78 0.78  0.84
#C4           V12          V10       1.00        0.58   0.58   0.58  0.58  0.58     1     1  0.58  0.58 0.65 0.65  0.73
#C5           V13          V11       1.00        0.54   0.54   0.54  0.54  0.54     1     1  0.54  0.54 0.62 0.62  0.70
#C6           V18          V17       1.00        0.45   0.45   0.45  0.45  0.45     1     1  0.45  0.45 0.53 0.53  0.62
#C7            V4           V1       0.99        0.47   0.47   0.48  0.47  0.48     1     1  0.47  0.48 0.55 0.55  0.64
#C8           V16          V14       0.98        0.41   0.43   0.41  0.43  0.41     1     1  0.43  0.41 0.50 0.49  0.58
#C9            C2           C3       0.93        0.78   0.84   0.84  0.84  0.84     2     2  0.72  0.72 0.86 0.86  0.90
#C10           C1          V22       0.91        0.56   0.67   0.56  0.68  0.56     2     1  0.51  0.56 0.71 0.63  0.75
#C11          V21          V24       0.87        0.45   0.51   0.53  0.51  0.53     1     1  0.51  0.53 0.56 0.58  0.62
#C12          C10          C11       0.86        0.58   0.74   0.62  0.72  0.62     3     2  0.49  0.45 0.76 0.67  0.79
#C13           C9           V8       0.84        0.64   0.90   0.64  0.88  0.64     4     1  0.69  0.64 0.90 0.68  0.90
#C14           C8          V15       0.84        0.41   0.58   0.41  0.58  0.41     2     1  0.41  0.41 0.61 0.48  0.63
#C15           C5           C4       0.82        0.59   0.70   0.74  0.70  0.73     2     2  0.54  0.58 0.72 0.74  0.80
#C16           V3           V2       0.81        0.32   0.41   0.38  0.41  0.38     1     1  0.41  0.38 0.45 0.43  0.48
#C17          C16           C7       0.81        0.45   0.48   0.64  0.48  0.64     2     2  0.32  0.47 0.55 0.64  0.67
#C18          C12          C17       0.81        0.59   0.79   0.67  0.73  0.62     5     4  0.43  0.34 0.79 0.70  0.83
#C19          V19           C6       0.80        0.40   0.40   0.62  0.40  0.62     1     2  0.40  0.45 0.47 0.64  0.64
#C20          C19          C14       0.77        0.49   0.64   0.64  0.57  0.58     3     3  0.38  0.37 0.66 0.65  0.74
#C21          C18          C20       0.74        0.58   0.83   0.74  0.74  0.66     9     6  0.35  0.32 0.82 0.72  0.86
#C22          C21          C13       0.70        0.62   0.86   0.90  0.73  0.78    15     5  0.29  0.64 0.86 0.78  0.90
#C23          C22          C15       0.65        0.55   0.90   0.79  0.77  0.74    20     4  0.31  0.49 0.90 0.65  0.91
#    beta rbar size
#C1  0.68 0.51    2
#C2  0.84 0.72    2
#C3  0.84 0.72    2
#C4  0.73 0.58    2
#C5  0.70 0.54    2
#C6  0.62 0.45    2
#C7  0.64 0.47    2
#C8  0.58 0.41    2
#C9  0.88 0.69    4
#C10 0.72 0.49    3
#C11 0.62 0.45    2
#C12 0.73 0.43    5
#C13 0.78 0.64    5
#C14 0.58 0.37    3
#C15 0.74 0.49    4
#C16 0.48 0.32    2
#C17 0.62 0.34    4
#C18 0.74 0.35    9
#C19 0.57 0.38    3
#C20 0.66 0.32    6
#C21 0.73 0.29   15
#C22 0.77 0.31   20
#C23 0.71 0.30   24
#
#$cor
#     [,1]
#[1,]    1
#
#$alpha
#[1] 0.91
#
#$size
#[1] 24
#
#$sorted
#$sorted$sorted
#                       item                content cluster loadings
#SeriesCompletion         23       SeriesCompletion       1     0.66
#WordClassification        8     WordClassification       1     0.63
#GeneralInformation        5     GeneralInformation       1     0.62
#PargraphComprehension     6  PargraphComprehension       1     0.62
#WordMeaning               9            WordMeaning       1     0.62
#ArithmeticProblems       24     ArithmeticProblems       1     0.62
#SentenceCompletion        7     SentenceCompletion       1     0.60
#Deduction                20              Deduction       1     0.59
#NumericalPuzzles         21       NumericalPuzzles       1     0.58
#ProblemReasoning         22       ProblemReasoning       1     0.58
#VisualPerception          1       VisualPerception       1     0.57
#StraightCurvedCapitals   13 StraightCurvedCapitals       1     0.57
#Code                     11                   Code       1     0.54
#NumberFigure             18           NumberFigure       1     0.51
#FigureRecognition        16      FigureRecognition       1     0.50
#Flags                     4                  Flags       1     0.46
#ObjectNumber             17           ObjectNumber       1     0.45
#CountingDots             12           CountingDots       1     0.44
#FigureWord               19             FigureWord       1     0.44
#Addition                 10               Addition       1     0.43
#WordRecognition          14        WordRecognition       1     0.41
#PaperFormBoard            3         PaperFormBoard       1     0.40
#NumberRecognition        15      NumberRecognition       1     0.38
#Cubes                     2                  Cubes       1     0.36
#
#
#$p.fit
#$p.fit$clusterfit
#[1] 0.78
#
#$p.fit$factorfit
#[1] 0.78
#
#
#$p.sorted
#$p.sorted$sorted
#                       item                content cluster loadings
#SeriesCompletion         23       SeriesCompletion       1     0.66
#WordClassification        8     WordClassification       1     0.63
#GeneralInformation        5     GeneralInformation       1     0.62
#PargraphComprehension     6  PargraphComprehension       1     0.62
#WordMeaning               9            WordMeaning       1     0.62
#ArithmeticProblems       24     ArithmeticProblems       1     0.62
#SentenceCompletion        7     SentenceCompletion       1     0.60
#Deduction                20              Deduction       1     0.59
#NumericalPuzzles         21       NumericalPuzzles       1     0.58
#ProblemReasoning         22       ProblemReasoning       1     0.58
#VisualPerception          1       VisualPerception       1     0.57
#StraightCurvedCapitals   13 StraightCurvedCapitals       1     0.57
#Code                     11                   Code       1     0.54
#NumberFigure             18           NumberFigure       1     0.51
#FigureRecognition        16      FigureRecognition       1     0.50
#Flags                     4                  Flags       1     0.46
#ObjectNumber             17           ObjectNumber       1     0.45
#CountingDots             12           CountingDots       1     0.44
#FigureWord               19             FigureWord       1     0.44
#Addition                 10               Addition       1     0.43
#WordRecognition          14        WordRecognition       1     0.41
#PaperFormBoard            3         PaperFormBoard       1     0.40
#NumberRecognition        15      NumberRecognition       1     0.38
#Cubes                     2                  Cubes       1     0.36
#
#
#$purified
#$purified$cor
#     [,1]
#[1,]    1
#
#$purified$sd
#[1] 13.79
#
#$purified$corrected
#     [,1]
#[1,] 0.91
#
#$purified$size
#[1] 24
#
#
### The function is currently defined as
function (r.mat,nclusters=1,alpha=3,beta=2,beta.size=4,alpha.size=3,correct=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2,labels=NULL,cut=0,n.iterations=0,title="ICLUST") {#should allow for raw data, correlation or covariances
 #ICLUST.options <- list(n.clus=1,alpha=3,beta=2,beta.size=4,alpha.size=3,correct=TRUE,reverse=TRUE,beta.min=.5,output=1,digits=2) 
	ICLUST.options <- list(n.clus=nclusters,alpha=alpha,beta=beta,beta.size=beta.size,alpha.size=alpha.size,correct=correct,reverse=reverse,beta.min=beta.min,output=output,digits=digits) 
	if(dim(r.mat)[1]!=dim(r.mat)[2]) {r.mat <- cor(r.mat,use="pairwise") }    #cluster correlation matrices, find correlations if not square matrix
	iclust.results <- ICLUST.cluster(r.mat,ICLUST.options)
	load <- cluster.loadings(iclust.results$clusters,r.mat,digits=digits)
	fits <- cluster.fit(r.mat,load$loadings,iclust.results$clusters,digits=digits)
	sorted <- ICLUST.sort(ic.load=load$loadings,labels=labels,cut=cut) #sort the loadings
	
	#now, iterate the cluster solution to clean it up (if desired)
	
		clusters <- iclust.results$clusters
		old.clusters <- clusters
		old.fit <- fits$clusterfit
		clusters <- factor2cluster(load$loadings,cut=cut,loading=FALSE)
		load <- cluster.loadings(clusters,r.mat,digits=digits)
		if (n.iterations > 0) {  #it is possible to iterate the solution to perhaps improve it 
		for (steps in 1:n.iterations) {   #
			load <- cluster.loadings(clusters,r.mat,digits=digits)
			clusters <- factor2cluster(load$loadings,cut=cut,loading=FALSE)
			if(dim(clusters)[2]!=dim(old.clusters)[2]) {change <- 999 
			load <- cluster.loadings(clusters,r.mat,digits=digits) } else {
			change <- sum(abs(clusters)-abs(old.clusters)) }  #how many items are changing?
			fit <- cluster.fit(r.mat,load$loadings,clusters,digits=digits)
		old.clusters <- clusters
		print(paste("iterations ",steps," change in clusters ", change, "current fit " , fit$clusterfit))
		if ((abs(change) < 1) | (fit$clusterfit <= old.fit)) {break}    #stop iterating if it gets worse or there is no change in cluster definitions
		old.fit <- fit$cluster.fit
					}
		}
	p.fit <- cluster.fit(r.mat,load$loadings,clusters,digits=digits)
	p.sorted <- ICLUST.sort(ic.load=load$loadings,labels=labels,cut=cut)
	purified <- cluster.cor(clusters,r.mat,digits=digits)
	list(title=title,clusters=iclust.results$clusters,corrected=load$corrected,loadings=load$loadings,fit=fits,results=iclust.results$results,cor=load$cor,alpha=load$alpha,size=load$size,sorted=sorted,p.fit = p.fit,p.sorted = p.sorted,purified=purified)
}   

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
