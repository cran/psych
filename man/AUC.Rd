\name{AUC}
\alias{AUC}
\alias{auc}
\alias{Specificity}
\alias{Sensitivity}


\title{Decision Theory measures of specificity, sensitivity, and d prime}
\description{In many fields, decisions and outcomes are categorical even though the underlying phenomenon are probably continuous.  E.g. students are accepted to graduate school or not, they finish or not. X-Rays are diagnosed as patients having cancer or not.   Outcomes of such decisions are usually labeled as Valid Positives, Valid Negatives, False Positives and False Negatives. In hypothesis testing, False Positives are known as Type I errors, while False Negatives are Type II errors.  The relationship between these four cells depends upon the correlation between the decision rule and the outcome as well as the level of evidence needed for a decision (the criterion).  Signal Detection Theory and Decision Theory have a number of related measures of performance (accuracy = VP + VN), Sensitivity (VP/(VP + FN)), Specificity (1 - FP), d prime (d'), and the area under the Response Operating Characteristic Curve (AUC). More generally, these are examples of correlations based upon dichotomous data.  \code{\link{AUC}} addresses some of these questions.  


}
\usage{
AUC(t=NULL,BR=NULL,SR=NULL,Phi=NULL,VP=NULL,labels=NULL,plot="b",zero=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{t}{a 4 x 1 vector or a 2 x2 table of TP, FP, FN, TN values (see below) May be 
 counts or proportions.} 
 \item{BR}{Base Rate of successful outcomes or actual symptom (if t is not specified)} 
 \item{SR}{Selection Rate for  candidates or diagnoses (if t is not specified)}
\item{Phi}{The Phi correlation coefficient between the predictor and the outcome variable
    (if t is not specified)}
\item{VP}{The number of Valid Positives (selected applicants who succeed; correct 
diagnoses).(if t and Phi are not specified)}
  \item{labels}{Names of variables  1 and 2}
 \item{plot}{"b" (both), "d" (decision theory), "a" (auc), or "n" neither}
 \item{zero}{If True, then the noise distribution is centered at zero}
}

\value{
  \item{phi}{Phi coefficient of the two by two table}
  \item{tetra}{Tetrachoric (latent) coefficient inferred from the two by two table}
   \item{r.bis}{Biserial correlation of continuous state of world with decision}
  \item{observed}{The observed input (as a check)}
  \item{probabilities}{Observed values/ total number of observations}
  \item{conditional}{prob / rowSums(prob)}
  \item{Accuracy}{percentage of True Positives + True Negatives}
  \item{Sensitivity}{VP/(VP + FN)}
  \item{Specificity}{VN/(FP + VN)}
  \item{d.prime}{difference of True Positives versus True Negatives}
  \item{beta}{ratio of ordinates at the decision point}
  
 
}

\details{The problem of making binary decisions about the state of the world is ubiquitous.  We see this in Null Hypothesis Significance Testing (NHST), medical diagnoses, and selection for occupations.  Variously known as NHST, Signal Detection Theory, clinical Assessment, or college admissions, all of these domains share the same two x two decision task.

Although the underlying phenomena are probably continuous, a typical decision or  diagnostic situation makes dichotomous decisions: Accept or Reject, correctly identified, incorrectly identified.   In Signal Detection Theory, the world has two states: Noise versus Signal + Noise.  The decision is whether there is  a signal or not.  

In diagnoses, it is whether to diagnose an illness or not given some noisy signal (e.g., an X-Ray, a set of diagnostic tests).

In college admissions, we accept some students and reject others.  Four-Five years later we observe who "succeeds" or graduates. 

All of these decisions lead to four cells based upon a two x two categorization.  Given the true state of the world is Positive or Negative, and a rater assigns positive or negative ratings, then the resulting two by two table has True Positives and True Negatives on the diagonal and False Positives and False Negatives off the diagonal. 

When expressed as percentages of the total, then  Base Rates (BR) depend upon the state of the world, but Selection Ratios (SR) are  under the control of the person making the decision and affect the number of False Positives and the number of Valid Positives.

Given a two x two table of counts or percentages \cr
\tabular{lllll}{
\tab       \tab Decide +  \tab Decide - \cr
\tab True + \tab VP \tab FN  \tab BR \cr
\tab True - \tab FP \tab VN  \tab 1- BR \cr
\tab \tab   SR \tab  1 - SR  \tab  (N)
}

Unfortunately, although this way of categorizing the data is typical in assessment (e.g., Wiggins 1973), and everything is expressed as percentages of the total, in some decision papers, VP are expressed as the ratio of VP to total positive decisions (e.g., Wickens, 1984).  This requires dividing through by the column totals (and represented as VP* and FP* in the table below).

The relationships implied by these data can be summarized as a \code{\link{phi}} or \code{\link{tetrachoric}}  correlation between the raters and the world, or as a decision process with several alternative measures:

 Sensitivity, Specificity, Accuracy, Area Under the Curve, and d' (d prime). These measures  may be defined as \cr

 \tabular{llll}{
\tab Measure    \tab Definition   \cr
\tab Sensitivity \tab VP/(VP+ FN) \cr
\tab Specificity\tab VN/(FP + VN)  \cr
\tab Accuracy \tab  VP + VN  \cr
\tab VP* \tab VP/(VP + FP) \cr
\tab FP* \tab (FP/(VP + FP \cr
\tab d'  \tab z(VP*) - z(FP*) \cr
\tab d' \tab sqrt(2) z(AUC) \cr
\tab beta \tab prob(X/S)/(prob(X/N) \cr
}

Although only one point is found, we can form a graphical display of VP versus FP as a smooth curve as a function of the decision criterion. The smooth curve assumes normality whereas the other merely are the two line segments between the points (0,0), (FP,VP), (1,1). The resulting correlation between the inferred continuous state of the world and the dichotomous decision process is a biserial correlation.

When using table input, the values can be counts and thus greater than 1 or merely probabilities which should add up to 1. Base Rates and Selection Ratios are proportions and thus less than 1. }

\author{William Revelle }

\references{

Metz, C.E. (1978) Basic principles of ROC analysis. Seminars in Nuclear Medicine, 8, 283-298.

Wiggins, Jerry S. (1973) Personality and Prediction: Principles of Personality Assessment. Addison-Wesley.

Wickens, Christopher D. (1984) Engineering Psychology and Human Performance.  Merrill.

 }
\seealso{ \code{\link{phi}}, \code{\link{phi2tetra}} ,\code{\link{Yule}}, \code{\link{Yule.inv}} \code{\link{Yule2phi}}, \code{\link{tetrachoric}} and \code{\link{polychoric}}, \code{\link{comorbidity}}}

\examples{
AUC(c(30,20,20,30))  #specify the table input
AUC(c(140,60,100,900))  #Metz example 1
AUC(c(80,120,40, 960)) #Metz example 2  Note how the accuracies are the same but d's differ
AUC(c(49,40,79,336)) #Wiggins p 249
AUC(BR=.05,SR=.254,Phi = .317) #Wiggins 251 extreme Base Rates





}

\keyword{multivariate }

