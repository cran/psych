\name{multilevel.reliability}
\alias{mlr}
\alias{multilevel.reliability}
\alias{mlArrange}
\alias{mlPlot}
\title{Find and plot various reliability/gneralizability coefficients for multilevel data
}
\description{
Various indicators of reliability of multilevel data (e.g., items over time nested within subjects) may be found using generalizability theory.  A basic three way anova is applied to the data from which variance components are extracted. Random effects for a nested design are found by lme.    These are, in turn, converted to several reliability/generalizability coefficients.  An optional call to lme4 to use lmer may be used for unbalanced designs with missing data. mlArrange is a  helper function to convert wide to long format.  Data can be rearranged from wide to long format, and multiple lattice plots of observations overtime for multiple variables and multiple subjects are created.
}
\usage{
mlr(x, grp = "id", Time = "time", items = c(3:5),alpha=TRUE,icc=FALSE, aov=TRUE,
      lmer=FALSE,lme = TRUE,long=FALSE,values=NA,na.action="na.omit",plot=FALSE,
        main="Lattice Plot by subjects over time")
mlArrange(x, grp = "id", Time = "time", items = c(3:5),extra=NULL)
mlPlot(x, grp = "id", Time = "time", items = c(3:5),extra=NULL, 
   col=c("blue","red","black","grey"),type="b",
    main="Lattice Plot by subjects over time",...)
multilevel.reliability(x, grp = "id", Time = "time", items = c(3:5),alpha=TRUE,icc=FALSE,
 aov=TRUE,lmer=FALSE,lme = TRUE,long=FALSE,values=NA,na.action="na.omit",
   plot=FALSE,main="Lattice Plot by subjects over time") #alias for mlr
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame with persons, time, and items.}

  \item{grp}{Which variable specifies people (groups)}

  \item{Time}{Which variable specifies the temporal sequence?}
  \item{items}{Which items should be scored?  Note that if there are multiple scales, just specify the items on one scale at a time.  An item to be reversed scored can be specified by a minus sign. If long format, this is the column specifying item number. }
  \item{alpha}{If TRUE, report alphas for every subject (default)}
  \item{icc}{If TRUE, find ICCs for each person -- can take a while}
  \item{aov}{if FALSE, and if icc  is FALSE, then just draw the within subject plots}
  \item{lmer}{Should we use the lme4 package and lmer or just do the ANOVA?  Requires the
       lme4 package to be installed.  Necessary to do crossed designs with missing data but takes a very long time.}
 \item{lme}{If TRUE, will find the nested components of variance.  Relatively fast.}
   \item{long}{Are the data in wide (default) or long format.}
   \item{values}{If the data are in long format, which column name (number) has the values to be analyzed?}
   \item{na.action}{How to handle missing data.  Passed to the lme function. }   
  \item{plot}{If TRUE, show a lattice plot of the data by subject}
  \item{extra}{Names or locations of extra columns to include in the long output.  These will be carried over from the wide form and duplicated for all items. See example.}
  \item{col}{Color for the lines in mlPlot.  Note that items are categorical and thus drawn in alphabetical order. Order the colors appropriately.}
  \item{type}{The standard type for lines (p,l, b)}
  \item{main}{The main title for the plot (if drawn)}
  \item{...}{Other parameters to pass to xyplot}
}
\details{
Classical reliabiiity theory estimates the amount of variance in a set of observations due to a true score that varies over subjects.  Generalizability theory extends this model to include other sources of variance, specifically, time.  The classic studies using this approach are people measured over multiple time points with multiple items.  Then the question is, how stable are various individual differences. Intraclass correlations (ICC) are found for each subject over items, and for each subject over time. Alpha reliabilities  are found for each subject for the items across time.   

More importantly, components of variance for people, items, time, and their interactions are found either by classical analysis of variance (aov) or by multilevel mixed effect modeling (lme).  These are then used to form several different estimates of generalizability.   Very thoughtful discussions of these procedure may be found in chapters by Shrout and Lane.  

The variance components are the Between Person Variance \eqn{\sigma^2_P}, the variance between items \eqn{\sigma^2_I}, over time \eqn{\sigma^2_T},  and their interactions. 

Then, \eqn{RKF} is the  reliability of average of all ratings across all items and  times (Fixed time effects). (Shrout and Lane, Equation 6): 

\deqn{R_{kF} = \frac{\sigma^2_P + \sigma^2_{PI}/n.I}{\sigma^2_P + \sigma^2_{PI}/n.I + \sigma^2_e/(n.I n.P}}{Rkf = (\sigma^2_P + \sigma^2_{PI}/n.I)/(\sigma^2_P + \sigma^2_{PI}/n.I + \sigma^2_e/(n.I n.P))} 

The generalizability of a single time point across all items (Random time effects) is just

\deqn{R_{1R} = \frac{\sigma^2_P + \sigma^2_{PI}/n.I}{\sigma^2_P + \sigma^2_{PI}/n.I + \sigma^2_T + \sigma^2_{PT}+ \sigma^2_e/(n.I)}}{R1R = (\sigma^2_P + \sigma^2_{PI}/n.I)/(\sigma^2_P +  \sigma^2_{PI}/n.I + \sigma^2_T + \sigma^2_{PT} \sigma^2_e/n.I)}
(Shrout and Lane equation 7 with a correction per Sean Lane.)

Generalizability of average time points across all items (Random effects). (Shrout and Lane, equation 8)
\deqn{R_{kR} = \frac{\sigma^2_P + \sigma^2_{PI}/n.I}{\sigma^2_P + \sigma^2_{PI}/n.I + \sigma^2_T/n.T + \sigma^2_{PT}/n.T+ \sigma^2_e/n.I}}{RkR = (\sigma^2_P + \sigma^2_{PI}/n.I)/(\sigma^2_P +  \sigma^2_{PI}/n.I + \sigma^2_T/n.T + \sigma^2_{PT}/n.T + \sigma^2_e/(n.I n.T))}

Generalizability of change scores (Shrout and Lane, equation  9)
\deqn{R_{C} = \frac{\sigma^2_{PT}}{\sigma^2_{PT} + \sigma^2_e/n.I}}{RC = (\sigma^2_PT)/(\sigma^2_PT +  \sigma^2_e/(n.I))}.


If the design may be thought of as fully crossed, then either aov or lmer can be used to estimate the components of variance.  With no missing data and a balanced design, these will give identical answers. However aov breaks down with missing data and seems to be very slow and very memory intensive for large problems ( 5,919  seconds for 209 cases with with 88 time points and three items on a Mac Powerbook with a 2.8 GHZ Intel Core I7). The slowdown probably is memory related, as the memory demands increased to 22.62 GB of compressed memory.   lmer will handle this design but is not nearly as slow  (242 seconds for the 209 cases with 88 time points and three items) as the aov approach.   
 
If the design is thought of as nested, rather than crossed, the components of variance are found using the lme function from nlme. This is very fast (114 cases with 88 time points and three items took 3.5 seconds). 


The nested design leads to the generalizability of K random effects Nested (Shrout and Lane, equation 10):

\deqn{R_{KRN} = \frac{\sigma^2_P }{\sigma^2_P + \sigma^2_{T(P)}/n.I + \sigma^2_e/(n.I n.P}}{RKkRN = (\sigma^2_P)/(\sigma^2_P + \sigma^2_{T(P)}/n.p + \sigma^2_e/(n.I n.T))} 
 
And, finally, to the reliability of between person differences, averaged over items.  (Shrout and Lane, equation 11).

\deqn{R_{CN} = \frac{\sigma^2_{T(P)} }{\sigma^2_{T(P)} + \sigma^2_e/(n.I}}{RCN = (\sigma^2_T(P)/(\sigma^2_T(P) + \sigma^2_e/(n.I))} 

Unfortunately, when doing the nested analysis, lme will sometimes issue an obnoxious error about failing to converge.  To fix this, turning off lme  and just using lmer seems to solve the problem (i.e., set lme=FALSE and lmer=TRUE).  (lme is part of core R and its namespace is automatically attached when loading \code{\link{psych}}). For many problems, lmer is not necessary and is thus not loaded.  However sometimes it is useful.  To use lmer it is necessary to have the lme4 package installed.  It will be automatically loaded if it is installed and requested. In the interests of making a 'thin' package, lmer is suggested,not required.

The input can either be in 'wide' or 'long' form.  If in wide form, then specify the grouping variable, the 'time' variable, and the the column numbers or names of the items. (See the first example).  If in  long format, then what is the column (name or number) of the dependent variable.  (See the second example.)


\code{\link{mlArrange}} takes a wide data.frame and organizes it into a `long' data.frame suitable for a lattice xyplot.  This is a convenient alternative to \code{\link{stack}}, particularly for unbalanced designs.  The wide data frame is reorganized into a long data frame organized by grp (typically a subject id), by Time (typically a time varying variable, but can be anything, and then stacks the items within each person and time.  Extra variables are carried over and matched to the appropriate grp  and Time. 

Thus, if we have N subjects over t time points for k items, in wide format for N * t rows where each row has k items and e extra pieces of information, we get a N x t * k row by 4 + e column dataframe.  The first four columns in the long output are id, time, values, and item names, the remaining columns are the extra values.  These  could be something such as a trait measure for each subject, or the situation in which the items are given.

\code{\link{mlArrange}} plots k items over the  t time dimensions for each subject. 
}

\value{
\item{n.obs}{Number of individuals}
\item{n.time}{Maximum number of time intervals}
\item{n.items}{Number of items}
 \item{components}{Components of variance associated with individuals, Time, Items, and their interactions.}
 \item{RkF }{Reliability of average of all ratings across all items and times (fixed effects).}
\item{R1R}{Generalizability of a single time point across all items (Random effects)}
\item{RkR}{Generalizability of average time points across all items (Random effects)}
\item{Rc}{Generalizability of change scores over time.}
\item{RkRn}{ Generalizability of between person differences averaged over time and items}  
\item{Rcn}{Generalizability of within person variations averaged over items (nested structure)}
\item{ANOVA}{The summary anova table from which the components are found (if done),}
\item{s.lmer}{The summary of the lmer analysis (if done),}
\item{s.lme}{The summary of the lme analysis (if done),}
\item{alpha}{Within subject alpha over items and time.}
\item{summary.by.person}{Summary table of ICCs organized by person,}
\item{summary.by.time}{Summary table of ICCs organized by time.}
\item{ICC.by.person}{A rather long list of ICCs by person.}
\item{ICC.by.time}{Another long list of ICCs, this time for each time period,}
\item{long}{The data (x) have been rearranged into long form for graphics or for further analyses using lme, lmer, or aov that require long form.}


}
\references{
Bolger, Niall and Laurenceau, Jean-Phillippe, (2013) Intensive longitudinal models.  New York.  Guilford Press. 

Cranford, J. A., Shrout, P. E., Iida, M., Rafaeli, E., Yip, T., & Bolger, N. (2006). A procedure for evaluating sensitivity to within-person change: Can mood measures in diary studies detect change reliably? Personality and Social Psychology Bulletin, 32(7), 917-929.

Revelle, W. and Condon, D.M. (2019) Reliability from alpha to omega: A tutorial.  Psychological Assessment, 31, 12, 1395-1411. https://doi.org/10.1037/pas0000754.  \url{https://personality-project.org/revelle/publications/rc.pa.19.pdf} Preprint available from personality-proejct. 

Revelle, W. and Wilt, J. (2017) Analyzing dynamic data: a tutorial. Personality and Individual Differences. DOI: 10.1016/j.paid.2017.08.020

Shrout, Patrick and Lane, Sean P (2012), Psychometrics.  In M.R. Mehl and T.S. Conner (eds)  Handbook of research methods for studying daily life, (p 302-320) New York. Guilford Press


}
\author{William Revelle}



\seealso{\code{\link{sim.multi}} and \code{\link{sim.multilevel}} to generate multilevel data, \code{\link{statsBy}} a for statistics for multi level analysis.
}
\examples{
#data from Shrout and Lane, 2012.

shrout <- structure(list(Person = c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L), Time = c(1L, 1L, 
1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 
4L, 4L), Item1 = c(2L, 3L, 6L, 3L, 7L, 3L, 5L, 6L, 3L, 8L, 4L, 
4L, 7L, 5L, 6L, 1L, 5L, 8L, 8L, 6L), Item2 = c(3L, 4L, 6L, 4L, 
8L, 3L, 7L, 7L, 5L, 8L, 2L, 6L, 8L, 6L, 7L, 3L, 9L, 9L, 7L, 8L
), Item3 = c(6L, 4L, 5L, 3L, 7L, 4L, 7L, 8L, 9L, 9L, 5L, 7L, 
9L, 7L, 8L, 4L, 7L, 9L, 9L, 6L)), .Names = c("Person", "Time", 
"Item1", "Item2", "Item3"), class = "data.frame", row.names = c(NA, 
-20L))

#make shrout super wide
#Xwide <- reshape(shrout,v.names=c("Item1","Item2","Item3"),timevar="Time", 
#direction="wide",idvar="Person")
#add more helpful Names
#colnames(Xwide ) <- c("Person",c(paste0("Item",1:3,".T",1),paste0("Item",1:3,".T",2), 
#paste0("Item",1:3,".T",3),paste0("Item",1:3,".T",4)))
#make superwide into normal form  (i.e., just return it to the original shrout data
#Xlong <-Xlong <- reshape(Xwide,idvar="Person",2:13)

#Now use these data for a multilevel repliability study, use the normal wide form output
mg <- mlr(shrout,grp="Person",Time="Time",items=3:5) 
#which is the same as 
#mg <- multilevel.reliability(shrout,grp="Person",Time="Time",items=
#         c("Item1","Item2","Item3"),plot=TRUE)
#to show the lattice plot by subjects, set plot = TRUE

#Alternatively for long input (returned in this case from the prior run)
mlr(mg$long,grp="id",Time ="time",items="items", values="values",long=TRUE)

#example of mlArrange
#First, add two new columns to shrout and 
#then convert to long output using mlArrange
total <- rowSums(shrout[3:5])
caseid <- rep(paste0("ID",1:5),4)
new.shrout <- cbind(shrout,total=total,case=caseid)
#now convert to long
new.long <- mlArrange(new.shrout,grp="Person",Time="Time",items =3:5,extra=6:7)
headTail(new.long,6,6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }% __ONLY ONE__ keyword per line
