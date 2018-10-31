\name{r.test}
\alias{r.test}
\title{Tests of significance for correlations}
\description{Tests the significance of a single correlation, the difference between two independent correlations, the difference between two dependent correlations sharing one variable (Williams's Test), or the difference between two dependent correlations with different variables (Steiger Tests).
}
\usage{
r.test(n, r12, r34 = NULL, r23 = NULL, r13 = NULL, r14 = NULL, r24 = NULL, 
          n2 = NULL,pooled=TRUE, twotailed = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size of first group }
  \item{r12}{Correlation to be tested}
  \item{r34}{Test if this correlation is different from r12, if r23 is specified, but r13 is not, then r34 becomes r13  }
  \item{r23}{if ra = r(12) and rb = r(13)  then test for differences of dependent correlations given r23}
  \item{r13}{implies  ra =r(12) and rb =r(34)  test for difference of dependent correlations }
  \item{r14}{implies   ra =r(12) and rb =r(34) }
  \item{r24}{ ra =r(12) and rb =r(34)}
  \item{n2}{n2 is specified in the case of two independent correlations. n2 defaults to n if if not specified }
  \item{pooled}{use pooled estimates of correlations}
  \item{twotailed}{ should a twotailed or one tailed test be used }
}
\details{Depending upon the input, one of four different tests of correlations is done.
1) For a sample size n, find the t value for a single correlation where \deqn{t = \frac{r * \sqrt(n-2)}{\sqrt(1-r^2)}
}{t = r* sqrt(n-2)/sqrt(1-r^2) } and 

\deqn{se = \sqrt{\frac{1-r^2}{n-2}}) }{se = sqrt((1-r^2)/(n-2))}.

2) For sample sizes of n and n2 (n2 = n if not specified) find the z of the difference between the z transformed correlations divided by the standard error of the difference of two z scores: 
\deqn{z = \frac{z_1  - z_2}{\sqrt{\frac{1}{(n_1 - 3) + (n_2 - 3)}}}}{t = (z_1  - z_2) * sqrt(1/((n_1)-3 + (n_2-3)))}.

3) For sample size n, and correlations r12, r13 and r23 test for the difference of two dependent correlations (r12 vs r13).  

4) For sample size n, test for the difference between two dependent correlations involving different variables.

Consider the correlations from Steiger (1980), Table 1:  Because these all from the same subjects, any tests must be of dependent correlations.  For dependent correlations, it is necessary to specify at least 3 correlations (e.g., r12, r13, r23) 

\tabular{lrrrrrrr}{
Variable \tab M1 \tab F1 \tab V1 \tab M2 \tab F2  \tab V2 \cr
M1  1.00 \cr
F1  \tab .10 \tab 1.00\cr
V1 \tab .40 \tab .50 \tab 1.00 \cr
M2 \tab .70 \tab .05 \tab .50 \tab 1.00 \cr
F2 \tab .05 \tab .70 \tab .50 \tab .50 \tab 1.00  \cr
V2 \tab .45 \tab .50 \tab .80 \tab .50 \tab .60 \tab 1.00 \cr
}


For clarity, correlations may be specified by value.  If specified by location and if doing the test of dependent correlations, if three correlations are specified, they are assumed to be in the order r12, r13, r23.

 Consider the examples from Steiger: 
 
 Case A: 
where Masculinity at time 1 (M1) correlates with Verbal Ability .5 (r12), femininity at time 1 (F1) correlates with Verbal ability  r13 =.4, and M1 correlates with F1 (r23= .1).  Then, given the correlations: r12 = .4, r13 = .5, and r23 = .1, t = -.89 for n =103, i.e.,
r.test(n=103, r12=.4, r13=.5,r23=.1)  

Case B: Test whether correlation between two variables (e.g., F and V) is the same over time (e.g. F1V1 = F2V2) 

r.test(n = 103, r12 = 0.5, r34 = 0.6, r23 = 0.5, r13 = 0.7, r14 = 0.5,  r24 = 0.8)



}
\value{
  
  \item{test}{Label of test done}
  \item{z}{z value for tests 2 or 4}
  \item{t}{t value for tests 1 and 3}
  \item{p}{probability value of z or t}

}
\references{

Cohen, J. and Cohen, P. and West, S.G. and Aiken, L.S. (2003) Applied multiple regression/correlation analysis for the behavioral sciences, L.Erlbaum Associates, Mahwah, N.J.

Olkin, I. and Finn, J. D. (1995). Correlations redux. Psychological Bulletin, 118(1):155-164. 

Steiger, J.H. (1980), Tests for comparing elements of a correlation matrix, Psychological Bulletin, 87, 245-251. 

Williams, E.J.  (1959) Regression analysis. Wiley, New York, 1959.


}
\author{William Revelle }
\note{Steiger specifically rejects using the Hotelling T test to test the difference between correlated correlations. Instead, he recommends Williams' test. (See also Dunn and Clark, 1971). These tests follow Steiger's advice.
The test of two independent correlations is just a z test of the difference of the  Fisher's z transformed correlations divided by the standard error of the difference.  (See Cohen et al, p 49).

One of the beautiful features of R is what works on single value works on vectors and matrices.  Thus, r.test can be used to test the pairwise diference of all the elements of a correlation matrix.  See the last example.

By default, the probabilities are reported to 2 decimal places.  This will, of course, sometimes lead to statements such as p < .1  when in fact p < .1001 or even more precisely p < .1000759.  To achieve the higher precision, use a print statement with the preferred number of digits.  See the next to last set of examples (courtesy of Julia Rohrer).
}
\seealso{ See also  \code{\link{corr.test}} which tests all the elements of a  correlation matrix, and \code{\link{cortest.mat}} to compare two matrices of correlations.  r.test  extends the tests in \code{\link{paired.r}},\code{\link{r.con}}}
\examples{

n <- 30 
r <- seq(0,.9,.1) 
rc <- matrix(r.con(r,n),ncol=2) 
test <- r.test(n,r)
r.rc <- data.frame(r=r,z=fisherz(r),lower=rc[,1],upper=rc[,2],t=test$t,p=test$p) 
round(r.rc,2) 

r.test(50,r)
r.test(30,.4,.6)       #test the difference between two independent correlations
r.test(103,.4,.5,.1)   #Steiger case A of dependent correlations 
r.test(n=103, r12=.4, r13=.5,r23=.1)  
#for complicated tests, it is probably better to specify correlations by name
r.test(n=103,r12=.5,r34=.6,r13=.7,r23=.5,r14=.5,r24=.8)   #steiger Case B 


##By default, the precision of p values is 2 decimals
#Consider three different precisions shown by varying the requested number of digits
r12 =  0.693458895410494 
r23 =  0.988475791500198 
r13 =  0.695966022434845 
print(r.test(n =  5105 ,  r12 =  r12 ,  r23 = r23 , r13 =  r13 )) #probability < 0.1
print(r.test(n =  5105 ,  r12 =  r12, r23 = r23 , r13 =  r13 ),digits=4) #p < 0.1001
print(r.test(n =  5105 , r12 =  r12, r23 = r23 , r13 =  r13 ),digits=8) #p< <0.1000759 



#an example of how to compare the elements of two matrices 
R1 <- lowerCor(bfi[1:200,1:5])  #find one set of Correlations
R2 <- lowerCor(bfi[201:400,1:5]) #and now another set sampled from the same population
test <- r.test(n=200, r12 = R1, r34 = R2)
round(lowerUpper(R1,R2,diff=TRUE),digits=2)  #show the differences between correlations
#lowerMat(test$p)  #show the p values of the difference between the two matrices
adjusted <- p.adjust(test$p[upper.tri(test$p)])
both <- test$p 
both[upper.tri(both)] <- adjusted
round(both,digits=2)  #The lower off diagonal are the raw ps, the upper the adjusted ps


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
