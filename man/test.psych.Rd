\name{test.psych}
\alias{test.psych}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Testing of functions in the psych package }
\description{Test to make sure the psych functions run on basic test data sets
}
\usage{
test.psych(first=1,last=5,short=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{first}{first=1: start with dataset first}
  \item{last}{last=5: test for datasets until last}
  \item{short}{short=TRUE - don't return any analyses}
}
\details{When modifying the psych package, it is useful to make sure that adding some code does not break something else.  This function tests the major functions on various standard data sets.
}
\value{
  \item{out }{if short=FALSE, then list of the output from all functions tested}

}

\author{ William Revelle}
\note{ Not for general consumption -- used to make sure functions throw error messages or correct for weird conditions.

The datasets tested are part of the standard R data sets and represent some of the basic problems encountered.}

\section{Warning }{Warning messages will be thrown by fa.parallel} 
}

\examples{
#test <- psych.test()
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
test.psych <- 
function(first=1,last=5,short=TRUE) {  
s1 <- USArrests         #  Violent Crime Rates by US State  (4 variables)
s2 <- attitude          #The Chatterjee-Price Attitude Data
s3 <- Harman23.cor$cov     #   Harman Example 2.3 8 physical measurements
s4 <- Harman74.cor$cov     #   Harman Example 7.4  24 mental measurements
s5 <- ability.cov$cov       #  8 Ability and Intelligence Tests 

#convert covariance to correlation
d5 <- diag(1/sqrt(diag(s5)))
s5 <- d5 %*% s5 %*% d5

datasets <- list(s1,s2,s3,s4,s5)
out <- list()

for (i in first:last) {
   test.data <- datasets[[i]]
 	pc <-   principal(test.data)
	pc2 <-    principal(test.data,2)
	fa2 <- factor.pa(test.data,2)
	fp <-    fa.parallel(test.data)
	ic <-   ICLUST(test.data)
	om <-  omega(test.data)
	fc <- factor.congruence(pc2,fa2)
	vss2 <- VSS(test.data)
	vsspc <- VSS(test.data,pc="pc")
	vss.scree <- VSS.scree(test.data)
	d <- describe(test.data)
	
	
	keys <- matrix(rep(0,dim(test.data)[2]*2),ncol=2)
	keys[,1] <- 1
	keys[1:3,2] <- 1
	if( dim(test.data)[1] != dim(test.data)[2]) {test.score <- score.items(keys,test.data)} else {test.score <- cluster.cor(keys,test.data)}
	
	out <- list(out,paste("test",i),pc,pc2,fa2,fp,ic,om,fc,vss2,vsspc,vss.scree,d,test.score)
  } #end loop
  
  #a few more tests
  
  simple <- item.sim(nvar=24)
  circ <-  circ.sim(nvar=24)
  fa.simple <- factor.pa(simple)
  psych.d <- psycho.demo()
  
 if (!short) { return(out)}

}#end function


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

