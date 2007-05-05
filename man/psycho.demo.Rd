\name{psycho.demo}
\alias{psycho.demo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create demo data for psychometrics }
\description{A not very interesting demo of what happens if bivariate continuous data are dichotomized.  Bascially a demo of r, phi, and polychor. }
\usage{
psycho.demo()
}
%- maybe also 'usage' for other objects documented here.
\details{Not one of my more interesting demonstrations.  See 
\url{http://personality-project.org/r/simulating-personality.html} and 
\url{http://personality-project.org/r/r.datageneration.html} for better demonstrations of data generation.
}
\value{a matrix of correlations)
  
}
\references{\url{http://personality-project.org/r/simulating-personality.html} and 
\url{http://personality-project.org/r/r.datageneration.html} for better demonstrations of data generation. }
\author{ William Revelle}

\seealso{\code{\link{VSS.simulate}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function() {    
 #simulate correlation matrix with variable cut points -- psychometric demo
 #make up some correlations with different cut points
cuts <-c(-2,-1,0,1,2)
nsub<-1000    #how many subjects
r<-.6         #population correlation of observed with theta
theta <-rnorm(nsub)  #make up some random normal theta scores
err<- rnorm(nsub)    #random normal error scores

obser<- theta*(r) + err*sqrt(1-r*r)  #observed = T + E

#convert to 0/1  with different values of cuts
trunc<- matrix(rep(obser,length(cuts)),ncol=length(cuts))  
for (i in 1:length(cuts)) {
   trunc[obser>cuts[i],i]<- 1
   trunc[obser< cuts[i],i]<- 0}
   
d.mat<- data.frame(theta,obser,trunc)  #combine into a data frame
trunc.cor<- cor(d.mat)                 #find the correlations
freq <- mean(d.mat)                    #find the frequencies of scores

#now, convert the upper diagonal to polychorics using John Fox's polychor and my phi2poly

for (i in 4:length(d.mat)) {
   for (j in 3:i) {
       trunc.cor[j,i]<- phi2poly(trunc.cor[i,j],freq[i],freq[j]) 
       }}

  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
