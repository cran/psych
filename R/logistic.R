"logistic" <-  
function(x,d=0, a=1,c=0,z=1) {c + (z-c)/(1+exp(a*(d-x)))}

"logit" <- 
function(p) {log(p/(1-p))}
#created March 20, 2010

#graded response model
"logistic.grm" <- 
function(x,d=0,a=1.5,c=0,z=1,r=2,s=c(-1.5,-.5,.5,1.5)){
if (r == 1) {p <- (1-logistic(x,d=s[1],a=a,c=c,z=z))} else {
if (r == (length(s)+1)) {p <- logistic(x,d=s[r-1],a=a,c=c,z=z) } else {
p <-  logistic(x,d=s[r-1],a=a,c=c,z=z) - logistic(x,d=s[r],a=a,c=c,z=z )
}}
p}