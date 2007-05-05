"phi2poly" <-
function(ph,cp,cc) {
	require(polycor)
     #ph is the phi coefficient
     #cp is the selection ratio of the predictor
     #cc is the success rate of the criterion
     r.marg<-rep(0,2)
     c.marg<- rep(0,2)
     p<-array(rep(0,4),dim=c(2,2))
     r.marg[1]<- cp
     r.marg[2]<- 1 -cp 
     c.marg[1]<- cc
     c.marg[2]<- 1-cc
     
     p[1,1]<- r.marg[1]*c.marg[1]+ ph*sqrt(prod(r.marg,c.marg))
     p[2,2]<- r.marg[2]*c.marg[2]+ ph*sqrt(prod(r.marg,c.marg))
     p[1,2]<- r.marg[1]*c.marg[2]- ph*sqrt(prod(r.marg,c.marg))
     p[2,1]<- r.marg[2]*c.marg[1]- ph*sqrt(prod(r.marg,c.marg))
     
     result<-polychor(p ) 
     return(result)}

