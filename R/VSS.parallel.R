"VSS.parallel" <-
function(ncases,nvariables,scree=FALSE,rotate="none")       #function call 
{  
 simdata=matrix(rnorm(ncases*nvariables),nrow=ncases,ncol=nvariables)    #make up simulated data
 if(scree) {VSS.scree(simdata) 
 testsim <- simdata}  else {testsim <- VSS(simdata,8,rotate)}
 return(testsim)
 }

