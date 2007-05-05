"VSS.parallel" <-
function(ncases,nvariables)       #function call 
{  
 simdata=matrix(rnorm(ncases*nvariables),nrow=ncases,ncol=nvariables)    #make up simulated data
 testsim <- VSS(simdata,8,"none")
 return(testsim)
 }

