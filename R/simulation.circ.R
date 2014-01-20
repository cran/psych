"simulation.circ" <-
function(samplesize=c(100,200,400,800), numberofvariables=c(16,32,48,72))  {
 ncases=length(samplesize)
 nvar <- length(numberofvariables)
  results <- matrix(NaN,ncol=ncases,nrow=nvar*ncases)
  results.ls <- list()
   case <- 1
  for (ss in 1:ncases) {
    for (nv in 1:nvar) {
    
      circ.data <- circ.sim(nvar=numberofvariables[nv],nsub=samplesize[ss])
      sim.data <-  circ.sim(nvar=numberofvariables[nv],nsub=samplesize[ss],circum=FALSE)
      elipse.data <- circ.sim(nvar=numberofvariables[nv],nsub=samplesize[ss],yloading=.4)
      r.circ<- cor(circ.data)
      r.sim <- cor(sim.data)
      r.elipse <- cor(elipse.data)
      pc.circ <- principal(r.circ,2)
      pc.sim <- principal(r.sim,2)
      pc.elipse <- principal(r.elipse,2)
      case <- case + 1
      results.ls[[case]] <- list(numberofvariables[nv],samplesize[ss],circ.tests(pc.circ)[1:4],circ.tests(pc.elipse)[1:4],circ.tests(pc.sim)[1:4])
     }
     }
     results.mat <- matrix(unlist(results.ls),ncol=14,byrow=TRUE)
    colnames(results.mat) <- c("nvar","n","c-gap","c-fisher","c-RT","c-VT","e-gap","e-fisher","e-RT","e-VT","s-gap","s-fisher","s-RT","s-VT")
    results.df <- data.frame(results.mat)
 	 return(results.df)
  }