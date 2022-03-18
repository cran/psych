faRotations <- function(loadings,r=NULL, rotate="oblimin" , hyper= .15, n.rotations =10,...) {

#adapted the multiple rotation technique from Niels Waller in fungible  faMain

#from fungible
 randStart <- function(dimension) {
        qr.Q(qr(matrix(rnorm(dimension^2), dimension, dimension)))}
        
 starts <- vector("list", n.rotations)  
  hyperplane <- fit <- complexity <- indetermin <-  list()
  fn <- NULL
  if(inherits(loadings, "fa")) {fn <- loadings$fn
  loadings <- loadings$loadings } 
  if(inherits(loadings, "principal")) {fn <- loadings$fn
      loadings <- loadings$loadings} 
  
  
 starts <- lapply(starts, function(x) randStart(dimension = ncol(loadings)))
 starts[[1]] <- diag(ncol(loadings))
original <- loadings
if(!is.null(r)) {r.inv <- solve(r)} else {r.inv <- NULL}
 for(i in 1:n.rotations) {
 initial <- starts[[i]]
 loadings <- original
 
  cl <- match.call()
 
 gpa <- c("Varimax","quartimax","bentlerQ","bentlerT","geominT","geominQ","targetT","targetQ","TargetT","quartimin","oblimin","simplemax") 
if (rotate=="varimax" |rotate=="Varimax" | rotate=="quartimax" | rotate =="bentlerT" | rotate =="geominT" | rotate =="targetT" | rotate =="bifactor"   | rotate =="TargetT"|
                       rotate =="equamax"| rotate =="varimin"|rotate =="specialT" | rotate =="Promax"  | rotate =="promax"| rotate =="cluster" |rotate == "biquartimin" |rotate == "TargetQ"  |rotate =="specialQ" ) {
Phi <- NULL 
switch(rotate,  #The orthogonal cases  for GPArotation + ones developed for psych
  varimax = {rotated <- stats::varimax(loadings)  #varimax is from stats, the others are from GPArotation 
   			         loadings <- rotated$loadings
   			         rot.mat <- rotated$rotmat},
   Varimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	       #varimax is from the stats package, Varimax is from GPArotations
   			#rotated <- do.call(rotate,list(loadings,...))
   			#rotated <- do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...))
   			rotated <- GPArotation::Varimax(loadings,Tmat=initial,...)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	quartimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	      
   			#rotated <- do.call(rotate,list(loadings))
   			rotated <- GPArotation::quartimax(loadings,Tmat=initial,...)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	bentlerT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	       
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::bentlerT(loadings,Tmat=initial,...)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	geominT	= {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	      
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::geominT(loadings,Tmat=initial,...)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	targetT = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   			rotated <- GPArotation::targetT(loadings,Tmat=initial,...)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   			
   	 bifactor = {rot <- bifactor(loadings,Tmat=initial)
   	             loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},  
   	 TargetT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	            rot <- GPArotation::targetT(loadings,Tmat=initial,...)
   	              loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},
   	equamax =  {rot <- equamax(loadings,Tmat=initial,...)
   	              loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))}, 
   	varimin = {rot <- varimin(loadings,Tmat=initial,...)
   	            loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},
   	specialT =  {rot <- specialT(loadings,Tmat=initial,...)
   	              loadings <- rot$loadings
   	              rot.mat <- t(solve(rot$Th))}, 
   	Promax =   {pro <- Promax(loadings)  #Promax without Kaiser normalization
     			loadings <- pro$loadings
     			 Phi <- pro$Phi 
     			 rot.mat <- pro$rotmat},
     promax =   {#pro <- stats::promax(loadings,...)   #from stats
                pro <- kaiser(loadings,rotate="Promax",...)   #calling promax will now do the Kaiser normalization before doing Promax rotation
     			 loadings <- pro$loadings
     			  rot.mat <- pro$rotmat
     			 # ui <- solve(rot.mat)
     			 # Phi <-  cov2cor(ui %*% t(ui))
     			  Phi <- pro$Phi 
     			},	
     cluster = 	 {loadings <- varimax(loadings,...)$loadings           			
								pro <- target.rot(loadings)
     			              	loadings <- pro$loadings
     			                Phi <- pro$Phi
     			                 rot.mat <- pro$rotmat},
     biquartimin =    {ob <- biquartimin(loadings,Tmat=initial,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				 rot.mat <- t(solve(ob$Th))}, 
     TargetQ  =  {ob <- TargetQ(loadings,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}, 
     specialQ = {ob <- specialQ(loadings,Tmat=initial,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				 rot.mat <- t(solve(pro$Th))})
     } else {
     

     #The following oblique cases all use GPArotation			                
     if (rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax" | rotate =="geominQ"  | rotate =="bentlerQ"  |rotate == "targetQ"  ) {
     				if (!requireNamespace('GPArotation')) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
     				    Phi <- NULL} else { 
     				      
     				             ob <- try(do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,Tmat=initial)))
     				               if(inherits(ob,as.character("try-error")))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
     				               ob <- Promax(loadings)}
     				                 
     				loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}
     		                             } else {message("Specified rotation not found, rotate='none' used")
     		                             Phi<- NULL}
     	 }
     	
     	 		
    signed <- sign(colSums(loadings))
    signed[signed==0] <- 1
    loadings <- loadings %*% diag(signed)  #flips factors to be in positive direction but loses the colnames
    if(!is.null(Phi)) {Phi <- diag(signed) %*% Phi %*% diag(signed) }  #added October 20, 2009 to correct bug found by Erich Studerus
 
 #fit statistics added 3/14/22
 
 if(!is.null(r.inv)) {  

           if(!is.null(Phi)) {indetermin[[i]] <- tr(t(loadings) %*% r.inv %*% loadings %*% Phi) } else {
                 indetermin[[i]] <- tr(t(loadings) %*% r.inv %*% loadings)
                  }
          }  else {indetermin[[i]] <- NA}
 hyperplane[[i]] <- mean(abs(loadings) < hyper)
complexity[[i]] <- mean((apply(loadings,1,function(x) sum(x^2)))^2/apply(loadings,1,function(x)sum(x^4)))
 fit[[i]]<- tr( var(loadings^2))
}


hyperplane<- unlist(hyperplane)
fit <- unlist(fit)
complexity<- unlist(complexity)
indetermin <- unlist(indetermin)
stats <- matrix(c(hyperplane,fit,complexity,indetermin),ncol=4)

colnames(stats)<- c("hyperplane","fit","complexity","indetermin")

best <- which(stats[,"hyperplane"]==max(stats[,"hyperplane"]))

if(length(best) > 1) best <- which(stats[best,"complexity"]==min(stats[best,"complexity"]))
  if(length(best)> 1) which(stats[best,"fit"]==min(stats[best,"fit"]))
if(length(best)> 1) best <- best[1]

if(rotate =="varimax") {rotate <- "Varimax"
    warning("Varimax from GPA rotation is used instead of varimax from stats")}

Tmat <- starts[[best]]
# now rerotate to the  best fitted solution
# need to consider those rotations that do not use GPArotation
 gpa <- c("Varimax","quartimax","bentlerQ","bentlerT","geominT","geominQ","targetT","targetQ","TargetT","quartimin","oblimin","simplimax")
 if(rotate %in% gpa) {

ob <- try(do.call(getFromNamespace(rotate,'GPArotation'),list(original,Tmat=Tmat)))
  if(inherits(ob,as.character("try-error")))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation") 
    rot.mat<- NULL} else {
loadings <- ob$loadings
Phi <- ob$Phi
 rot.mat <- t(solve(ob$Th))}
 } else {
 #do the other rotations
 switch(rotate, 
   	 bifactor = {rot <- bifactor(loadings,Tmat=Tmat)
   	             loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},
   	  	equamax =  {rot <- equamax(loadings,Tmat=Tmat)
   	              loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))}, 
   	varimin = {rot <- varimin(loadings,Tmat=Tmat)
   	            loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},
   	specialT =  {rot <- specialT(loadings,Tmat=Tmat)
   	              loadings <- rot$loadings
   	              rot.mat <- t(solve(rot$Th))}, 
   	Promax =   {pro <- Promax(loadings)  #Promax without Kaiser normalization
     			loadings <- pro$loadings
     			 Phi <- pro$Phi 
     			 rot.mat <- pro$rotmat},
     promax =   {#pro <- stats::promax(loadings,...)   #from stats
                pro <- kaiser(loadings,rotate="Promax",...)   #calling promax will now do the Kaiser normalization before doing Promax rotation
     			 loadings <- pro$loadings
     			  rot.mat <- pro$rotmat
     			 # ui <- solve(rot.mat)
     			 # Phi <-  cov2cor(ui %*% t(ui))
     			  Phi <- pro$Phi 
     			},	
     cluster = 	 {loadings <- varimax(loadings,...)$loadings           			
								pro <- target.rot(loadings)
     			              	loadings <- pro$loadings
     			                Phi <- pro$Phi
     			                 rot.mat <- pro$rotmat},
     biquartimin =    {ob <- biquartimin(loadings,Tmat=Tmat)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				 rot.mat <- t(solve(ob$Th))}, 
     TargetQ  =  {ob <- TargetQ(loadings,Tmat=Tmat,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}, 
     specialQ = {ob <- specialQ(loadings,Tmat=Tmat)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				 rot.mat <- t(solve(pro$Th))}) 
 }
    
 result <- list(rotation.stats=stats,loadings=loadings,Phi=Phi,rot.mat=rot.mat, fn=fn, Call=cl)
 class(result)<- c("psych","fa")
 return(result)
  }
  
  
  