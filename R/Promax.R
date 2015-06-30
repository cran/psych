"Promax" <- 
function (x, m = 4) 
{
 if(!is.matrix(x) & !is.data.frame(x) )  {
        if(!is.null(x$loadings)) x <- as.matrix(x$loadings)
      } else {x <- x}   
    if (ncol(x) < 2) 
        return(x)
    dn <- dimnames(x)
   xx <- stats::varimax(x)
    x <- xx$loadings
    Q <- x * abs(x)^(m - 1)
    U <- lm.fit(x, Q)$coefficients
    d <- try(diag(solve(t(U) %*% U)),silent=TRUE)
     if(class(d)=="try-error") {warning("Factors are exactly uncorrelated and the model produces a singular matrix. An approximation is used")
        ev <- eigen(t(U) %*% U)
        ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
        UU <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
       diag(UU)  <- 1
       d <- diag(solve(UU))}
    U <- U %*% diag(sqrt(d))
    dimnames(U) <- NULL
    z <- x %*% U
    U <- xx$rotmat %*% U
    ui <- solve(U)
    Phi <- ui %*% t(ui)
    dimnames(z) <- dn
    class(z) <- "loadings"
    result <- list(loadings = z, rotmat = U,Phi = Phi)
    class(result) <- c("psych","fa")
    return(result)
}
#obviously a direct copy of the promax function, with the addition of returning the angles between factors
#based upon a suggestion to the R-help news group by Ulrich Keller and John Fox. 

#added May 31st following suggestions to R-Help by Gunter Nickel
"equamax" <- function(L, Tmat=diag(ncol(L)),  eps=1e-5, maxit=1000) {
kappa=ncol(L)/(2*nrow(L))
 if(requireNamespace('GPArotation')) {GPArotation::cfT(L, Tmat=diag(ncol(L)),  eps=eps, maxit=maxit)}  else {stop("biquartimin requires GPArotation")}}


#based completely on the GPArotation  GPForth function
#modified to call the varimin function which is derived from the varimax function

varimin <- function(L, Tmat = diag(ncol(L)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000) {
     if(requireNamespace('GPArotation')) {GPArotation::GPForth(A=L,Tmat = diag(ncol(L)), normalize = normalize, eps = eps, 
    maxit = maxit, method = "varimin") }  else {stop("biquartimin requires GPArotation")}}
    
 vgQ.varimin <- 
function (L) 
{
    QL <- sweep(L^2, 2, colMeans(L^2), "-")
    list(Gq = L * QL, f = sqrt(sum(diag(crossprod(QL))))^2/4, 
        Method = "varimin")
}

specialT <- specialQ <- function(L, Tmat = diag(ncol(L)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000) {
     write("A dummy function that can be replaced with either an orthogonal (specialT) or oblique (specialQ) call.  You will need to supply it")
     list(NA)
     }
  
#a general function to call a number of different rotation functions 
#meant to simplify code in fa, principal, faBy, but perhaps not ready for prime time  
#not yet included in the public package 
"faRotate" <-
function(loadings,rotate="oblimin",...)  {
 if (rotate=="varimax" |rotate=="Varimax" | rotate=="quartimax" | rotate =="bentlerT" | rotate =="geominT" | rotate =="targetT" | rotate =="bifactor"   | rotate =="TargetT"|
                       rotate =="equamax"| rotate =="varimin"|rotate =="specialT" | rotate =="Promax"  | rotate =="promax"| rotate =="cluster" |rotate == "biquartimin" |rotate == "TargetQ"  |rotate =="specialQ" ) {
Phi <- NULL 
switch(rotate,  #The orthogonal cases  for GPArotation + ones developed for psych
  varimax = {rotated <- stats::varimax(loadings)  #varimax is from stats, the others are from GPArotation 
   			         loadings <- rotated$loadings},
   Varimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	       #varimax is from the stats package, Varimax is from GPArotations
   			#rotated <- do.call(rotate,list(loadings,...))
   			#rotated <- do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...))
   			rotated <- GPArotation::Varimax(loadings)
   			loadings <- rotated$loadings} ,
   	quartimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	      
   			#rotated <- do.call(rotate,list(loadings))
   			rotated <- GPArotation::quartimax(loadings)
   			loadings <- rotated$loadings} ,
   	bentlerT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	       
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::bentlerT(loadings)
   			loadings <- rotated$loadings} ,
   	geominT	= {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	      
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::geominT(loadings)
   			loadings <- rotated$loadings} ,
   	targetT = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	     
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::targetT(loadings)
   			loadings <- rotated$loadings} ,
   			
   	 bifactor = {loadings <- bifactor(loadings)$loadings},    #the next four solutions were not properly returning the values
   	 TargetT =  {loadings <- TargetT(loadings,...)$loadings},
   	equamax =  {loadings <- equamax(loadings)$loadings}, 
   	varimin = {loadings <- varimin(loadings)$loadings},
   	specialT =  {loadings <- specialT(loadings)$loadings}, 
   	Promax =   {pro <- Promax(loadings)
     			 loadings <- pro$loadings
     			  Phi <- pro$Phi },
     promax =   {pro <- Promax(loadings)
     			 loadings <- pro$loadings
     			  Phi <- pro$Phi },	
     cluster = 	 {loadings <- varimax(loadings)$loadings           			
								pro <- target.rot(loadings)
     			              	loadings <- pro$loadings
     			                Phi <- pro$Phi},
     biquartimin =    {ob <- biquartimin(loadings,)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi}, 
     TargetQ  =  {ob <- TargetQ(loadings,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi}, 
     specialQ = {ob <- specialQ(loadings,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi})
     } else {
     #The following oblique cases all use GPArotation			                
     if (rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax" | rotate =="geominQ"  | rotate =="bentlerQ"  |rotate == "targetQ"  ) {
     				if (!requireNamespace('GPArotation')) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
     				    Phi <- NULL} else { 
     				      
     				             ob <- try(do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...)))
     				               if(class(ob)== as.character("try-error"))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
     				               ob <- Promax(loadings)}
     				                 
     				loadings <- ob$loadings
     				 Phi <- ob$Phi}
     		                             } else {message("Specified rotation not found, rotate='none' used")}
     	 } 
     return(list(loadings=loadings,Phi=Phi))
     }	 
     

