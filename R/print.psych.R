#reorganized May 25, 2009 to call several print functions (psych.print.x where x = {fa, omega, stats, vss}
#reorganized, January 18, 2009 to make some what clearer

"print.psych" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

iclust <- omega <- omegaSem <- vss <- scores <- mchoice <- fac.pa <- principal <- gutt <- sim <- alpha <- describe <- corr.test <- r.test <- cortest <-  cluster.cor <- cluster.loadings <-comorbid <- kappa <- thurstone <- stats <- ICC <- mat.reg <- parallel <- tetra <- poly <-  irt.fa <- irt.poly <- mardia <- partial <- extension <- circ  <- schmid <-  FALSE
#first, figure out which psych function was called
if(length(class(x)) > 1)  {
   if(class(x)[2] =='sim')  sim <- TRUE
   if(class(x)[2] =='vss')  vss <- TRUE
   if(class(x)[2] =='iclust')  iclust <- TRUE
   if(class(x)[2] =='omega')  omega <- TRUE
    if(class(x)[2] =='omegaSem')  omegaSem <- TRUE
   if(class(x)[2] =='fa')  fac.pa <- TRUE
  if(class(x)[2] == "schmid")         schmid <- TRUE
   if(class(x)[2] =='principal') principal <- fac.pa <- TRUE
   if(class(x)[2] == 'alpha') alpha <- TRUE
   if(class(x)[2] == 'describe') describe <- TRUE
   if(class(x)[2] == "corr.test") corr.test <- TRUE
   if(class(x)[2] == "r.test") r.test <- TRUE
   if(class(x)[2] == "cortest") cortest <- TRUE
   if(class(x)[2] == "score.items") scores <- TRUE
    if(class(x)[2] == "mchoice") mchoice <- TRUE
   if(class(x)[2] == "cluster.cor") cluster.cor <- TRUE
   if(class(x)[2] == "cluster.loadings") cluster.loadings <- TRUE
   if(class(x)[2] == "guttman") gutt <- TRUE
   if(class(x)[2] == "thurstone") thurstone <- TRUE
   if(class(x)[2] == "ICC") ICC <- TRUE
   if(class(x)[2] == "stats") stats <- TRUE
   if(class(x)[2] == "mat.regress") mat.reg <- TRUE
   if(class(x)[2] == "parallel") parallel <- TRUE
   if(class(x)[2] == "kappa")    kappa <- TRUE 
   if(class(x)[2] == "comorbid")    comorbid <- TRUE
   if(class(x)[2] == "tetra")    tetra <- TRUE
   if(class(x)[2] == "poly")    poly <- TRUE
   if(class(x)[2] == "irt.fa")    irt.fa <- TRUE
   if(class(x)[2] == "irt.poly")    irt.poly <- TRUE
   if(class(x)[2] == "mardia")    mardia <- TRUE
  if(class(x)[2] == "partial.r")    partial <- TRUE
  if(class(x)[2] == "extension")    extension <- TRUE
  if(class(x)[2] == "circ")         circ <- TRUE
     } 
else {     
#these next test for non-psych functions that may be printed using print.psych.fa
if(!is.null(x$communality.iterations)) {fac.pa <- TRUE}
if(!is.null(x$uniquenesses)) {fac.pa <- TRUE}
if(!is.null(x$rotmat)) {fac.pa <- TRUE}
if(!is.null(x$Th)) {fac.pa <- TRUE}
}  

## the following functions have their own print function
 if(omega)  print.psych.omega(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(omegaSem)print.psych.omegaSem(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(schmid) print.psych.schmid(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(vss) print.psych.vss(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(fac.pa) print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(iclust) print.psych.iclust(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(stats) print.psych.stats(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(extension) print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)

## 
##Now, for the smaller print jobs, just do it here.
if(all) {class(x) <- "list"
         print(x) } else {    #find out which function created the data and then print accordingly


if(describe) {if  (length(dim(x))==1) {class(x) <- "list"
              attr(x,"call") <- NULL
              print(round(x,digits=digits))
                  } else  {class(x) <- "data.frame" 
            print(round(x,digits=digits)) }
         }
         
          
if(corr.test) {cat("Call:")
              print(x$Call)
              cat("Correlation matrix \n")
              print(round(x$r,digits))
              cat("Sample Size \n")
              print(x$n)
              cat("Probability value \n")
             print(round(x$p,digits))
         }

if(r.test) {cat("Correlation tests \n")
            cat("Call:")
              print(x$Call)
              cat( x$Test,"\n")
              if(!is.null(x$t)) {cat(" t value" ,round(x$t,digits),"   with probability <", signif(x$p,digits) )}
              if(!is.null(x$z)) {cat(" z value" ,round(x$z,digits),"   with probability ",  round(x$p,digits) )}               
              if(!is.null(x$ci)) {cat("\n and confidence interval ",round(x$ci,digits) ) }
         }
         
 if(cortest) {cat("Tests of correlation matrices \n")
            cat("Call:")
            print(x$Call)
            cat(" Chi Square value" ,round(x$chi,digits)," with df = ",x$df, "  with probability <", signif(x$p,digits) )
         }

if(partial) {cat("partial correlations \n")
            print(round(unclass(x),digits))
         }
         
if(circ)    {cat("Tests of circumplex structure \n")
            cat("Call:")
             print(x$Call)
             res <- data.frame(x[1:4]) 
             print(res,digits=2)
             }

             
 if(scores) {
    cat("Call: ")
    print(x$Call)
    if(x$raw) {
	cat("\n(Unstandardized) Alpha:\n") } else {cat("\n(Standardized) Alpha:\n") }
	
	print(x$alpha,digits=digits)
	if(!is.null(x$alpha.ob)) {cat("\nStandardized Alpha of observed scales:\n")
	print(x$alpha.ob,digits=digits)}
  	cat("\nAverage item correlation:\n")
  	print(x$av.r,digits=digits)
	cat("\n Guttman 6* reliability: \n")
	print(x$G6,digits=digits)     
	if(iclust) {cat("\nOriginal Beta:\n")
	 print(x$beta,digits) }	          
 	
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 	if(!is.null(x$alpha.ob)) {cat("\nNote that these are the correlations of the complete scales based on the correlation matrix,\n not the observed scales based on the raw items.\n")}
	 	
	 print(x$corrected,digits) 
	 
	  if(!is.null(x$item.cor) ) {
	   cat("\nItem by scale correlations:\n corrected for item overlap and scale reliability\n" )
	   
	 print(round(x$item.corrected,digits=digits)) } 
	 if(!is.null(x$response.freq)) {
	 cat("\nNon missing response frequency for each item\n")
	 print(round(x$response.freq,digit=digits))}
	
  }
  
   if(mchoice) {
    cat("Call: ")
    print(x$Call)
	cat("\n(Unstandardized) Alpha:\n")
	print(x$alpha,digits=digits)
  	cat("\nAverage item correlation:\n")
  	print(x$av.r,digits=digits)
	
	
	
	 if(!is.null(x$item.stats)) {
	 cat("\nitem statistics \n")
	 print(round(x$item.stats,digit=digits))}
  }
  
  
if(cluster.cor) {
    cat("Call: ")
    print(x$Call)
	cat("\n(Standardized) Alpha:\n")
	print(x$alpha,digits)
    cat("\n(Standardized) G6*:\n")
    print(x$G6,digits)
  	cat("\nAverage item correlation:\n")
	print(x$av.r,digits)
	cat("\nNumber of items:\n")
	print(x$size)
	
	      
 	 # cat("\nScale intercorrelations:\n")
	 # print(x$cor,digits=digits)
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	  }
	
##	
if(cluster.loadings) {
    cat("Call: ")
    print(x$Call)
	cat("\n(Standardized) Alpha:\n")
	print(x$alpha,digits)
		cat("\n(Standardized) G6*:\n")
	print(x$G6,digits)
  	cat("\nAverage item correlation:\n")
	print(x$av.r,digits)
	cat("\nNumber of items:\n")
	print(x$size)
    cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	 
   cat("\nItem by scale intercorrelations\n corrected for item overlap and scale reliability\n")
	 print(x$loadings,digits) 
  #cat("\nItem by scale Pattern matrix\n")
	# print(x$pattern,digits) 
	
  }
####
if(alpha) {
cat("\nReliability analysis ",x$title," \n")
cat("Call: ")
print(x$call)
cat("\n ")
print(x$total,digits=digits)
cat("\n Reliability if an item is dropped:\n")
    print(x$alpha.drop,digits=digits)
cat("\n Item statistics \n")
   print(x$item.stats,digits=digits)
 if(!is.null(x$response.freq)) {
	 cat("\nNon missing response frequency for each item\n")
	 print(round(x$response.freq,digit=digits))}
}


####
 if(gutt) {
  cat("Call: ")
    print(x$Call)
 cat ("\nAlternative estimates of reliability\n")
 cat("Beta = ", round(x$beta,digits), " This is an estimate of the worst split half reliability")  
 cat ("\nGuttman bounds \nL1 = ",round(x$lambda.1,digits), "\nL2 = ", round(x$lambda.2,digits), "\nL3 (alpha) = ", round(x$lambda.3,digits),"\nL4 (max) = " ,round(x$lambda.4,digits), "\nL5 = ", round(x$lambda.5,digits), "\nL6 (smc) = " ,round(x$lambda.6,digits), "\n")
 cat("TenBerge bounds \nmu0 = ",round(x$tenberge$mu0,digits), "mu1 = ", round(x$tenberge$mu1,digits), "mu2 = " ,round(x$tenberge$mu2,digits), "mu3 = ",round(x$tenberge$mu3,digits) , "\n")
 cat("\nalpha of first PC = ",round( x$alpha.pc,digits), "\nestimated greatest lower bound based upon communalities= ", round(x$glb,digits),"\n")
 #cat("\nbeta estimated by first and second PC = ", round(x$beta.pc,digits), " This is an exploratory statistic \n")
 } 
  
 ## 
 if(thurstone) {
cat("Thurstonian scale (case 5) scale values ")
cat("\nCall: ")
     print(x$Call)
print(x$scale)
cat("\n Goodness of fit of model  ", round(x$GF,digits))
 }
 ##
 if(sim) { if(is.matrix(x)) {x <-unclass(x) 
             round(x,digits) } else {
              cat("Call: ")
              print(x$Call)
             cat("\n $model (Population correlation matrix) \n")
             print(x$model,digits)
             if(!is.null(x$reliability)) { cat("\n$reliability (population reliability) \n")
                print(x$reliability,digits) } 
             if(!is.null(x$N) && !is.null(x$r)) {
             cat("\n$r  (Sample correlation matrix  for sample size = ",x$N,")\n")
             print(x$r,digits)}
             }
             
 }
            
             
 if(ICC) {cat("Call: ")
              print(x$Call)
            cat("\nIntraclass correlation coefficients \n")
            print(x$results,digits=digits)
            cat("\n Number of subjects =", x$n.obs, "    Number of Judges = ",x$n.judge)

   }
  
   if(kappa) {if(is.null(x$cohen.kappa)) {
            cat("Call: ")
            print(x$Call)
            
            cat("\nCohen Kappa and Weighted Kappa correlation coefficients and confidence boundaries \n")
           
            print(x$confid,digits=digits)
            cat("\n Number of subjects =", x$n.obs,"\n")} else {
            cat("\nCohen Kappa (below the diagonal) and Weighted Kappa (above the diagonal) \nFor confidence intervals and detail print with all=TRUE\n")
            print(x$cohen.kappa,digits) 
            }
   }
   
  if(comorbid) {cat("Call: ")
              print(x$Call)
            cat("Comorbidity table \n")
            print(x$twobytwo,digits=digits)
            cat("\nimplies phi = ",round(x$phi,digits), " with Yule = ", round(x$Yule,digits), " and tetrachoric correlation of ", round(x$tetra$rho,digits))
            cat("\nand normal thresholds of ",round(-x$tetra$tau,digits))
          
   }
   
  if(tetra) {cat("Call: ")
              print(x$Call)
            cat("tetrachoric correlation \n")
            if(!is.null(x$twobytwo)) {
              print(x$twobytwo,digits=digits)
              cat("\n implies tetrachoric correlation of ",round(x$rho,digits))} else {
            
            print(x$rho,digits)
            cat("\n with tau of \n")
            print(x$tau,digits)
          }
   }
  
  
    if(poly) {cat("Call: ")
              print(x$Call)
            cat("Polychoric correlations \n")
            if(!is.null(x$twobytwo)) {
              print(x$twobytwo,digits=digits)
              cat("\n implies tetrachoric correlation of ",round(-x$rho,digits))} else {
            
            print(x$rho,digits)
            cat("\n with tau of \n")
            print(x$tau,digits)
          }
   }
  if(mat.reg) { cat("Call: ")
              print(x$Call)
            cat("\nMultiple Regression from matrix input \n")
           cat("\nBeta weights \n")
           print(round(x$beta,digits))
           cat("\nMultiple R \n") 
           print(round(x$R,digits))
            cat("\nMultiple R2 \n") 
            print(x$R2,digits)
              if(!is.null(x$se)) {
               cat("\n SE of Beta weights \n")
           print(round(x$se,digits))
           cat("\n t of Beta Weights \n") 
           print(round(x$t,digits))
            cat("\nProbability of t < \n") 
            print(signif(x$Probability,digits))
            cat("\n Shrunken R2 \n")
           print(x$shrunkenR2,digits)
           cat("\nStandard Error of R2  \n") 
           print(x$seR2,digits)
            cat("\nF \n") 
            print(round(x$F,digits))
             cat("\nProbability of F < \n") 
           print(signif(x$probF,digits+1))
            cat("\n degrees of freedom of regression \n") 
            print(x$df)
           }

   }
 if(parallel) { cat("Call: ")
              print(x$Call) 
              if(!is.null(x$fa.values) & !is.null(x$pc.values) ) {
                  parallel.df <- data.frame(fa=x$fa.values,fa.sim =x$fa.sim$mean,pc= x$pc.values,pc.sim =x$pc.sim$mean)
fa.test <- which(!(x$fa.values > x$fa.sim$mean))[1]-1
pc.test <- which(!(x$pc.values > x$pc.sim$mean))[1] -1
cat("Parallel analysis suggests that ")
cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
                  cat("\n Eigen Values of \n")
                  
                  colnames(parallel.df) <- c("Original factors","Simulated data","Original components", "simulated data")
                  print(round(parallel.df[1:max(fa.test,pc.test),],digits))} else {
              if(!is.null(x$fa.values)) {cat("\n eigen values of factors\n")
              print(round(x$fa.values,digits))}
               if(!is.null(x$fa.sim)){cat("\n eigen values of simulated factors\n") 
               print(round(x$fa.sim$mean,digits))}
              if(!is.null(x$pc.values)){cat("\n eigen values of components \n")
              print(round(x$pc.values,digits))}
             
              if(!is.null(x$pc.sim)) {cat("\n eigen values of simulated components\n") 
              print(round(x$pc.sim$mean,digits))}
              
              }
         }
  
  if(irt.fa) {
   cat("Item Response Analysis using Factor Analysis = ")
   cat("\nCall: ")
   print(x$Call)
    nf <- length(x$irt$difficulty)
    for(i in 1:nf) {temp <- data.frame(discrimination=x$irt$discrimination[,i],location=x$irt$difficulty[[i]])
    cat("\nItem discrimination and location for factor ",colnames(x$irt$discrimination)[i],"\n")
    print(round(temp,digits))}
 # print(round(x$coefficients,digits))
  }
  
    if(irt.poly) {
   cat("Item Response Analysis using Factor Analysis = ")
   cat("\nCall: ")
   print(x$Call)
    
    nf <- length(x$irt$difficulty)
    for(i in 1:nf) {temp <- data.frame(discrimination=x$irt$discrimination[,i],location=x$irt$difficulty[[i]])
    cat("\nItem discrimination and location for factor ",colnames(x$irt$discrimination)[i],"\n")
    print(round(temp,digits))}
  #print(round(x$coefficients,digits))
  }
  
  if(mardia) {
   cat("Call: ")
     print(x$Call) 
     cat("\nMardia tests of multivariate skew and kurtosis\n")
     cat("Use describe(x) the to get univariate tests")
      cat("\nn.obs =",x$n.obs,"  num.vars = ",x$n.var,"\n")
     cat("b1p = ",round(x$b1p,digits),"  skew = ",round(x$skew,digits ), " with probability = ", signif(x$p.skew,digits)) 
     cat("\n small sample skew = ",round(x$small.skew,digits ), " with probability = ", signif(x$p.small,digits)) 
     cat("\nb2p = ", round(x$b2p,digits),"  kurtosis = ",round(x$kurtosis,digits)," with probability = ",signif(x$p.kurt,digits ))
  }
  
  
  }     #end of the not list condition
  
}  #end function
