#reorganized May 25, 2009 to call several print functions (psych.print.x where x = {fa, omega, stats, vss}
#reorganized, January 18, 2009 to make some what clearer

"print.psych" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

iclust <- omega <- vss <- scores <- fac.pa <- principal <- gutt <- sim <- alpha <- describe <- corr.test <- r.test <- cortest <-  cluster.cor <- cluster.loadings <- thurstone <- stats <- ICC <- mat.reg <-  FALSE
#first, figure out which psych function was called
if(length(class(x)) > 1)  {
   if(class(x)[2] =='sim')  sim <- TRUE
   if(class(x)[2] =='vss')  vss <- TRUE
   if(class(x)[2] =='iclust')  iclust <- TRUE
   if(class(x)[2] =='omega')  omega <- TRUE
   if(class(x)[2] =='fa')  fac.pa <- TRUE
   if(class(x)[2] =='principal') principal <- fac.pa <- TRUE
   if(class(x)[2] == 'alpha') alpha <- TRUE
   if(class(x)[2] == 'describe') describe <- TRUE
   if(class(x)[2] == "corr.test") corr.test <- TRUE
   if(class(x)[2] == "r.test") r.test <- TRUE
   if(class(x)[2] == "cortest") cortest <- TRUE
   if(class(x)[2] == "score.items") scores <- TRUE
   if(class(x)[2] == "cluster.cor") cluster.cor <- TRUE
   if(class(x)[2] == "cluster.loadings") cluster.loadings <- TRUE
   if(class(x)[2] == "guttman") gutt <- TRUE
   if(class(x)[2] == "thurstone") thurstone <- TRUE
   if(class(x)[2] == "ICC") ICC <- TRUE
   if(class(x)[2] == "stats") stats <- TRUE
   if(class(x)[2] == "mat.regress") mat.reg <- TRUE
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
 if(vss) print.psych.vss(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(fac.pa) print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(iclust) print.psych.iclust(x,digits=digits,all=all,cut=cut,sort=sort,...)
 if(stats) print.psych.stats(x,digits=digits,all=all,cut=cut,sort=sort,...)

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
                                
 if(scores) {
    cat("Call: ")
    print(x$Call)
	cat("\n(Unstandardized) Alpha:\n")
	print(x$alpha,digits=digits)
  	cat("\nAverage item correlation:\n")
  	print(x$av.r,digits=digits)
	cat("\n Guttman 6* reliability: \n")
	print(x$G6,digits=digits)     
	if(iclust) {cat("\nOriginal Beta:\n")
	 print(x$beta,digits) }	          
 	
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	 
	  if(!is.null(x$item.cor) ) {
	   cat("\nItem by scale correlations:\n corrected for item overlap and scale reliability\n" )
	 print(round(x$item.corrected,digits=digits)) } 
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
}


####
 if(gutt) {
  cat("Call: ")
    print(x$Call)
 cat ("\nAlternative estimates of reliability\n")
 cat("Beta = ", x$beta, " This is an estimate of the worst split half reliability")  
 cat ("\nGuttman bounds \nL1 = ",x$lambda.1, "\nL2 = ", x$lambda.2, "\nL3 (alpha) = ", x$lambda.3,"\nL4 (max) = " ,x$lambda.4, "\nL5 = ", x$lambda.5, "\nL6 (smc) = " ,x$lambda.6, "\n")
 cat("TenBerge bounds \nmu0 = ",x$tenberge$mu.0, "mu1 = ", x$tenberge$mu1, "mu2 = " ,x$tenberge$mu2, "mu3 = ",x$tenberge$mu3 , "\n")
 cat("\nalpha of first PC = ", x$alpha.pc, "\nestimated greatest lower bound = ", x$lambda.4,"\n")
 cat("\nbeta estimated by first and second PC = ", round(x$beta.pc,digits), " This is an exploratory statistic \n")
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
            print(x$results)
            cat("\n Number of subjects =", x$n.obs, "    Number of Judges = ",x$n.judge)

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
           print(signif(x$probF,digits))
            cat("\n degrees of freedom of regression \n") 
            print(x$df)
           }

   }
  
  }     #end of the not list condition
  
}  #end function
