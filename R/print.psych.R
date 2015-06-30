#reorganized May 25, 2009 to call several print functions (psych.print.x where x = {fa, omega, stats, vss}
#reorganized, January 18, 2009 to make some what clearer
#added the switch capability, August 25, 2011 following suggestions by Joshua Wiley

"print.psych" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,short=TRUE,lower=TRUE,...) { 

if(length(class(x)) > 1)  { value <- class(x)[2] } else {
#these next test for non-psych functions that may be printed using print.psych.fa
if((!is.null(x$communality.iterations)) | (!is.null(x$uniquenesses)) | (!is.null(x$rotmat)) |  (!is.null(x$Th)) ) {value <- fa }
}
if(all) value <- "all"
if(value == "score.items") value <- "scores"
if(value =="set.cor") value <- "setCor"
switch(value,

## the following functions have their own print function
 extension = { print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)},
 extend = {print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)},
 fa =  {print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)},
 fa.ci = { print.psych.fa.ci(x,digits=digits,all=all,... )},
 iclust= { print.psych.iclust(x,digits=digits,all=all,cut=cut,sort=sort,...)},
 omega = { print.psych.omega(x,digits=digits,all=all,cut=cut,sort=sort,...)},
 omegaSem=  {print.psych.omegaSem(x,digits=digits,all=all,cut=cut,sort=sort,...)},
  principal ={print.psych.fa(x,digits=digits,all=all,cut=cut,sort=sort,...)},
  schmid = { print.psych.schmid(x,digits=digits,all=all,cut=cut,sort=sort,...)},
  stats = { print.psych.stats(x,digits=digits,all=all,cut=cut,sort=sort,...)},
  vss= { print.psych.vss(x,digits=digits,all=all,cut=cut,sort=sort,...)},
  cta = {print.psych.cta(x,digits=digits,all=all,...)},


##Now, for the smaller print jobs, just do it here.
all=  {class(x) <- "list"
         print(x) },    #find out which function created the data and then print accordingly
         
alpha = {
	cat("\nReliability analysis ",x$title," \n")
	cat("Call: ")
	print(x$call)
	cat("\n ")
	print(x$total,digits=digits)
	if(!is.null(x$total$ase)){ cat("\n lower alpha upper     95% confidence boundaries\n")
	cat(round(c(x$total$raw_alpha - 1.96* x$total$ase, x$total$raw_alpha,x$total$raw_alpha +1.96* x$total$ase),digits=digits) ,"\n")}
	if(!is.null(x$boot.ci)) {cat("\n lower median upper bootstrapped confidence intervals\n",round(x$boot.ci,digits=digits))}
	cat("\n Reliability if an item is dropped:\n")
    print(x$alpha.drop,digits=digits)
	cat("\n Item statistics \n")
   print(x$item.stats,digits=digits)
	 if(!is.null(x$response.freq)) {
	 cat("\nNon missing response frequency for each item\n")
	 print(round(x$response.freq,digits=digits))}
},
bestScales ={ 
     df <- data.frame(correlation=x$r,n.items = x$n.items)
    cat("The items most correlated with the criteria yield r's of \n")
    print(round(df,digits=digits)) 
    if(length(x$value)> 0) {cat("\nThe best items, their correlations and content  are \n")
    print(x$value) } else {cat("\nThe best items and their correlations are \n")
     for(i in 1:length(x$short.key)) {print(round(x$short.key[[i]],digits=digits))} 
     } 
      },
      

bifactor = {
   cat("Call: ")
     print(x$Call) 
 	cat("Alpha:                ",round(x$alpha,digits),"\n") 
 	cat("G.6:                  ",round(x$G6,digits),"\n")
 	cat("Omega Hierarchical:   " ,round(x$omega_h,digits),"\n")
 #	cat("Omega H asymptotic:   " ,round(x$omega.lim,digits),"\n")
 	cat("Omega Total           " ,round(x$omega.tot,digits),"\n")
 	print(x$f,digits=digits,sort=sort)
 	},
         
circ =     {cat("Tests of circumplex structure \n")
            cat("Call:")
             print(x$Call)
             res <- data.frame(x[1:4]) 
             print(res,digits=2)
             },
             
circadian =  {if(!is.null(x$Call)) {cat("Call: ")
                 print(x$Call)}
      	cat("\nCircadian Statistics :\n")
      	
    if(!is.null(x$F)) {
        cat("\nCircadian F test comparing groups :\n")
         print(round(x$F,digits)) 
         if(short)  cat("\n To see the pooled and group statistics, print with the short=FALSE option")
         }	
         
	if(!is.null(x$pooled) && !short) { cat("\nThe pooled circadian statistics :\n")
	print(  x$pooled)}
	
	if(!is.null(x$bygroup) && !short) {cat("\nThe  circadian statistics by group:\n")
	print(x$bygroup)}
	#if(!is.null(x$result)) print(round(x$result,digits))
	if(!is.null(x$phase.rel)) {
	     cat("\nSplit half reliabilities are split half correlations adjusted for test length\n")
	     x.df <- data.frame(phase=x$phase.rel,fits=x$fit.rel)
	     print(round(x.df,digits)) }
	if(is.data.frame(x)) {class(x) <- "data.frame"
            print(round(x,digits=digits)) }	
	},
	
cluster.cor =  {
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
	cat("\nSignal to Noise ratio based upon average r and n \n")
	print(x$sn,digits=digits)
	
	      
 	 # cat("\nScale intercorrelations:\n")
	 # print(x$cor,digits=digits)
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	  },
	

cluster.loadings =  {
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
  },


comorbid = {cat("Call: ")
              print(x$Call)
            cat("Comorbidity table \n")
            print(x$twobytwo,digits=digits)
            cat("\nimplies phi = ",round(x$phi,digits), " with Yule = ", round(x$Yule,digits), " and tetrachoric correlation of ", round(x$tetra$rho,digits))
            cat("\nand normal thresholds of ",round(-x$tetra$tau,digits))
          
   },



     
cor.ci = {cat("Call:")
              print(x$Call) 
               cat("\n Coefficients and bootstrapped confidence intervals \n")
               lowerMat(x$rho)
   phis <- x$rho[lower.tri(x$rho)]
   cci <- data.frame(lower.emp =x$ci$low.e, lower.norm=x$ci$lower,estimate =phis ,upper.norm= x$ci$upper, upper.emp=x$ci$up.e,p = x$ci$p)
   rownames(cci) <- rownames(x$ci)
   cat("\n scale correlations and bootstrapped confidence intervals \n")
   
   print(round(cci,digits=digits))
              
            },  
            
cor.cip = {class(x) <- NULL
   cat("\n High and low confidence intervals \n")
   print(round(x,digits=digits))
              
            },       
     
corr.test = {cat("Call:")
              print(x$Call)
              cat("Correlation matrix \n")
              print(round(x$r,digits))
              cat("Sample Size \n")
              print(x$n)
              if(x$sym) {cat("Probability values (Entries above the diagonal are adjusted for multiple tests.) \n")} else {
                 if (x$adjust != "none" )  {cat("Probability values  adjusted for multiple tests. \n")}}
             print(round(x$p,digits))
             cat("\n To see confidence intervals of the correlations, print with the short=FALSE option\n")
             if(!short) {cat("\n Confidence intervals based upon normal theory.  To get bootstrapped values, try cor.ci\n")
             print(round(x$ci,digits)) }
         },     

corr.p = {cat("Call:")
              print(x$Call)
              cat("Correlation matrix \n")
              print(round(x$r,digits))
              cat("Sample Size \n")
              print(x$n)
              if(x$sym) {cat("Probability values (Entries above the diagonal are adjusted for multiple tests.) \n")} else {
                 if (x$adjust != "none" )  {cat("Probability values  adjusted for multiple tests. \n")}}
             print(round(x$p,digits))
             cat("\n To see confidence intervals of the correlations, print with the short=FALSE option\n")
              if(!short) {cat("\n Confidence intervals based upon normal theory.  To get bootstrapped values, try cor.ci\n")
             print(round(x$ci,digits)) }

         },     
         
cortest= {cat("Tests of correlation matrices \n")
            cat("Call:")
            print(x$Call)
            cat(" Chi Square value" ,round(x$chi,digits)," with df = ",x$df, "  with probability <", signif(x$p,digits),"\n" )
            if(!is.null(x$z)) cat("z of differences = ",round(x$z,digits),"\n")
         },



         

describe= {if  (length(dim(x))==1) {class(x) <- "list"
              attr(x,"call") <- NULL
              print(round(x,digits=digits))
                  } else  {class(x) <- "data.frame" 
            print(round(x,digits=digits)) }
         },
         
describeData = {if  (length(dim(x))==1) {class(x) <- "list"
              attr(x,"call") <- NULL
              print(round(x,digits=digits))
                  } else  {
                  cat('n.obs = ', x$n.obs, "of which ", x$complete.cases,"  are complete cases.   Number of variables = ",x$nvar," of which all are numeric ",x$all.numeric," \n")
                 
            print(x$variables) }
         },
 
faBy = { cat("Call: ")
        print(x$Call)
        cat("\n Factor analysis by Groups\n")
        cat("\nAverage standardized loadings (pattern matrix) based upon correlation matrix for all cases as well as each group\n")
        cat("\nlow and high ", x$quant,"% quantiles\n")
        print(x$faby.sum,digits)
        if(!short) {  
        print(x$mean.loading,digits=digits)
        cat("\n Average factor intercorrelations for all cases and  each group\n")
        print(x$mean.Phi,digits=2)
        cat("\nStandardized loadings (pattern matrix) based upon correlation matrix for all cases as well as each group\n")
        print(x$loadings,digits=digits)
        cat("\n With factor intercorrelations for all cases and for each group\n")
        print(x$Phi,digits=2)
        if(!is.null(x$fa)) {
        cat("\nFactor analysis results for each group\n")
        print(x$faby.sum,digits)
        }}
        
       
# cat("Factor Analysis with confidence intervals using method = ",x$f$fm )
#   # cat("\nCall: ")
#
#   
#   nfactors <-dim(x$loadings)[2]
#    c("\n Confidence intervals\n")
#    if(is.null(x[["ci"]])) { lc <- lci <- data.frame(unclass(x$loadings),x$low,x$highr)} else { lc <- lci <- data.frame(unclass(x$loadings),x$ci)}
#   
#   
#    for(i in 1:nfactors) {
#    lci[,(i-1)*3 +2 ] <- lc[,i] 
#    lci[,(i-1)*3 +1 ] <- lc[,i+nfactors] 
#    lci[,(i-1)*3 +3 ] <- lc[,i+nfactors*2] 
#     }
#     colnames(lci) <- paste(rep(c("low","coeff","upper"),nfactors),sep="")
#    for(i in 1:nfactors) { colnames(lci)[(i-1)*3+2] <- colnames(x$loadings)[i] }
#    
#    cat("\n Coefficients and bootstrapped confidence intervals \n")
#    print (round(lci,digits=digits))
#   if(!is.null(x$Phi)) { phis <- x$Phi[lower.tri(x$Phi)]
#    cci <- data.frame(lower=x$cis$ci.rot$lower,estimate = phis,upper= x$cis$ci.rot$upper)
#     cnR <- abbreviate(colnames(x$loadings),minlength=5) 
#       k <- 1
#      for(i in 1:(nfactors-1)) {for (j in (i+1):nfactors) {
#       rownames(cci)[k] <- paste(cnR[i],cnR[j],sep="-")
#       k<- k +1 }}  #added 10/4/14
#    cat("\n Interfactor correlations and bootstrapped confidence intervals \n")
#    print(cci,digits=digits)}
#   
#}
        
},
    
guttman =  {
  cat("Call: ")
    print(x$Call)
 cat("\nAlternative estimates of reliability\n")
# cat("Beta = ", round(x$beta,digits), " This is an estimate of the worst split half reliability")  
 cat("\nGuttman bounds \nL1 = ",round(x$lambda.1,digits), "\nL2 = ", round(x$lambda.2,digits), "\nL3 (alpha) = ", round(x$lambda.3,digits),"\nL4 (max) = " ,round(x$lambda.4,digits), "\nL5 = ", round(x$lambda.5,digits), "\nL6 (smc) = " ,round(x$lambda.6,digits), "\n")
 cat("TenBerge bounds \nmu0 = ",round(x$tenberge$mu0,digits), "mu1 = ", round(x$tenberge$mu1,digits), "mu2 = " ,round(x$tenberge$mu2,digits), "mu3 = ",round(x$tenberge$mu3,digits) , "\n")
 cat("\nalpha of first PC = ",round( x$alpha.pc,digits), "\nestimated greatest lower bound based upon communalities= ", round(x$glb,digits),"\n")
 cat("\nbeta found by splitHalf  = ", round(x$beta,digits),"\n")
 } ,
 
 ICC =  {cat("Call: ")
              print(x$Call)
            cat("\nIntraclass correlation coefficients \n")
            print(x$results,digits=digits)
            cat("\n Number of subjects =", x$n.obs, "    Number of Judges = ",x$n.judge)

   },
   
iclust.sort = {
     nvar <- ncol(x$sort) 
     x$sort[4:nvar] <- round(x$sort[4:nvar],digits)
     print(x$sort)
       },

   
irt.fa = {
   cat("Item Response Analysis using Factor Analysis \n")
   cat("\nCall: ")
   print(x$Call)
   if (!is.null(x$plot)) print(x$plot)
    if(!short) {
    nf <- length(x$irt$difficulty)
    for(i in 1:nf) {temp <- data.frame(discrimination=x$irt$discrimination[,i],location=x$irt$difficulty[[i]])
    cat("\nItem discrimination and location for factor ",colnames(x$irt$discrimination)[i],"\n")
    print(round(temp,digits))}

    cat("\n These parameters were based on the following factor analysis\n")
    print(x$fa)
    } else {summary(x$fa)}
  },
  
irt.poly =  {
   cat("Item Response Analysis using Factor Analysis  \n")
   cat("\nCall: ")
   print(x$Call)
    if (!is.null(x$plot)) print(x$plot) #this calls the polyinfo print function below
     if(!short) {
    nf <- length(x$irt$difficulty)
    for(i in 1:nf) {temp <- data.frame(discrimination=x$irt$discrimination[,i],location=x$irt$difficulty[[i]])
    cat("\nItem discrimination and location for factor ",colnames(x$irt$discrimination)[i],"\n")
    print(round(temp,digits))}
   
    cat("\n These parameters were based on the following factor analysis\n")
    print(x$fa)
    } else {summary(x$fa) }
  },
 
kappa = {if(is.null(x$cohen.kappa)) {
            cat("Call: ")
            print(x$Call)
            
            cat("\nCohen Kappa and Weighted Kappa correlation coefficients and confidence boundaries \n")
           
            print(x$confid,digits=digits)
            cat("\n Number of subjects =", x$n.obs,"\n")} else {
            cat("\nCohen Kappa (below the diagonal) and Weighted Kappa (above the diagonal) \nFor confidence intervals and detail print with all=TRUE\n")
            print(x$cohen.kappa,digits) 
            }
   },
   
mardia =  {
   cat("Call: ")
     print(x$Call) 
     cat("\nMardia tests of multivariate skew and kurtosis\n")
     cat("Use describe(x) the to get univariate tests")
      cat("\nn.obs =",x$n.obs,"  num.vars = ",x$n.var,"\n")
     cat("b1p = ",round(x$b1p,digits),"  skew = ",round(x$skew,digits ), " with probability = ", signif(x$p.skew,digits)) 
     cat("\n small sample skew = ",round(x$small.skew,digits ), " with probability = ", signif(x$p.small,digits)) 
     cat("\nb2p = ", round(x$b2p,digits),"  kurtosis = ",round(x$kurtosis,digits)," with probability = ",signif(x$p.kurt,digits ))
  },

  
mchoice =  {
    cat("Call: ")
    print(x$Call)
	cat("\n(Unstandardized) Alpha:\n")
	print(x$alpha,digits=digits)
  	cat("\nAverage item correlation:\n")
  	print(x$av.r,digits=digits)
	 if(!is.null(x$item.stats)) {
	 cat("\nitem statistics \n")
	 print(round(x$item.stats,digits=digits))}
  },  
  
mediate  = { cat("\nMediation analysis \n")

cat("Call: ")
    print(x$Call)
    dv <- x$names[1]
    iv <- rownames(x$direct)
    mv <- x$names[-c(1:(length(iv)+1))]
    cat("\nThe DV (Y) was ", dv,". The IV (X) was ", iv,". The mediating variable(s) = ", mv,".")
    if(!is.null(x$mod)) cat("  The moderating variable(s) = ", x$names[x$mod])
    
   for(i in 1:length(iv)) { cat("\n\nTotal Direct effect(c) of ",iv[i], " on ", dv," = ",round(x$total[i],digits), "  S.E. = ", round(x$se.bx[i],digits), " t direct = ",round(x$tt[i],digits), "  with probability = ", signif(x$probt[i],digits))
    cat("\nDirect effect (c') of ",iv[i],  " on ", dv," removing ", mv ," = ",round(x$direct[i],digits), "  S.E. = ", round(x$se.beta[i],digits), " t direct = ",round(x$t[i],digits), "  with probability = ", signif(x$prob[i],digits))
     
   if(is.null(x$mod)) { cat("\nIndirect effect (ab) of ",iv[i], " on ", dv," through " ,mv , "  = ", round(x$indirect[i],digits),"\n")
   cat("Mean bootstrapped indirect effect = ",round(x$mean.boot[i],digits), " with standard error = ",round(x$sd.boot[i],digits), " Lower CI = ",round(x$ci.quant[1,i],digits), "   Upper CI = ", round(x$ci.quant[2,i],digits))}
     }
   if(is.null(x$mod)) {
     cat("\nSummary of a, b, and ab estimates and ab confidence intervals\n")
     if(!any(is.na(x$abc))) {print(round(x$abc,digits))} else {
      cat("\n 'a'  paths \n")
      print(round(x$a,digits))
      cat("\n'b' paths \n")
      print(round(x$b,digits))
      cat("\n'ab' paths \n")
      print(round(x$ab,digits))
     }
    
    cat("\nratio of indirect to total effect=  ", round(x$ratit,digits))
    cat("\nratio of indirect to direct effect= ", round(x$ratid,digits))
    }  else {
    cat("\nEffect of interaction of ",iv[1], " with ", iv[2] , "  = ", round(x$direct[3],digits),"  S.E. = ", round(x$se.beta[3,1],digits), " t direct = ",round(x$t[3,1],digits), "  with probability = ", signif(x$prob[3,1],digits))
    cat("\nIndirect effect due to interaction  of ",iv[1], " with ", iv[2] , "  = ", round(x$int.ind,digits))
    cat("\nMean bootstrapped indirect interaction effect = ",round(x$mean.boot[1],digits), " with standard error = ",round(x$sd.boot[1],digits), " Lower CI = ",round(x$ci.quant[1],digits), "   Upper CI = ", round(x$ci.quant[2,i],digits))
    cat("\nSummary of a, b, and ab estimates and ab confidence intervals\n")
   if(!is.na(x$abc)) {print(round(x$abc,digits))} else {
      print(round(x$a,digits))
      print(round(x$b,digits))
      print(round(x$ab,digits))
     }
   
    cat("\nR2 of model = ", round(x$R2,digits)) 
   }},
      
    
mixed= { cat("Call: ")
    print(x$Call)
    if(is.null(x$rho)) {if(lower) {print(lowerMat(x,digits=digits))} else {print(x,digits)} } else {
    if(lower) {if(length(x$rho)>1) print(lowerMat (x$rho),digits=digits)} else {print(x$rho,digits)}}
   },
   
 
 paired.r = {cat("Call: ")
            print(x$Call)
        print(x$test)
        if(is.null(x$z)) {cat("t =",round(x$t,digits))  
         } else {cat("z =",round(x$z,digits)) }
         cat("  With probability = ",round(x$p,digits))
         },

  

         
parallel= {
cat("Call: ")
              print(x$Call) 
              if(!is.null(x$fa.values) & !is.null(x$pc.values) ) {
                  parallel.df <- data.frame(fa=x$fa.values,fa.sim =x$fa.sim,pc= x$pc.values,pc.sim =x$pc.sim)
		fa.test <- x$nfact
		pc.test <- x$ncomp
		
		cat("Parallel analysis suggests that ")
		cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
                  cat("\n Eigen Values of \n")
                  
                  colnames(parallel.df) <- c("Original factors","Simulated data","Original components", "simulated data")}
              if(is.na(fa.test) ) fa.test <- 0
              if(is.na(pc.test)) pc.test <- 0
               if(!any(is.na(parallel.df)))  {print(round(parallel.df[1:max(fa.test,pc.test),],digits))}  else {
              if(!is.null(x$fa.values)) {cat("\n eigen values of factors\n")
              print(round(x$fa.values,digits))}
               if(!is.null(x$fa.sim)){cat("\n eigen values of simulated factors\n") 
               print(round(x$fa.sim,digits))}
              if(!is.null(x$pc.values)){cat("\n eigen values of components \n")
              print(round(x$pc.values,digits))}
             
              if(!is.null(x$pc.sim)) {cat("\n eigen values of simulated components\n") 
              print(round(x$pc.sim,digits=digits))}
            }  
    },
         
  partial.r = {cat("partial correlations \n")
            print(round(unclass(x),digits))
         }, 
         
 phi.demo = {print(x$tetrachoric)
            cat("\nPearson (phi) below the diagonal, phi2tetras above the diagonal\n")
             print(round(x$phis,digits))
             cat("\nYule correlations")
             print(x$Yule)
},
         
 poly= {cat("Call: ")
              print(x$Call)
            cat("Polychoric correlations \n")
            if(!is.null(x$twobytwo)) {
              print(x$twobytwo,digits=digits)
              cat("\n implies tetrachoric correlation of ",round(-x$rho,digits))} else {
            
            if(lower) {lowerMat (x$rho,digits) } else {print(x$rho,digits)}
            cat("\n with tau of \n")
            print(x$tau,digits)
            }
   },

polydi= {cat("Call: ")
              print(x$Call)
            cat("Correlations of polytomous with dichotomous\n")
            print(x$rho,digits)
            cat("\n with tau of \n")
            print(x$tau,digits)
   },
   
 polyinfo =  {cat("Item Response Analysis using Factor Analysis  \n")
   cat("\n Summary information by factor and item")
     names(x$sumInfo ) <- paste("Factor",1:length(x$sumInfo))
    for (f in 1:length(x$sumInfo)) {
     cat("\n Factor = ",f,"\n")
     temp <- x$sumInfo[[f]]
     temps <- rowSums(temp)
     if(sort) {ord <- order(temps,decreasing=TRUE)
              temp <- temp[ord,]
              temps <- temps[ord]}
     temp <- temp[temps > 0,]
    
     summary <- matrix(c(colSums(temp),sqrt(1/colSums(temp)),1-1/colSums(temp)),nrow=3,byrow=TRUE)
     rownames(summary) <-c("Test Info","SEM", "Reliability")
     temp <- rbind(temp,summary)
     if(ncol(temp) == 61) {print(round(temp[,seq(1,61,10)],digits=digits)) } else {print(round(temp,digits=digits))} #this gives us info at each unit 
     }
 if(!short) {    
     cat("\n Average information (area under the curve) \n")
     AUC <-x$AUC
     max.info <-x$max.info
     if(dim(AUC)[2]==1) {item <- 1:length(AUC) } else {item <- 1:dim(AUC)[1]}
     if(sort) {
 		#first sort them into clusters
  		#first find the maximum for each row and assign it to that cluster
  		 cluster <- apply(AUC,1,which.max)
 		 ord <- sort(cluster,index.return=TRUE)
  		 AUC <-  AUC[ord$ix,,drop=FALSE]
  		 max.info <-  max.info[ord$ix,,drop=FALSE]
      #now sort column wise
      #now sort the AUC that have their highest AUC on each cluster

  		items <- table(cluster)   #how many items are in each cluster?
  		first <- 1
    	for (i in 1:length(items)) {# i is the factor number
		if(items[i] > 0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(AUC[first:last,i]),decreasing=TRUE,index.return=TRUE)
				
   				AUC[first:last,] <- AUC[item[ord$ix+first-1],]
   				max.info[first:last,] <- max.info[item[ord$ix+first-1],]
   				rownames(AUC)[first:last] <- rownames(max.info)[first:last]  <- rownames(AUC)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
         }    #end of sort 		 
     print(AUC,digits=digits)
     cat("\nMaximum value is at \n")
     print(max.info,digits=digits)
    }
    },

overlap =  {
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
	cat("\nSignal to Noise ratio based upon average r and n \n")
	print(x$sn,digits=digits)
	
	 cat("\nScale intercorrelations corrected for item overlap and attenuation \n adjusted for overlap correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	
	 if(short) {cat("\n In order to see the item by scale loadings and frequency counts of the data\n print with the short option = FALSE") } else {
	  if(!is.null(x$item.cor) ) {
	   cat("\nItem by scale correlations:\n corrected for item overlap and scale reliability\n" )
	   print(round(x$item.cor,digits=digits)) } 
	 }
	
	  },
	   
                  
r.test =  {cat("Correlation tests \n")
            cat("Call:")
              print(x$Call)
              cat( x$Test,"\n")
              if(!is.null(x$t)) {cat(" t value" ,round(x$t,digits),"   with probability <", signif(x$p,digits) )}
              if(!is.null(x$z)) {cat(" z value" ,round(x$z,digits),"   with probability ",  round(x$p,digits) )}               
              if(!is.null(x$ci)) {cat("\n and confidence interval ",round(x$ci,digits) ) }
         },   
         
residuals = { 
   if (lower) {lowerMat (x,digits=digits)} else {print(x,digits)}
	},   
	
	
scree = {
       cat("Scree of eigen values \nCall: ") 
       print(x$Call) 
       if(!is.null(x$fv)) {cat("Eigen values of factors ")
              print(round(x$fv,digits))}
       if (!is.null(x$pcv)) {cat("Eigen values of Principal Components")
              print(round(x$pcv,digits))}
         },  
scores =  {
    cat("Call: ")
    print(x$Call)
    if(x$raw) {
	cat("\n(Unstandardized) Alpha:\n") } else {cat("\n(Standardized) Alpha:\n") }
	print(x$alpha,digits=digits) 
	if(!is.null(x$ase)) {cat("\nStandard errors of unstandardized Alpha:\n")
	rownames(x$ase) <- "ASE  "
	print(x$ase,digit=digits) }
	if(!is.null(x$alpha.ob)) {cat("\nStandardized Alpha of observed scales:\n")
	print(x$alpha.ob,digits=digits)}
  	cat("\nAverage item correlation:\n")
  	print(x$av.r,digits=digits)
	cat("\n Guttman 6* reliability: \n")
	print(x$G6,digits=digits)   
	cat("\nSignal/Noise based upon av.r : \n")
	print(x$sn,digits=digits)  
	#if(iclust) {cat("\nOriginal Beta:\n")
	# print(x$beta,digits) }	          
 	
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 	if(!is.null(x$alpha.ob)) {cat("\nNote that these are the correlations of the complete scales based on the correlation matrix,\n not the observed scales based on the raw items.\n")}
	 	
	 print(x$corrected,digits) 
	 if(short) {cat("\n In order to see the item by scale loadings and frequency counts of the data\n print with the short option = FALSE") } else {
	  if(!is.null(x$item.cor) ) {
	   cat("\nItem by scale correlations:\n corrected for item overlap and scale reliability\n" )
	   
	 print(round(x$item.corrected,digits=digits)) } 
	 if(!is.null(x$response.freq)) {
	 cat("\nNon missing response frequency for each item\n")
	 print(round(x$response.freq,digits=digits))}
	 }
	
  },
  
setCor= { cat("Call: ")
              print(x$Call)
            if(x$raw) {cat("\nMultiple Regression from raw data \n")} else {
            cat("\nMultiple Regression from matrix input \n")}
           cat("\nBeta weights \n")
           print(round(x$beta,digits))
           cat("\nMultiple R \n") 
           print(round(x$R,digits))
            cat("multiple R2 \n") 
            print(x$R2,digits)
             cat("\n Unweighted multiple R \n") 
           print(round(x$ruw,digits))
             cat(" Unweighted multiple R2 \n") 
           print(round(x$ruw^2,digits))
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
            
            if(!is.null(x$cancor)) {
            cat("\nVarious estimates of between set correlations\n")
            cat("Squared Canonical Correlations \n")
            print(x$cancor2,digits=digits)
            if(!is.null(x$Chisq)) {cat("Chisq of canonical correlations \n")
            print(x$Chisq,digits=digits)}
            cat("\n Average squared canonical correlation = ",round(x$T,digits=digits))
          
            cat("\n Cohen's Set Correlation R2 = ",round(x$Rset,digits=digits))
            #print(x$Rset,digits=digits)
           if(!is.null(x$Rset.shrunk)){ cat("\n Shrunken Set Correlation R2 = ",round(x$Rset.shrunk,digits=digits))
          
            cat("\n F and df of Cohen's Set Correlation ",round(c(x$Rset.F,x$Rsetu,x$Rsetv), digits=digits))}
             cat("\nUnweighted correlation between the two sets = ",round(x$Ruw,digits)) 
           
           }


   },
   
   
sim =  { if(is.matrix(x)) {x <-unclass(x) 
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
 },
 
 smoother = {x <- unclass(x)
             print(x)
             },
             
 
 split ={ cat("Split half reliabilities  ")
   cat("\nCall: ")
      print(x$Call)
   cat("\nMaximum split half reliability (lambda 4) = ",round(x$maxrb,digits=digits))
   cat("\nGuttman lambda 6                          = ",round(x$lambda6,digits=digits))
   cat("\nAverage split half reliability            = ",round(x$meanr,digits=digits))
   cat("\nGuttman lambda 3 (alpha)                  = ",round(x$alpha,digits=digits))
   cat("\nMinimum split half reliability  (beta)    = ",round(x$minrb,digits=digits))
      
 },
 
statsBy ={
 cat("Statistics within and between groups  ")
   cat("\nCall: ")
   print(x$Call)
 cat("Intraclass Correlation 1 (Percentage of variance due to groups) \n")
 print(round(x$ICC1,digits))
cat("Intraclass Correlation 2 (Reliability of group differences) \n")
print(round(x$ICC2,digits))
cat("eta^2 between groups  \n")
print(round(x$etabg^2,digits))
if(!short) {cat("Correlation between groups \n")
lowerMat(x$rbg)
cat("Correlation within groups \n")
lowerMat(x$rwg)
}
},

            
   
 tetra = {cat("Call: ")
              print(x$Call)
            cat("tetrachoric correlation \n")
            if(!is.null(x$twobytwo)) {
              print(x$twobytwo,digits=digits)
              cat("\n implies tetrachoric correlation of ",round(x$rho,digits))} else {
           if(is.matrix(x$rho) &&  lower) {lowerMat (x$rho,digits)} else { print(x$rho,digits)}
            cat("\n with tau of \n")
            print(x$tau,digits)
          }
   },
   
thurstone = {
	cat("Thurstonian scale (case 5) scale values ")
	cat("\nCall: ")
     print(x$Call)
	print(x$scale)
	cat("\n Goodness of fit of model  ", round(x$GF,digits))
 },
         

   


       
KMO = {cat("Kaiser-Meyer-Olkin factor adequacy")
     cat("\nCall: ")
   print(x$Call)
   cat("Overall MSA = ",round(x$MSA,digits))
   cat("\nMSA for each item = \n")
   print(round(x$MSAi,digits))
   },
   
yule = {cat("Yule and Generalized Yule coefficients")
     cat("\nCall: ")
   print(x$Call)
   cat("\nYule coefficient \n")
   print(round(x$rho,digits))
   cat("\nUpper and Lower Confidence Intervals = \n")
   print(round(x$ci,digits))
   },
   
Yule = {cat("Yule and Generalized Yule coefficients")
   
   cat("\nLower CI  Yule coefficient Upper CI \n")
   print(round(c(x$lower,x$rho,x$upper),digits))
}
  )   #end of switch 
}  #end function
