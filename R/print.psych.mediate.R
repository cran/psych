"print.psych.mediate" <- function(x,digits=2,short=TRUE) {
 cat("\nMediation/Moderation Analysis \nCall: ")
    print(x$Call)
    dv <- x$var.names[["DV"]]
   # iv <- x$var.names[["IV"]]
    mv <- x$var.names[["med"]]
    mod <- x$var.names[["mod"]]
   # dv <- x$names[1]
    iv <- rownames(x$direct)
    if(iv[1] == "Intercept") iv <- iv[-1]
    niv <- length(iv)
    nmed <- length(mv)
    ndv <- length(dv)
    nz <- length(x$var.names[["z"]])
   
   # if(dim(x$a)) {mv <- names(x$a)} else {mv <- colnames(x$a)
    cat("\nThe DV (Y) was ", dv,". The IV (X) was ", iv,". The mediating variable(s) = ", mv,".")
   if(!is.null(x$mod)) cat("  The moderating variable(s) = ",mod)
   if(!is.null(x$var.names$z))  cat(" Variable(s)  partialled out were", x$var.names[["z"]])
  
  
  
  if(!is.null(mv)) {
   for(j in 1:ndv) { 
   for(i in 1:niv) { cat("\n\nTotal effect(c) of ",iv[i], " on ", dv[j]," = ",round(x$direct[i+1,j],digits), "  S.E. = ", round(x$total.reg$se[i+1,j],digits), " t  = ",round(x$total.reg$t[i+1,j],digits)," df= ",x$total.reg$df, "  with p = ", signif(x$total.reg$prob[i+1,j],digits))
    cat("\nDirect effect (c') of ",iv[i],  " on ", dv[j]," removing ", mv ," = ",round(x$cprime.reg$beta[i+1,j],digits), "  S.E. = ", round(x$cprime.reg$se[i+1,j],digits), " t  = ",round(x$cprime.reg$t[i+1,j],digits), " df= ", x$cprime.reg$df, "  with p = ", signif(x$cprime.reg$prob[i+1,j],digits))
     
   if(is.null(x$mod)) { cat("\nIndirect effect (ab) of ",iv[i], " on ", dv[j]," through " ,mv , "  = ", round(x$ab[i,j],digits),"\n")
   cat("Mean bootstrapped indirect effect = ",round(x$boot$mean[i],digits), " with standard error = ",round(x$boot$sd[i],digits), " Lower CI = ",round(x$boot$ci[1,i],digits), "   Upper CI = ", round(x$boot$ci[2,i],digits))}
     }
     
     F <-  x$cprime.reg$df * x$cprime.reg$R2[j]/(((nrow(x$cprime.reg$beta)-1) * (1-x$cprime.reg$R2[j])))
      pF <-  -expm1(pf(F,nrow(x$cprime.reg$beta),x$cprime.reg$df,log.p=TRUE)) 
      cat("\nR =", round(sqrt(x$cprime.reg$R2[j]),digits),"R2 =", round(x$cprime.reg$R2[j],digits),  "  F =", round(F,digits), "on",nrow(x$cprime.reg$beta)-1, "and", x$cprime.reg$df,"DF   p-value: ",signif(pF,digits+1), "\n") 
    
     
     }
     if(short) {cat("\n To see the longer output, specify short = FALSE in the print statement or ask for the summary")} else {
   
 if(is.null(x$mod)) {

 
    
    cat("\n\n Full output  \n")
    
      summary.psych.mediate(x)
    
    if(FALSE) {
    cat("\nDirect effect estimates (traditional regression)    (c') \n")
     for(j in 1:ndv) {
     
    if (niv==1) { dfd <- round(data.frame(direct=x$cprime.reg$beta[,j],se = x$cprime.reg$se[,j],t=x$cprime.reg$t[,j],df=x$cprime.reg$df),digits)
     dfdp <- cbind(dfd,p=signif(x$cprime.reg$prob[,j],digits=digits+1)) } else {
     dfd <- round(data.frame(direct=x$cprime.reg$beta[1:(niv+1+nmed),j],se = x$cprime.reg$se[1:(niv+1+nmed),j],t=x$cprime.reg$t[1:(niv+1+nmed),j],df=x$cprime.reg$df),digits)
     dfdp <- cbind(dfd,p=signif(x$cprime.reg$prob[1:(niv+1+nmed),j],digits=digits+1))
     }
      colnames(dfdp) <- c(dv[j],"se","t","df","Prob")
     
   print(dfdp)
     F <-  x$cprime.reg$df * x$cprime.reg$R2[j]/(((nrow(x$cprime.reg$beta)-1) * (1-x$cprime.reg$R2[j])))
      pF <-  -expm1(pf(F,nrow(x$cprime.reg$beta)-1,x$cprime.reg$df,log.p=TRUE)) 
      cat("\nR =", round(sqrt(x$cprime.reg$R2[j]),digits),"R2 =", round(x$cprime.reg$R2[j],digits),  "  F =", round(F,digits), "on",nrow(x$cprime.reg$beta)-1, "and", x$cprime.reg$df,"DF   p-value: ",signif(pF,digits+1), "\n") 
    
   
   }
 
     cat("\n Total effect estimates (c) \n")
      
        for(j in 1:ndv) {

    dft <- round(data.frame(direct=x$total.reg$beta[,j],se = x$total.reg$se[,j],t=x$total.reg$t[,j],df=x$total.reg$df),digits)
    dftp <- cbind(dft,p = signif(x$total.reg$prob[,j],digits=digits+1))
    colnames(dftp) <- c(dv[j],"se","t","df","Prob")
    rownames(dftp) <- rownames(x$total.reg$beta)
     print(dftp)
    }

    
     
    cat("\n 'a'  effect estimates \n")

  if(niv==1) {
    	dfa <- round(data.frame(a = x$a.reg$beta[1,1:nmed],se = x$a.reg$se[1,1:nmed],t = x$a.reg$t[1,1:nmed],df= x$a.reg$df),digits)
    	dfa <- cbind(dfa,p=signif(x$a.reg$prob[1,1:nmed],digits=digits+1))
    	if(NROW(dfa) ==1) {rownames(dfa) <- rownames(x$a.reg$beta)[-1]
    	colnames(dfa) <-  c(colnames(x$a.reg$beta),"se","t","df", "Prob")} else {
    	rownames(dfa) <- colnames(x$a.reg$beta)
    	colnames(dfa) <-  c(rownames(x$a.reg$beta),"se","t","df", "Prob")}
    	
    	print(dfa)}  else {
    	
     	for (i in 1:nmed) {
     	dfa <- round(data.frame(a = x$a.reg$beta[,i],se = x$a.reg$se[,i],t = x$a.reg$t[,i],df= x$a.reg$df),digits)
    	dfa <- cbind(dfa,p=signif(x$a.reg$prob[,i],digits=digits+1))
     	rownames(dfa) <-rownames(x$a.reg$beta)
     	colnames(dfa) <-  c(colnames(x$a.reg$beta)[i],"se","t","df","Prob") 
     	print(dfa) }
     	
     	}
     	        
      cat("\n 'b'  effect estimates \n")
      for (j in 1:ndv) {
      if(niv==1) {
     dfb <- round(data.frame(direct=x$b.reg$beta[-(1:niv),j],se = x$b.reg$se[-(1:niv),j],t=x$b.reg$t[-(1:niv),j], df=x$b.reg$df),digits)
     dfb <- cbind(dfb,p=signif(x$b.reg$prob[-(1:niv),j],digits=digits+1))} else {
      dfb <- round(data.frame(direct=x$b.reg$beta[-(1:niv),j],se = x$b.reg$se[-(1:niv),j],t=x$b.reg$t[-(1:niv),j],df=x$b.reg$df),digits)
     dfb <- cbind(dfb,p=signif(x$b.reg$prob[-(1:niv),j],digits=digits+1))}
     rownames(dfb) <- rownames(x$b.reg$beta)[-(1:niv)]
     colnames(dfb) <-  c(dv[j],"se","t","df", "Prob")
      print(dfb)
      }
 
      cat("\n 'ab'  effect estimates \n")
     
 for (j in 1:ndv) {
     
      dfab  <-round(data.frame(indirect = x$ab[,j],boot = x$boot$mean[,j],sd=x$boot$sd[,j],
                           lower=x$boot$ci[1,1:niv],
                           upper=x$boot$ci[2,1:niv]),digits)
      rownames(dfab) <- rownames(x$ab)
      colnames(dfab)[1] <- dv[j]
      print(round(dfab,digits))
      }
      
        if(nmed > 1) {
    cat("\n 'ab' effects estimates for each mediator \n")
    for (j in 1:nmed) {
        dfab  <-round(data.frame(indirect = x$all.ab[,j],boot = x$boot$mean[,j+ndv],sd=x$boot$sd[,j+ndv],
                           lower=x$boot$ci[1,(j*niv +1):(j*niv +niv)],
                           upper=x$boot$ci[2,(j*niv +1):(j*niv +niv)]),digits)
      rownames(dfab) <- rownames(x$ab)
      colnames(dfab)[1] <- mv[j]
      print(round(dfab,digits))
      }
     
    
    }
}
    }  else {
    cat("\n\nEffect of interaction of ",iv[1], " with ", iv[2] , "  = ", round(x$direct[3],digits),"  S.E. = ", round(x$direct.reg$se[3,1],digits), " t  = ",round(x$direct.reg$t[3,1],digits), "  with p = ", signif(x$direct.reg$prob[3,1],digits))
    cat("\nIndirect effect due to interaction  of ",iv[1], " with ", iv[2] , "  = ", round(x$indirect,digits))
    cat("\nMean bootstrapped indirect interaction effect = ",round(x$boot$mean[1],digits), " with standard error = ",round(x$boot$sd[1],digits), " Lower CI = ",round(x$boot$ci.ab[1],digits), "   Upper CI = ", round(x$boot$ci.ab[2,i],digits))
    cat("\nSummary of a, b, and ab estimates and ab confidence intervals\n")
    } 
    }
    } else {#This is a pure regression  model, just show it
  summary(x)
   #  for(i in 1:ndv) {cat("\n DV = ",colnames(x$total.reg$beta)[i], "\n")
#           result.df <- data.frame( round(x$total.reg$beta[,i],digits),round(x$total.reg$se[,i],digits),round(x$total.reg$t[,i],digits),signif(x$total.reg$prob[,i],digits))
#               colnames(result.df) <- c("slope","se", "t", "p")              
#              print(result.df)
#             # cat("\nWith R2 = ", round(x$cpime.reg$R2[i], digits))
#               F <-  x$cprime.reg$df * x$cprime.reg$R2[i]/(((nrow(x$cprime.reg$beta)-1) * (1-x$cprime.reg$R2[i])))
#       pF <-  -expm1(pf(F,nrow(x$cprime.reg$beta),x$cprime.reg$df,log.p=TRUE)) 
#       cat("\nR =", round(sqrt(x$cprime.reg$R2[i]),digits),"R2 =", round(x$cprime.reg$R2[i],digits),  "  F =", round(F,digits), "on",nrow(x$cprime.reg$beta)-1, "and", x$total.reg$df,"DF   p-value: ",signif(pF,digits+1), "\n") 
#     
       }
     
      }
