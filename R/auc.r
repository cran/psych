#Written December 7, 2017 to show various decision processes.
 
"AUC" <- function(t=NULL,BR=NULL,SR=NULL,Phi=NULL,VP=NULL,labels=NULL,plot="b",zero=TRUE) {
col <- c("blue","red")
alpha <- .5
 col <- adjustcolor(col,alpha.f =alpha)     
if(!is.null(t) & is.null(BR)) {stopifnot(prod(dim(t)) ==4 || length(t) ==4)
  if(is.vector(t)) t <- matrix(t,2,	byrow=TRUE)
  } else {if(!is.null(t)) { phi <- SR
  SR <- BR
  BR <- t
 } 
  if(!is.null(Phi)) { 
  VP <- BR * SR + Phi *sqrt(BR * (1-BR) * SR *(1-SR))}
  t <- matrix(c(VP,SR -VP, BR-VP, 1- SR - BR + VP),2,2)
  }
 

  phi <- phi(t)
  tetra <- tetrachoric(t)
  
  colnames(t) <- c("Predicted.Pos","Predicted.Neg")
  rownames(t) <- c("True.Pos","True.Neg")
  observed <- t
  
    p1 <- t(t ( t)/rowSums(t))
    t <- t / sum(t)
   rs <- rowSums(t)
   cs <- colSums(t)
  p <- (t)/rs    #this has converted everything to conditional probabilities
  
 
  Sensitivity <- p[1,1]   #Sensitivity  VP
  Specificity <- p[2,2]    #VN
  FP <- 1 -  Specificity  #really? Yes, because we are talking conditional probabilities  
  PD <- rowSums(t)
  p1 <- p
  FN <- 1- Specificity   #really? Yes because it is conditional
  Accuracy <- Sensitivity *PD[1] + Specificity * (1-PD[1])
  phi <- phi(t)
 
  
 
   q <- qnorm(p1)   #based upon the conditional probabilities
   #spell these out for readability
   zVP <- q[1,1]
   zFP <- q[2,1]
    criterion <- -zFP
   
   zVN <- q[2,2]
  d.prime <-  zVP - zFP
  beta <- pnorm(d.prime - criterion)*rs[1]/(pnorm(criterion)*rs[2])
  xmax <- max(4,d.prime+3)
 x <- seq(-3,xmax,.1)
 
 y2 <- dnorm(x , -(!zero)*d.prime) * rs[2]    #noise
 y <- dnorm(x, zero * d.prime ) * rs[1]  #signal + noise
 
  if(plot == "b") {op <- par(mfrow=c(2,1))}
  if((plot =="b") | (plot =="a")) {
   plot(Sensitivity ~ FP,xlim=c(0,1),ylim=c(0,1),ylab="Valid Positives",xlab="False Positives",main="Valid Positives as function of False Positives")
  segments(0,0,FP,Sensitivity)
  segments(FP,Sensitivity,1,1)
  segments(0,0,1,1)}
 
 fpx <-  pnorm(x-qnorm(Specificity))
 vpx <-  pnorm(x+ qnorm(Sensitivity))
 fpx.diff <- diff(fpx)
 
 lower.sum <- sum(fpx.diff * vpx[-1])
 upper.sum <- sum(fpx.diff * vpx[-length(vpx)])
 auc <- (lower.sum + upper.sum)/2
  if((plot =="b") | (plot =="a")) { points(vpx ~ fpx,typ="l",lty="dashed")
   }
 
   if((plot =="b") | (plot =="d")) {
  plot(y ~ x, ylim=c(0,.4),ylab="Probability of observation",main="Decision Theory",type="l")
  points(y2 ~ x,lty="dashed",typ="l")
  
 # curve(dnorm(x,q[1,1]) * cs[1],-3,3,ylim=c(0,.4),ylab="Probability of observation",main="Decision Theory")
 # curve(dnorm(x,-q[2,2]) * cs[2], add=TRUE,lty="dashed")
  #segments((!zero)*(q[2,2]),0,(!zero)*(q[2,2]),dnorm(q[2,2])*cs[2])
  x1 <- x[x >= (criterion -(!zero) *d.prime)] 
  x1r <-rev(x1)
  y1 <- y[x >=  (criterion - (!zero) *d.prime)]
  y2c <- y2 [x >= (criterion -(!zero) *d.prime)]
  y1r <- rep(0,length(y1)) 
  polygon(c(x1,x1r),c(y1,y1r),col= col[1])
  polygon(c(x1,x1r),c(y2c,y1r),col= col[2])
  }
 if(plot =="b")  par(op)
 
 result<- list(observed=observed,probabilities=t,conditional=p1,q=q,Accuracy = Accuracy,Sensitivity=Sensitivity,Specificity=Specificity,AUC = auc,d.prime = d.prime,beta = beta, criterion=criterion,phi=phi,tetrachoric=tetra$rho)
class(result) <- c("psych","auc")
return(result)
   }
   
print.psych.auc <- function(x,digits=2) {
cat('Decision Theory and Area under the Curve\n')
cat('\nThe original data implied the following 2 x 2 table\n')
print(x$probabilities,digits=digits)
cat('\nConditional probabilities of \n') 
print(x$conditional,digits=digits)
cat('\nAccuracy = ',round(x$Accuracy,digits=digits),' Sensitivity = ',round(x$Sensitivity,digits=digits), '  Specificity = ',round(x$Specificity,digits=digits),  '\nwith Area Under the Curve = ', round(x$AUC,digits=digits) )
cat('\nd.prime = ',round(x$d.prime,digits=digits), ' Criterion = ',round(x$criterion,digits=digits), ' Beta = ', round(x$beta,digits=digits))
cat('\nObserved Phi correlation = ',round(x$phi,digits=digits), '\n Inferred latent (tetrachoric) correlation  = ',round(x$tetrachoric,digits=digits))
}


