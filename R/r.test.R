 "r.test" <- 
 function(n,r12, r34=NULL, r23=NULL,r13=NULL,r14=NULL,r24=NULL,n2=NULL,pooled=TRUE, twotailed=TRUE) {
  cl <- match.call()
  if(is.null(r34) & is.null(r13) & is.null(r23)) {  #test for significance of r
     
  
     t <- r12*sqrt(n-2)/sqrt(1-r12^2) 
     p <- 1-pt(abs(t),n-2) 
     if(twotailed) p <- 2*p
     ci <- r.con(r12,n)
      result <-  list(Call=cl,Test="Test of significance of a  correlation",t=t,p=p,ci=ci)
  } else {if(is.null(r23)) { #compare two independent correlations
        xy.z <- 0.5*log((1+r12)/(1-r12))
        xz.z <- 0.5*log((1+r34)/(1-r34))
        if(is.null(n2)) n2 <- n
        se.diff.r <- sqrt(1/(n-3) + 1/(n2-3))
        diff <- xy.z - xz.z
        z <- abs(diff/se.diff.r)
         p <- (1-pnorm(z ))
          if(twotailed) p <- 2*p
      result <-  list(Call=cl,Test="Test of difference between two independent correlations",z=z,p=p)
                             }  else { if (is.null(r14)) {#compare two dependent correlations case 1
      
        
        #here we do two tests of dependent correlations
       #figure out whether correlations are being specified by name or order
        if(!is.null(r34)) {if(is.null(r13)) {r13 <- r34} }
       diff <- r12-r13
       determin=1-r12*r12 - r23*r23 - r13*r13 + 2*r12*r23*r13
       av=(r12+r13)/2
       cube= (1-r23)*(1-r23)*(1-r23)
       t2 = diff * sqrt((n-1)*(1+r23)/(((2*(n-1)/(n-3))*determin+av*av*cube)))
       p <- pt(abs(t2),n-3,lower.tail=FALSE)  #changed to n-3 on 30/11/14
        if(twotailed) p <- 2*p
        #the call is ambiguous, we need to clarify it
        cl <- paste("r.test(n = ",n, ",  r12 = ",r12,",  r23 = ",r23,",  r13 = ",r13, ")")
      result <- list(Call=cl,Test="Test of difference between two correlated  correlations",t=t2,p=p)                                  
      } else { #compare two dependent correlations, case 2
       z12 <- fisherz(r12)
    z34 <- fisherz(r34)
    pooledr <- (r12+r34)/2	
    if (pooled) { r1234=  1/2 * ((r13 - pooledr*r23)*(r24 - r23*pooledr) + (r14 - r13*pooledr)*(r23 - pooledr*r13)   +(r13 - r14*pooledr)*(r24 - pooledr*r14) + (r14 - pooledr*r24)*(r23 - r24*pooledr))
               z1234 <- r1234/((1-pooledr^2)*(1-pooledr^2))} else {
 
     r1234=  1/2 * ((r13 - r12*r23)*(r24 - r23*r34) + (r14 - r13*r34)*(r23 - r12*r13)   +(r13 - r14*r34)*(r24 - r12*r14) + (r14 - r12*r24)*(r23 - r24*r34))
     z1234 <- r1234/((1-r12^2)*(1-r34^2))}
      ztest  <- (z12-z34)* sqrt(n-3) /sqrt(2*(1-z1234))
      z <- ztest
       p <- (1-pnorm(abs(z) ))
          if(twotailed) p <- 2*p
       result <-  list(Call=cl,Test="Test of difference between two dependent correlations",z=z,p=p)
      
                              }
           }
    } 
    class(result) <- c("psych", "r.test")
    return(result)
   }