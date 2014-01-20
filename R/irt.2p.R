"irt.2p" <- 
function(delta,beta,items) {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
 irt.2par <- function(x,delta,beta,scores) {
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  
  

 num <- dim(items)[1]
 fit <- matrix(NaN,num,2)
 total <- rowMeans(items,na.rm=TRUE)
 count <- rowSums(!is.na(items))
 
 for (i in 1:num) {
 	
 	if (count[i]>0)  {myfit <- optimize(irt.2par,c(-4,4),beta=beta,delta=delta,scores=items[i,]) #how to do an apply?
     		fit[i,1] <- myfit$minimum    
     		fit[i,2] <- myfit$objective  #fit of optimizing program
     		} else {
         	fit[i,1] <- NA
  			fit[i,2] <- NA 
    			}    #end if else
    }  #end loop 
	irt.2p <-data.frame(total,theta=fit[,1],fit=fit[,2],count)}
