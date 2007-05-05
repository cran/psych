#steps towards an IRT program
#we find the difficulties using irt.item.diff.rasch
#now estimate the thetas 
#Then, to find the person parameters, use optimize 

"irt.1p" <- 
function(delta,items) {
#
#the basic one parameter model (aka a Rasch model)
 irt <- function(x,delta,scores) {
  fit <- -1*(log(scores/(1+exp(delta-x)) + (1-scores)/(1+exp(x-delta))))
mean(fit,na.rm=TRUE)
  }
 #
 delta<- delta
 items <-items
 num <- dim(items)[1]
 fit <- matrix(NA,num,2)
 total <- rowMeans(items,na.rm=TRUE)
 count <- rowSums(!is.na(items))
 
 for (i in 1:num) {
 	
 	if (count[i]>0)  {myfit <- optimize(irt,c(-4,4),delta=delta,scores=items[i,]) #how to do an apply?
     		fit[i,1] <- myfit$minimum    
     		fit[i,2] <- myfit$objective  #fit of optimizing program
     		} else {
         	fit[i,1] <- NA
  			fit[i,2] <- NA 
    			}    #end if else
    }  #end loop 
irt.1p <-data.frame(total,theta=fit[,1],fit=fit[,2],count)}