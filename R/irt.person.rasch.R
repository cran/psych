#steps towards an IRT program
#we find the difficulties using ir.item.diff.rasch
#now estimate the thetas 
#Then, to find the person parameters, use optimize 

"irt.person.rasch" <- 
function(diff,items) {
#
#the basic one parameter model
 irt <- function(x,diff,scores) {
  fit <- -1*(log(scores/(1+exp(diff-x)) + (1-scores)/(1+exp(x-diff))))
mean(fit,na.rm=TRUE)
  }
 #
 diff<- diff
 items <-items
 num <- dim(items)[1]
 fit <- matrix(NA,num,2)
 total <- rowMeans(items,na.rm=TRUE)
 count <- rowSums(!is.na(items))
 
 for (i in 1:num) {
 	
 	if (count[i]>0)  {myfit <- optimize(irt,c(-4,4),diff=diff,scores=items[i,]) #how to do an apply?
     		fit[i,1] <- myfit$minimum    
     		fit[i,2] <- myfit$objective  #fit of optimizing program
     		} else {
         	fit[i,1] <- NA
  			fit[i,2] <- NA 
    			}    #end if else
    }  #end loop 
irt.person.rasch <-data.frame(total,theta=fit[,1],fit=fit[,2],count)}