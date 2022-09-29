#thurstonian scaling
#Thurstone case V  (assumption of equal and uncorrelated error variances)
#Version of March 22, 2005
#revised April 4, 2008
#Revised March 29, 2009 to be much cleaner code
#Do a Thurstone scaling (Case 5) of either a square choice matrix or a rectangular rank order matrix
#Need to add output options to allow users to see progress
#data are either a set of rank orders or
#a set of of choices (columns chosed over row)
#output is a list of 
# scale scores
#goodness of fit estimate    (1 - sumsq(error)/sumsq(original)   
# the model choice matrix
#the error (residual) matrix
#the original choice matrix 
#modified 5/10/22 to handle case of missing data in ranks

"thurstone" <- 
function(x, ranks = FALSE,digits=2) {     #the main routine
 cl <- match.call()
    if (ranks) {choice <- choice.mat(x)   #convert rank order information to choice matrix 
       } else {if (is.matrix(x)) choice <- x
              choice <- as.matrix(x)}
        
        scale.values <- colMeans(qnorm(choice)) - min(colMeans(qnorm(choice)))  #do the Thurstonian scaling 
        model <- pnorm(-scale.values %+% t(scale.values))
		error <- model - choice
		fit <- 1-(sum(error*error)/sum(choice*choice))
       result <- list(scale=round(scale.values,digits), GF= fit, choice=choice,residual=error,Call=cl)
       class(result) <- c("psych","thurstone")
    return(result)}
    
#the functions used by thurstone (local to thurstone)


     
#if we have rank order data, then convert them to a choice matrix
#convert a rank order matrix into a hoice matrix
#needs to handle missing data  fixed 5/10/22
orders.to.choice <- function(x,y) {      #does one subject (row) at a time
   nvar <- length(x)          #changed to length 11/8/20
   
   for (i in 1:nvar) {
      for (j in 1:nvar) {
      if(!is.na(x[j]) & !is.na(x[i])) {
        if (x[j]< x[i] ) {y[i,j] <- y[i,j]+1
         } else if  (x[j] == x[i]) {
                       y[j,i] <- y[j,i]+.5    #take ties into account
                       y[i,j] <- y[i,j] + .5
              } else  y[j,i] <- y[j,i]+1        
          }}} 
          return(y)}
          
          
#  makes repeated calls to orders.to.choice  -- can we vectorize this?        
choice.mat <-function(x) {
   nsubs<- dim(x)[1]
   nvar <- dim(x)[2]
   y <- matrix(0,ncol=nvar,nrow=nvar)
   for (k in 1:nsubs) {
      y <-orders.to.choice(x[k,],y) }     #is there a way to vectorize this?
     d <- diag(y)    
     y<- y/(2*d)
     lower <- 1/(4*nsubs)     #in case of 0 or 1, we put limits 
     upper <- 1- lower
     y[y<lower] <-lower    
     y[y>upper] <- upper
    return(y) }
    
    
#irt type data
#subjects endorse or do not endorse an item  
item.to.choice <- function(x) {
    nsubs<- dim(x)[1]
   nvar <- dim(x)[2]
   count <- t(x) %*% x
   dx <- diag(count)
   y <- dx - count
   diag(y) <- dx
   y <- y/nsubs
   
   for (k in 1:nsubs) {
     temp <- x[k,]
        for (i in 1:nvar) {
        for (j in 1:nvar) {
        if(temp[j]>0) count [i,j] <- count[i,j]+ temp[i]
          } 
          }
          
         }     #is there a way to vectorize this?
     d <- diag(y)    
     y<- y/(2*d)
     lower <- 1/(4*nsubs)     #in case of 0 or 1, we put limits 
     upper <- 1- lower
     y[y<lower] <-lower    
     y[y>upper] <- upper
    return(y) }
    


#do hare ranking    
hare <- function(x, n=1) {
#loop over candidates
#find the person(s) with the fewest 1 place votes

while(NCOL(x) > n) {
#for(iter in (1:(NCOL(x)-n))) {  #continue until n +1 are left
count.first <- colSums(apply(x,2,function(x) x==1), na.rm=TRUE)
min.count <- min(count.first)
which.min <- which(count.first == min.count)
which.min <- which.min[1] #just do one candidate at a time
#reallocate their votes
who.min <- which(x[,which.min]==1, arr.ind=TRUE)
#the row is the person we want 
if(NROW(who.min )> 0) {
if(!is.null(dim(who.min))) who.min <- who.min[,1]

for(i in (1:length(who.min)))  {
 x[who.min[i],] <- x[who.min[i],] -1} 
 }  #reallocate their votes
x <- x[,-which.min,drop=FALSE]   #drop the lowest vote getter
   #do it again until the NCOL = the number we want

 }
 if(NCOL(x) > 1){
 count.first <- colSums(apply(x,2,function(x) x==1), na.rm=TRUE) 
  
 result <- sort(count.first,decreasing=TRUE)[1:n]} else {result <- colSums(x==1,na.rm=TRUE)
}
 
return(result)}   #show the winners



#do hare ranking    a different way  
hare2 <- function(x, n=1) {
#loop over candidates
#find the person(s) with the fewest 1 place votes

while(NCOL(x) > n) {
#for(iter in (1:(NCOL(x)-n))) {  #continue until n +1 are left

temp <- apply(x,1,function(x) which(x==1))
tt <- table(unlist(temp))
min.tt <- min(tt)
which.min <- which(tt==min.tt)
who.min <- which(x[,which.min]==1, arr.ind=TRUE)

#the row is the person we want 
if(NROW(who.min )> 0) {
if(!is.null(dim(who.min))) who.min <- who.min[,1]

for(i in (1:length(who.min)))  {
 x[who.min[i],] <- x[who.min[i],] -1} 
 }  #reallocate their votes
x <- x[,-which.min,drop=FALSE]   #drop the lowest vote getter
   #do it again until the NCOL = the number we want

 }
 if(NCOL(x) > 1){
 count.first <- colSums(apply(x,2,function(x) x==1), na.rm=TRUE) 
  
 result <- sort(count.first,decreasing=TRUE)[1:n]} else {result <- colSums(x==1,na.rm=TRUE)
}
 
return(result)}   #show the winn
 
# count first place votes

    
