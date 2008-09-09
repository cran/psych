#thurstonian scaling
#Thurstone case V  (assumption of equal and uncorrelated error variances)
#Version of March 22, 2005
#revised April 4, 2008
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

"thurstone" <- 
function(x, ranks = FALSE,digits=2) {     #the main routine
    if (ranks) {choice <- choice.mat(x)  
       } else {if (is.matrix(x)) choice <- x
              choice <- as.matrix(x)}
            
        scale.values <- scale.vect(choice)  #convert rank order information to choice matrix 
     
       fit <- thurstone.fit(scale.values,choice,digits)        #see how well this model fits
       thurstone <- list("scale"=round(scale.values,2), "GF"= fit$Fit, "residual"=fit$residual,"data"= choice)
    return(thurstone)}
    
#the functions used by thurstone (local to thurstone)

 #this next function does almost all of the work, by transforming choices to normal deviates
 scale.mat <- function(x) {scale.mat <- qnorm(x)}   #convert choice matrix to normal deviates
 
 #find the average normal deviate score for each item, and subtract from the minimum
 scale.vect <- function(x)           #do the actual thurstonian scaling
     {nvar <- dim(x)[2]
     score <- colSums(scale.mat(x))   #find the summed scores
     minscore <-min(score)            #rescale to 0 as minumum
     relative <- score - minscore
     relative <- relative/nvar
     return(relative) }         #the thurstone scale values
     
 #next two functions are used for finding goodness of fit   -- not a very powerful measure  
 thurstone.model <- function (x,square) {
    if (!square) {      #returns the lower diagonal distances
    return(1-pnorm(dist(x,diag=TRUE))) } else {
    return(1-pnorm(dist(x,diag=TRUE,upper=TRUE)))}
   }
    
 thurstone.fit <- function (s,d,digits) {   #s is the thurstone scale values, d is the data matrix
     model <- thurstone.model(s,FALSE)
     dm <- as.dist(d)              #makes the distance matrix a lower diagonal
     error <- model - dm         #error = model - data  or more conventionally data = model+error
     fit1 <- 1- sum(error*error)/sum(dm*dm)   #a typical fit measure -- squared error/squared original
    
     model <- thurstone.model(s,TRUE)
     error <- model - dm  
     fit2 <- 1- sum(error*error)/sum(d*d)   #a typical fit measure -- squared error/squared original
     thurstone.fit <- list("Fit"=c(round(fit1,digits),round(fit2,digits)),"model"=round(model,2),"residual"=round(error,2))
     
     }
     
#if we have rank order data, then convert them to a choice matrix
#convert a rank order matrix into a choice matrix
orders.to.choice <- function(x,y) {      #does one subject (row) at a time
   nvar <-dim(x)[2]
   for (i in 1:nvar) {
      for (j in 1:nvar) {
        if (x[j]< x[i] ) {y[i,j] <- y[i,j]+1
         } else if  (x[j] == x[i]) {
                       y[j,i] <- y[j,i]+.5    #take ties into account
                       y[i,j] <- y[i,j] + .5
              } else  y[j,i] <- y[j,i]+1        
          }} 
          return(y)}
          
#  makes repeated calls to orders.to.choice  -- can we vectorize this?        
choice.mat <-function(x) {
   nsubs<- dim(x)[1]
   nvar <- dim(x)[2]
   y <- matrix(0,ncol=5,nrow=5)
   for (k in 1:nsubs) {
      y <-orders.to.choice(x[k,],y) }     #is there a way to vectorize this?
     d <- diag(y)    
     y<- y/(2*d)
     lower <- 1/(4*nsubs)     #in case of 0 or 1, we put limits 
     upper <- 1- lower
     y[y<lower] <-lower    
     y[y>upper] <- upper
    return(y) }
    
