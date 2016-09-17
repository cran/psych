#revised August 31, 2012 to allow for estimation of all 0s or all 1s
#modified March 10, 2016 to allow for quasi Bayesian estimates using normal theory
#which doesn't seem to do what I want to do, so we are not doing it
#June 30-July 7 , 2016  Trying to make it faster by parallelizing the code and
#reducing the number of items to score when using keys.
#uses local minima -- problematic for small number of items 
#corrected various problems with finding total (sum) scores  7/4/16
#starting to make parallel for speed
#seems to result in at least a factor of 2 improvement
#when using stats from fa, the discrim parameter are not necessarily in the order from a keyslist
#this is only a problem if selecting items from a larger set of items
#fixed this August 21, 2016
#the irt.2 function (dichotomous items) iis much slower than the polytomous solution
#probably because we took parallelization one step too far
#I have not removed that extra level

#
#the function to do 2  parameter dichotomous IRT
"score.irt.2" <- 
function(stats,items,keys=NULL,cut=.3,bounds=c(-4,4),mod="logistic") {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
#This does the normal fit
#has several internal functions

irt.2par.norm <-  function(x,delta,beta,scores) {
  fit <- -1*(log(scores*(1-pnorm(beta*(delta-x))) + (1-scores)*(1-pnorm(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  } 
#This does the logistic fit 
 irt.2par <- function(x,delta,beta,scores) { 
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  
###  
##the next two are the parallelized functions
#parallelize by subject seems most helpful
bySubject <- function(i,count.f,total.f,items.f,discrim.f,diffi.f) {

 	#First we consider the case of all right or all wrong
 	#but we also need to consider the person with no data!
 	if (count.f[i] > 0) {if((sum(items.f[i,],na.rm=TRUE) ==0 )  | (prod(items.f[i,],na.rm=TRUE) ==  1 )) { 
 	if(sum(items.f[i,],na.rm=TRUE) ==0 ) {
 		if(mod=="logistic") {p <- log(1-(1-items.f[i,])/(1+exp(discrim.f*(diffi.f))) )} else {  #logistic
 				p <- log(1-(pnorm(items.f[i,]*discrim.f*diffi.f))) }  #normals
  	 pall <- exp(sum(p,na.rm=TRUE))
   
   #	theta[i] <- qnorm(pnorm(qnorm(pall))/2)  #the z value of 1/2 the quantile value of pall
   theta <- qnorm(pnorm(qnorm(pall)))  #the z value of the quantile value of pall
   	fit <- 0 
 	 # cat ("\nThe case of all wrong",i,theta[i])
 	} else {  #the case of all right	
 		if(mod == "logistic") {	p <- log((items.f[i,])/(1+exp(discrim.f*(diffi.f))) )} else  { 
 	                        	p <- log((items.f[i,])*(1 - pnorm(1- discrim.f*(diffi.f)) ))  } 
 	
 	 pall <- exp(sum(p,na.rm=TRUE))
  # theta <- qnorm(1-pnorm(qnorm(pall))/2)  #the z value of 1/2 the quantile value of pall
   theta <- qnorm(1-pnorm(qnorm(pall))) #or, perhaps just hte z value of the quantile value of pall
   fit <- 0 
       }
 	} else { #cat("the normal case",i )
 	
 	 	beta=discrim.f[!is.na(items.f[i,])]
 	   delta=diffi.f[!is.na(items.f[i,])]
 	   scores=t(items.f[i,!is.na(items.f[i,])]) #make this numeric in the case of weird (highly missing) data
 	if(mod=="logistic") {
  myfit <- optimize(irt.2par,bounds,beta=beta,delta=delta,scores=scores)  #how to do an apply?
  #myfit <- optimize(irt.2par,bounds,beta=discrim.f,delta=diffi.f,scores=items.f[i,])
    } else {myfit <- optimize(irt.2par.norm,bounds,beta=beta,delta=delta,scores=scores)} #do a normal fit function   # do this one as well
     		
     		theta <- myfit$minimum    
     		fit <- myfit$objective  #fit of optimizing program 
     		}} else  {#cat("\nno items for subject",i)
             total.f[i]  <- NA
     		 theta <- NA
  		     fit  <- NA 
 	                    
 	        }  #end if count ... else  
 	   return(list(theta,total.f[i],fit) )                            	          
    }  #end bySubject
 
 #parallelize by factor   
#this is the the one to use when parallelized  
bigFunction <- function(f,n.obs,stats,items,keys=NULL,cut=.3,bounds=c(-5,5),mod="logistic") {
 nf <- length(stats$difficulty)
diff <- stats$difficulty[[f]]
cat <- dim(diff)[2]

if(nf < 2) {#discrim <- drop(stats$discrimination)
            discrim <- stats$discrimination  # although I need to check this with keys 
             if(!is.null(keys)) {discrim <- discrim * abs(keys)}
			 } else {discrim <- stats$discrimination[,f] 
			if(!is.null(keys)) {discrim <- discrim * abs(keys[,f])
                    }}

###

fit <- rep(NA,n.obs)
theta <- rep(NA,n.obs)


if(is.null(keys)) {#drop the items with discrim < cut
                       
	items.f <- items[,(abs(discrim[,f]) > cut) ,drop=FALSE]  #get rid of the those items that are not being processed for this factor
	diffi.f <- diff[(abs(discrim[,f]) > cut)]   #and the redundant diffi
	discrim.f <- discrim[(abs(discrim[,f]) >cut),drop=FALSE ]  #and get rid of the unnecessary discrim values

     } else { #the  case of scoring with a keys vector
                   
	items.f <- items[,(abs(keys[,f]) > 0) ,drop=FALSE]  #get rid of the those items that are not being processed for this factor
	discrim.f <- discrim[(abs(keys[,f]) > 0),drop=FALSE ]  #and get rid of the unnecessary discrim values
	diffi.f <- diff[(abs(keys[,f]) > 0)]   #and the redundant diffi
}

diffi.vect <- as.vector(t(diffi.f)) 

#discrim.F.vect <- rep(discrim.f,each=cat)
#discrim.f <- discrim[(abs(discrim > cut)),drop=FALSE] 
discrim.F.vect <- as.vector(t(discrim.f))

if(is.matrix(discrim)) discrim.F.vect <- drop(discrim.F.vect)
 total <- rowMeans(t(t(items.f)*sign(discrim.F.vect)),na.rm=TRUE)
 count <- rowSums(!is.na(items.f))
#We can speed this up somewhat if we don't try to fit items with 0 discrim (i.e., items that do not load on the factor or are not keyed)

#do it for all subject for this factor
#now, lets parallelize this for each subject as well
#this is probably a bad idea, for it leads to an amazing amount of overhead in terms of memory and processes
#mapply for debugging,  mcmapply for parallel
#lets just try mapply to see if it gets around the problem
#actually, making this mcmapply and the call to bigFunction mapply seems to be the solution
#especially when we are doing scoreIrt.1pl or scoreIrt.2pl which is already doing the parallelsim there

subjecttheta <-mcmapply(bySubject,c(1:n.obs),MoreArgs = list(count,total,items.f,discrim.f,diffi.f))  #returns a list of theta and fit

subjecttheta <- matrix(unlist(subjecttheta),ncol=3,byrow=TRUE)
theta <- subjecttheta[,1]
total <- subjecttheta[,2]
fit <- subjecttheta[,3]
    theta [theta < bounds[1]] <- bounds[1]
    theta[theta > bounds[2]] <- bounds[2]
   #  if((!is.null(keys)) & (all(keys[,f] == -1) || (sign(cor(discrim,keys[,f],use="pairwise")) < 0) )) {theta <- -theta
   #if((!is.null(keys)) & (all(keys[,f] == -1)  )) {theta <- -theta
    #                             total <- -total}
            
  
 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]
 nvar <- dim(items)[2]  
# scores <- matrix(NaN,nrow=n.obs,ncol=nf*3)

#scores <- list(nf*3)

scores <- list(theta,total,fit)
return(scores)
} # end of bigFunction  

#now do the scoring one factor at a time but allowing multiple cores

#we now start score.irt.2 proper
#this finds scores using multiple cores if they are available

 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]

#this parallels by factor  which in turn is parallelized by subject in bySubject
#use mapply for debugging, mcmapply for parallel processing
#since we are already parallelizing by scale when we call scoreIrt.1pl or .2pl, this is not necessary to parallelize
scores <-  mapply(bigFunction,c(1:nf),MoreArgs=list(n.obs=n.obs,items=items,stats=stats,keys=keys, cut=cut, bounds=bounds, mod=mod))

nf <- length(stats$difficulty)
scores <- matrix(unlist(scores),ncol=nf*3)
scores <- scores[,c(seq(1,nf*3,3),seq(2,nf*3+1,3),seq(3,nf*3 +2,3))]
colnames(scores) <- paste(rep(c("theta","total","fit"),each=nf),1:nf,sep="")
return(scores)  
                    
}#end of score.irt.2
#############

very.close <- function(x,y,tolerance = .Machine$double.eps) {
abs(x-y) < tolerance}

"score.irt.poly" <- 
function(stats,items,keys=NULL,cut=.3,bounds=c(-4,4),mod="logistic") {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
#created July 4, 2011
 irt.2par.poly <- function(x,delta,beta,scores) {
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  
####The function that is  parallelized
big.poly <- function(f,n.obs,stats,items,keys=NULL,cut=.3,bounds=c(-5,5),mod="logistic") {
nf <- ncol(stats$discrimination)
#for (f in 1:nf) { #do it for every factor/scale  
diff <- stats$difficulty[[f]]

if(nf < 2) {discrim <- stats$discrimination
              if(!is.null(keys)) {discrim <- discrim * abs(keys)
                              } } else {discrim <- stats$discrimination[,f] 
                                if(!is.null(keys)) {discrim <- discrim * abs(keys[,f]) 
                  }  
                    }
cat <- dim(diff)[2]
 total <- rep(NA,n.obs)
 fit <- rep(NA,n.obs)
 theta <- rep(NA,n.obs)
 item.f <- t(items)
 item.f[abs(discrim) < cut] <- NA  #this does not change the item, just the temp version of the item
 item.f <- t(item.f)
 
###
if(!is.null(keys)) {item.f <- item.f[,(abs(keys[,f] )> 0) ,drop=FALSE]  #get rid of the those items that are not being processed for this factor
discrim.f <- discrim[(abs(keys[,f]) > 0),drop=FALSE ]  #and get rid of the unnecessary discrim values
#diffi.f <- diff[(abs(keys[,f]) > 0),drop=FALSE]   #and the redundant diffi
diffi.f <- diff[(abs(keys[,f]) > 0),byrows=TRUE]

#diffi.vect <- as.vector(t(diffi.f)) 
#diffi.vect <- as.vector(diffi.f) 
diffi.vect <- as.vector(t(diff[(abs(keys[,f]) > 0),byrows=TRUE]))
discrim.F.vect <- rep(discrim.f,each=cat)
} else { discrim.f <- discrim
diffi.f <- diff
diffi.vect <- as.vector(t(diff) )
discrim.F.vect <- rep(discrim.f,each=cat)
}

## notice that this is vectorized and does it for all subjects 
#seem to have solved the problem of missing items which are reversed.

 total <- rowMeans(t(t(item.f )* as.vector(sign(discrim.f))),na.rm=TRUE) #fixed 11/11/11 to be as.vector

# total.positive <- rowMeans(t(t(item.f)* as.vector(sign(discrim.f) > 0)),na.rm=TRUE)
# total.negative <- rowMeans(t(t(item.f)* as.vector(sign(discrim.f) < 0)),na.rm=TRUE)
##

 num.keyed <- rowSums(!is.na(item.f))
 num.reversed <- rowSums(!is.na(item.f[,discrim.f < 0,drop=FALSE]))

 total <- total + num.reversed * (max.item- min.item+1)/num.keyed  + min.item 
 total[is.nan(total)] <- NA
 count <- rowSums(!is.na(item.f))


#but now, we need to the next step one at a time (I think)
 for (subj in 1:n.obs) {	
 	if (count[subj]> 0)  { #just do cases where we actually have data
 	       newscore <-  NULL
 	       score <- item.f[subj,] #just the items to be scored
 	        if((very.close(total[subj],min.item)) | (very.close(total [subj],(max.item+min.item)) )){ # first check for all lowest responses or all highest responses  
 	           	 if(very.close(total[subj],min.item)) { # The case of all wrong
 	            	if(mod=="logistic") { p <- log(1-1/(1+exp(abs(discrim.f[!is.na(item.f[subj,])])*diffi.f[!is.na(item.f[subj,]),-1])))} else {  #logistic   
 					p <- log(1-(pnorm(1- abs(discrim.f[!is.na(item.f[subj,])])*diffi.f[!is.na(item.f[subj,]),-1])) ) }  #normals
  	        pall <- exp(sum(p,na.rm=TRUE)) 
  	         theta[subj] <- qnorm(pnorm(qnorm(pall))) 
  	         fit[subj] <- 0} else {
 	       if(very.close(total [subj],(max.item+min.item))) {#all right
 	            if(mod == "logistic") {	 p <- log(1/(1+exp(abs(discrim.f[!is.na(item.f[subj,])])*diffi.f[!is.na(item.f[subj,]),])))} else  { 
 	                        	p <- log((1)*(1 - pnorm(1- abs(discrim.f[!is.na(item.f[subj,])])*diffi.f[!is.na(item.f[subj,]),])) )  } 

 	           pall <- exp(sum(p,na.rm=TRUE))
              theta[subj] <- qnorm(1-pnorm(qnorm(pall)))  #the z value of  the quantile value of pall
            fit[subj] <- 0 } 
 	       }} else {  #just process those items where we have some responses

           for (i in 1:ncol(item.f)) {  #Treat the items as a series of 1 or 0 responses  - but note that cat = max - min
                 if(is.na(score[i])) {newscore <- c(newscore,rep(NA,cat)) } else {
                   if(very.close(score[i],( cat))) {newscore <- c(newscore,rep(1,cat))
                    } else {
                        newscore <- c(newscore,rep(1,score[i]),rep(0,cat-score[i])) }
                 
                      }}
 	        myfit <- optimize(irt.2par.poly,bounds,beta=discrim.F.vect,delta=diffi.vect,scores=newscore) 
     		theta[subj] <- myfit$minimum    
     		fit[subj] <- myfit$objective  #fit of optimizing program
 
     		}} else {
         	fit[subj] <- NA
  			theta[subj] <- NA 
    			}    #end if else
    } 
   
   if((!is.null(keys)) & (all(keys[,f] == -1) )) {theta <- -theta
                                 total <- -total}
    theta[theta < bounds[1]] <- bounds[1]
    theta[theta > bounds[2]] <- bounds[2]
    
    scores <- list(theta,total, fit)
    return(scores)
    
    } #end of big function
    

##the start of the irt.poly.function after setting up the various subfunctions
 min.item <- min(items,na.rm=TRUE)  #uses local minima --probably  problematic for small number of items
 items <- items - min.item #this converts scores to positive values from 0  up 
 max.item <- max(items,na.rm=TRUE)  #we use this when reverse score -- but note that this is not the original max value.  We will adjust this in total
 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]
 nvar <- dim(items)[2]

#mcmapply for parallel, mappy for debugging
scores   <-  mcmapply(big.poly,1:nf,MoreArgs=list(n.obs=n.obs,stats=stats,items=items,keys=keys,cut=.3,bounds=bounds,mod=mod))
 

scores <- matrix(unlist(scores),ncol=nf*3)
scores <- scores[,c(seq(1,nf*3,3),seq(2,nf*3+1,3),seq(3,nf*3 +2,3))]
colnames(scores) <- paste(rep(c("theta","total","fit"),each=nf),1:nf,sep="")

return(scores)
} #end of score.irt.poly 

"score.irt" <- function(stats=NULL,items,keys=NULL,cut=.3,bounds=c(-4,4),mod="logistic")	{
   message("score.irt has been replaced by scoreIrt, please change your call")
   scoreIrt(stats=stats,items=items,keys=keys,cut=cut,bounds=bounds,mod=mod) }
   
"scoreIrt" <- function(stats=NULL,items,keys=NULL,cut=.3,bounds=c(-4,4),mod="logistic")	  {
#added the tau option in switch in case we have already done irt.tau   6/29/16
#we need to adjust the discrimination order from irt.fa to match the order of the items
if(!is.null(keys) && is.list(keys)){   select <- sub("-","",unlist(keys))
              items <- items[select] 
              keys <- make.keys(items,keys)}
              
if(!is.null(keys) && (is.vector(keys)))  keys <- matrix(keys)
if (length(class(stats)) > 1) {
    if(!is.null(keys) && is.vector(keys)) keys <- as.matrix(keys)
    switch(class(stats)[2],
    irt.poly = {scores <- score.irt.poly(stats$irt,items,keys,cut,bounds=bounds,mod=mod) },
    irt.fa =   {scores <- score.irt.2(stats$irt,items,keys,cut,bounds=bounds,mod=mod)},
    fa = {tau <- irt.tau(items)  #this is the case of a factor analysis to be applied to irt
          nf <- dim(stats$loadings)[2]
          diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-stats$loadings[,i]^2) }
     
         discrim <- stats$loadings/sqrt(1-stats$loadings^2)
   	class(diffi) <- NULL
   	class(discrim) <- NULL
   	new.stats <- list(difficulty=diffi,discrimination=discrim)
    scores <- score.irt.poly(new.stats,items,keys,cut,bounds=bounds)},
    
    tau = {tau <- stats   #get the tau stats from a prior run
         if(is.matrix(keys)) {nf <- dim(keys)[2]} else {nf <-1}
         diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau }
         discrim <- keys
  		 class(diffi) <- NULL
   		class(discrim) <- NULL
   		new.stats <- list(difficulty=diffi,discrimination=discrim)
   if(dim(tau)[2] ==1)  {scores <- score.irt.2(stats=new.stats,items=items,keys=keys,cut=cut,bounds=bounds)} else {
                         scores <- score.irt.poly(stats=new.stats,items=items,keys=keys,cut=cut,bounds=bounds)}
                         }   
    )
       #we should have a null case
    } else {#input is a keys matrix
         tau <- irt.tau(items)  #this is the case of a using a scoring matrix to be applied to irt
         if(is.matrix(keys)) {nf <- dim(keys)[2]} else {nf <-1}
         diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau }
         if(!is.null(keys)) {discrim <- keys} else {stop("I am sorry, you specified  tau  but not  keys.")}
   class(diffi) <- NULL
   class(discrim) <- NULL
   new.stats <- list(difficulty=diffi,discrimination=discrim)
   if(dim(tau)[2] ==1)  {scores <- score.irt.2(stats=new.stats,items=items,keys=keys,cut=cut,bounds=bounds)} else {
                         scores <- score.irt.poly(stats=new.stats,items=items,keys=keys,cut=cut,bounds=bounds)}
                         }
scores <- data.frame(scores)
if(!is.null(keys)) {colnames(scores) <-c( paste(colnames(keys),"theta",sep="-"),paste(colnames(keys),"total",sep="-"),paste(colnames(keys),"fit",sep="-"))}
return(scores)
}


#find tau from dichotomous or polytomous data without bothering to find the correlations
#useful for score.irt
#modified July 14, 2016 to speed up significantly by dropping the xt <- table(x)  line

"irt.tau" <- function(x) {  

x <-as.matrix(x)
nvar <- dim(x)[2]
xmin <- min(x,na.rm=TRUE)
xmax <- max(x,na.rm=TRUE)
nvalues <- xmax-xmin +1
if(nvalues > 10) stop("You have more than 10 categories for your items, polychoric is probably not needed")
#xt <- table(x)   #this can take a long time for sapa data

#nvalues <- length(xt)  #find the number of response alternatives 
if(nvalues ==2) {tau <- -qnorm(colMeans(x,na.rm=TRUE))
   tau <- as.matrix(tau)
   rownames(tau) <- colnames(x)}  else {
if(nvalues > 10) stop("You have more than 10 categories for your items, polychoric is probably not needed")
#xmin <- min(x,na.rm=TRUE)
xfreq <- apply(x- xmin+ 1,2,tabulate,nbins=nvalues)
n.obs <- colSums(xfreq)
xfreq <- t(t(xfreq)/n.obs)
tau <- qnorm(apply(xfreq,2,cumsum))[1:(nvalues-1),]  #these are the normal values of the cuts
  
if(!is.matrix(tau)) tau <- matrix(tau,ncol=nvar)
rownames(tau) <- paste0(xmin:(xmax-1))
colnames(tau) <- colnames(x)
if(dim(tau)[1] < dim(tau)[2]) tau <- t(tau) #rows are variables, columns are subjects
}
class(tau) <- c("psych","tau")  #added the tau class so score.irt can use the tau values
return(tau)
}


#added August 6, 2012
"irt.responses" <-
function(theta,items, breaks = 11,show.missing=FALSE,show.legend=TRUE,legend.location="topleft",colors=NULL,...) {
#if(is.null(colors)) colors =c("gray0", "blue3", "red3", "darkgreen", "gold2", "gray50", "cornflowerblue", "mediumorchid2")
if(is.null(colors)) colors =c("black", "blue", "red", "darkgreen", "gold2", "gray50", "cornflowerblue", "mediumorchid2")
#legend.location <- c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center","none")
#uniqueitems <- unique(as.vector(unlist(items)))
item.counts <- names(table(as.vector(unlist(items))))
uniqueitems <- as.numeric(item.counts)
nalt <- length(uniqueitems) + 1 #include the missing value from response.frequencies
nvar <- ncol(items) 
theta.min <- min(theta,na.rm=TRUE)
theta.max <- max(theta,na.rm=TRUE)
binrange <- cut(theta, breaks = breaks)
binnums <- as.numeric(binrange)
items <- as.matrix(items)
stats <- by(items,binnums,function(x) response.frequencies(x,uniqueitems=uniqueitems))
stats.m <- unlist(stats)
stats.m <- matrix(stats.m,ncol=nvar*nalt,byrow=TRUE)
theta <- seq(theta.min,theta.max,length.out=breaks)
for (i in 1:nvar) {
plot(theta,stats.m[,i],ylim=c(0,1),typ="l",xlab="theta",ylab="P(response)",main=paste(colnames(items)[i]),col=colors[1],...)
	for(j in 1:(nalt-2+show.missing)) {
			points(theta,stats.m[,i+nvar*j],typ="l",lty=(j+1),col=colors[j+1 ],...)
		}
	 if(show.legend) { legend(legend.location, paste(item.counts[1:(nalt-1 + show.missing)]), text.col = colors[1:(nalt-1+show.missing)], lty = 1:(nalt-1+show.missing), ncol=4,bty="n")}
	}}





#developed based on suggestions and code by David Condon
#scores multiple scales with full 2pl parameters
#gets around the problem of tau differences for 0/1 and 1/6 scales.
#Requires finding the correlation matrix for each scale, rather than taking advantage of a prior correlation matrix
scoreIrt.2pl <- function(itemLists,items,correct=.5,messages=FALSE,cut=.3,bounds=c(-4,4),mod="logistic") {
   nvar <- length(itemLists)
   select <- sub("-","",unlist(itemLists)) #select just the items that will be scored
    select <- select[!duplicated(select)]
   items <- items[select]  #this should reduce memory load
   #we turn off the sorting option in irt.fa so that the item discriminations match the scoring order
   #small function is called using parallel processing
   smallFunction <- function(i,selection,correct,cut=cut,bounds=bounds,mod=mod) {
      select <- sub("-","",selection[[i]])
      selectedItems <- as.matrix(items[select])
      if(!messages) {suppressMessages(stats <- irt.fa(selectedItems,correct=correct,plot=FALSE,sort=FALSE))} else {
                              stats <- irt.fa(selectedItems,correct=correct,plot=FALSE,sort=FALSE)}
      scores <- scoreIrt(stats,selectedItems,cut=cut,bounds=bounds,mod=mod)
      scores <- scores$theta
  }
   #use mapply for debugging, mcmapply for parallel processing
   #items is global and not passed to save memory
   scoresList <- mcmapply(smallFunction,c(1:nvar),MoreArgs=list(selection=itemLists,correct=correct,cut=cut,bounds=bounds,mod=mod))

   colnames(scoresList) <- names(itemLists)
   return(scoresList)
   }
   

 #A perhaps more robust way of calling scoreIrt is to find tau just for a few items at a time and then scoring one scale at a time.
 #the alternative is use scoreIrt for all of them at once with a keys function.
   
scoreIrt.1pl <- function(keys.list,items,correct=.5,messages=FALSE,cut=.3,bounds=c(-4,4),mod="logistic") {
   select <- sub("-","",unlist(keys.list))
      select <- select[!duplicated(select)]
   items <- items[select] 
  nf <- length(keys.list)
  fix <-  is.numeric(keys.list[[1]])
  
   smallFunction <- function(i,keys.list,correct,cut=cut,bounds=bounds,mod=mod) {
          list.i <- keys.list[[i]]
           keys <- rep(1,length(list.i))
            neg <- grep("-", list.i)
            keys[neg] <- -1
            select <- sub("-", "", list.i)
           # select <- colnames(items)[select]
      selectedItems <- as.matrix(items[select])
      stats <- irt.tau(selectedItems)
      scores <- scoreIrt(stats,selectedItems,keys,cut=cut,bounds=bounds,mod=mod)
     # stats <- irt.tau(items[select])
     # scores <- scoreIrt(stats,items[select],keys,cut=cut,bounds=bounds,mod=mod)

      scores <- scores[,1]
  }
   #use mapply for debugging, mcmapply for parallel processing
   #items are global and not passed
   scoresList <- mcmapply(smallFunction,c(1:nf),MoreArgs=list(keys.list=keys.list,correct=correct,cut=cut,bounds=bounds,mod=mod))
   
   colnames(scoresList) <- names(keys.list)
   return(scoresList)
   }
   

   

      
   

 
 
 
 
 
 
 
 
 