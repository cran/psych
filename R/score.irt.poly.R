"score.irt.2" <- 
function(stats,items,keys=NULL,cut=.3,bounds=c(-3.5,3.5)) {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
 irt.2par <- function(x,delta,beta,scores) {
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  
 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]
 nvar <- dim(items)[2]  
scores <- matrix(NA,nrow=n.obs,ncol=nf*3)
 #now do the scoring one factor at a time
 
for (f in 1:nf) {
diff <- stats$difficulty[[f]]
if(nf < 2) {discrim <- drop(stats$discrimination)
            if(!is.null(keys)) discrim <- discrim * abs(keys)
} else {discrim <- stats$discrimination[,f] 
if(!is.null(keys)) discrim <- discrim * abs(keys[,f])
}

fit <- rep(NA,n.obs)
theta <- rep(NA,n.obs)
items.f <-t(items)
items.f[abs(discrim) < cut] <- NA
items.f <- t(items.f)

 total <- rowMeans(t(t(items.f)*sign(discrim)),na.rm=TRUE)
 count <- rowSums(!is.na(items.f))
 
 for (i in 1:n.obs) {
 	
 	if (count[i]>0)  {myfit <- optimize(irt.2par,bounds,beta=discrim,delta=diff,scores=items.f[i,]) #how to do an apply?
     		theta[i] <- myfit$minimum    
     		fit[i] <- myfit$objective  #fit of optimizing program
     		} else {
            theta[i] <- NA
  		      fit[i]  <- NA 
    			}    #end if else
    }  #end loop 
	scores[,f] <- theta
	scores[,nf + f] <- total
	scores[,2*nf + f] <- fit
	}#end of f loop
	colnames(scores) <- paste(rep(c("theta","total","fit"),each=nf),1:nf,sep="")
	return(scores)
}
	

"score.irt.poly" <- 
function(stats,items,keys=NULL,cut=.3,bounds=c(-3.5,3.5)) {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
#created July 4, 2011
 irt.2par.poly <- function(x,delta,beta,scores) {
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  

 min.item <- min(items,na.rm=TRUE)
 items <- items - min.item #this converts scores to positive values
 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]
 nvar <- dim(items)[2]
 
scores <- matrix(NA,nrow=n.obs,ncol=nf*3)

for (f in 1:nf) {
diff <- stats$difficulty[[f]]
if(nf < 2) {discrim <- stats$discrimination} else {discrim <- stats$discrimination[,f] }
cat <- dim(diff)[2]
if(!is.null(keys)) discrim <- discrim * abs(keys[,f])
diff.vect <- as.vector(t(diff)) 
discrim.vect <- rep(discrim,each=cat)


 fit <- rep(NA,n.obs)
 theta <- rep(NA,n.obs)
 item.f <- t(items)
 item.f[abs(discrim) < cut] <- NA
 item.f <- t(item.f)
 total <- rowMeans(t(t(item.f )* sign(discrim)),na.rm=TRUE)
 count <- rowSums(!is.na(item.f))
 
 for (subj in 1:n.obs) {
 	
 	if (count[subj]>0)  {
 	       newscore <-  NULL
 	       score <- item.f[subj,]
 	     
           for (i in 1:nvar) {
                 if(is.na(score[i])) {newscore <- c(newscore,rep(NA,cat)) } else {
                    if(score[i] == cat) {newscore <- c(newscore,rep(1,score[i]))
                    } else {
           newscore <- c(newscore,rep(1,score[i]),rep(0,cat-score[i])) }}
                     }
 	        myfit <- optimize(irt.2par.poly,bounds,beta=discrim.vect,delta=diff.vect,scores=newscore) 
     		theta[subj] <- myfit$minimum    
     		fit[subj] <- myfit$objective  #fit of optimizing program
     		} else {
         	fit[subj] <- NA
  			theta[subj] <- NA 
    			}    #end if else
    }  #end loop 
	scores[,f] <- theta
	scores[,nf + f] <- total
	scores[,2*nf + f] <- fit
	}
	colnames(scores) <- paste(rep(c("theta","total","fit"),each=nf),1:nf,sep="")
	return(scores)
	}

"score.irt" <- function(stats=NULL,items,keys=NULL,cut=.3,bounds=c(-3.5,3.5)) {
if (length(class(stats)) > 1) {
    switch(class(stats)[2],
    irt.poly = {scores <- score.irt.poly(stats$irt,items,keys,cut,bounds=bounds) },
    irt.fa =   {scores <- score.irt.2(stats$irt,items,keys,cut,bounds=bounds)},
    fa = {tau <- irt.tau(items)  #this is the case of a factor analysis to be applied to irt
          nf <- dim(stats$loadings)[2]
          diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-stats$loadings[,i]^2) }
     
         discrim <- stats$loadings/sqrt(1-stats$loadings^2)
   class(diffi) <- NULL
   class(discrim) <- NULL
   new.stats <- list(difficulty=diffi,discrimination=discrim)
    scores <- score.irt.poly(new.stats,items,keys,cut,bounds=bounds)})
       #we should have a null case
    } else {#input is a keys matrix
         tau <- irt.tau(items)  #this is the case of a using a scoring matrix to be applied to irt
         if(is.matrix(keys)) {nf <- dim(keys)[2]} else {nf <-1}
         diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau }
         discrim <- keys
   class(diffi) <- NULL
   class(discrim) <- NULL
   new.stats <- list(difficulty=diffi,discrimination=discrim)
   if(dim(tau)[2] ==1)  {scores <- score.irt.2(stats=new.stats,items=items,cut=cut,bounds=bounds)} else {
                         scores <- score.irt.poly(stats=new.stats,items=items,cut=cut,bounds=bounds)}
                         }
scores <- data.frame(scores)
return(scores)
}



"irt.tau" <- function(x) {
x <-as.matrix(x)
nvar <- dim(x)[2]
xt <- table(x)
nvalues <- length(xt)  #find the number of response alternatives 
if(nvalues ==2) {tau <- -qnorm(colMeans(x,na.rm=TRUE))
   tau <- as.matrix(tau)
   rownames(tau) <- colnames(x)}  else {
if(nvalues > 10) stop("You have more than 10 categories for your items, polychoric is probably not needed")
xmin <- min(x,na.rm=TRUE)
xfreq <- apply(x- xmin+ 1,2,tabulate,nbins=nvalues)
n.obs <- colSums(xfreq)
xfreq <- t(t(xfreq)/n.obs)
tau <- qnorm(apply(xfreq,2,cumsum))[1:(nvalues-1),]  #these are the normal values of the cuts
  
if(!is.matrix(tau)) tau <- matrix(tau,ncol=nvar)
rownames(tau) <- names(xt)[1:(nvalues-1)]
colnames(tau) <- colnames(x)
if(dim(tau)[1] < dim(tau)[2]) tau <- t(tau) #rows are variables, columns are subjects
}
return(tau)
}