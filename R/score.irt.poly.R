#revised August 31, 2012 to allow for estimation of all 0s or all 1s
"score.irt.2" <- 
function(stats,items,keys=NULL,cut=.3,bounds=c(-5,5),mod="logistic") {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
#This does the normal fit
irt.2par.norm <-  function(x,delta,beta,scores) {
  fit <- -1*(log(scores*(1-pnorm(beta*(delta-x))) + (1-scores)*(1-pnorm(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  } 
#This does the logistic fit 
 irt.2par <- function(x,delta,beta,scores) {
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  
 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]
 nvar <- dim(items)[2]  
scores <- matrix(NaN,nrow=n.obs,ncol=nf*3)
 #now do the scoring one factor at a time
 
for (f in 1:nf) {
diff <- stats$difficulty[[f]]
if(nf < 2) {discrim <- drop(stats$discrimination)
            if(!is.null(keys)) {discrim <- discrim * abs(keys)
                                }
} else {discrim <- stats$discrimination[,f] 
if(!is.null(keys)) {discrim <- discrim * abs(keys[,f])
                    }}

fit <- rep(NA,n.obs)
theta <- rep(NA,n.obs)
items.f <-t(items)
items.f[abs(discrim) < cut] <- NA
items.f <- t(items.f)

if(is.matrix(discrim)) discrim <- drop(discrim)
 total <- rowMeans(t(t(items.f)*sign(discrim)),na.rm=TRUE)
 count <- rowSums(!is.na(items.f))

 for (i in 1:n.obs) {
 	progressBar(i,n.obs,"score.IRT")
 	if (count[i]>0) {if((sum(items.f[i,],na.rm=TRUE) ==0 )  | (prod(items.f[i,],na.rm=TRUE) ==  1 )) { 
 	if(sum(items.f[i,],na.rm=TRUE) ==0 ) {
 		if(mod=="logistic") {p <- log(1-(1-items.f[i,])/(1+exp(discrim*(diff))) )} else {  #logistic
 				p <- log(1-(pnorm(items.f[i,]*discrim*diff))) }  #normals
  	 pall <- exp(sum(p,na.rm=TRUE))
   
   	theta[i] <- qnorm(pnorm(qnorm(pall))/2)  #the z value of 1/2 the quantile value of pall
   	fit[i] <- 0 
 	 # cat ("\nThe case of all wrong",i,theta[i])

 	} else {
 		if(mod == "logistic") {	p <- log((items.f[i,])/(1+exp(discrim*(diff))) )} else  { 
 	                        	p <- log((items.f[i,])*(1 - pnorm(1- discrim*(diff)) ))  } 
 	
 	 pall <- exp(sum(p,na.rm=TRUE))
   theta[i] <- qnorm(1-pnorm(qnorm(pall))/2)  #the z value of 1/2 the quantile value of pall
   fit[i] <- 0 


 	#cat("\nthe case of all right",i,theta[i]) 
 }
 	} else { #cat("the normal case",i )
 	if(mod=="logistic") {
  myfit <-optimize(irt.2par,bounds,beta=discrim,delta=diff,scores=items.f[i,])  #how to do an apply?
    } else {myfit <-optimize(irt.2par.norm,bounds,beta=discrim,delta=diff,scores=items.f[i,])} #do a normal fit function
     		theta[i] <- myfit$minimum    
     		fit[i] <- myfit$objective  #fit of optimizing program 
     		}} else  {#cat("\nno items",i)

     		 theta[i] <- NA
  		      fit[i]  <- NA 
 	                    
 	        }  #end if count ... else          
 	                     	          
    }  #end loop 
    theta [theta < bounds[1]] <- bounds[1]
    theta[theta > bounds[2]] <- bounds[2]
     if((!is.null(keys)) & (all(keys[,f] == -1) || (sign(cor(discrim,keys[,f])) < 0) )) {theta <- -theta
                                 total <- -total}
	scores[,f] <- theta
	scores[,nf + f] <- total
	scores[,2*nf + f] <- fit
	} #end of f loop
	
	colnames(scores) <- paste(rep(c("theta","total","fit"),each=nf),1:nf,sep="")
	return(scores)
}
	

"score.irt.poly" <- 
function(stats,items,keys=NULL,cut=.3,bounds=c(-5,5)) {
#find the person parameters in a 2 parameter model we use deltas and betas from irt.discrim and irt.person.rasch
#find the person  parameter
#created July 4, 2011
 irt.2par.poly <- function(x,delta,beta,scores) {
  fit <- -1*(log(scores/(1+exp(beta*(delta-x))) + (1-scores)/(1+exp(beta*(x-delta)))))
  mean(fit,na.rm=TRUE)
  }  

 min.item <- min(items,na.rm=TRUE)
 items <- items - min.item #this converts scores to positive values
 max.item <- max(items,na.rm=TRUE)  #we use this when reverse score
 nf <- length(stats$difficulty)
 n.obs <- dim(items)[1]
 nvar <- dim(items)[2]
 
scores <- matrix(NaN,nrow=n.obs,ncol=nf*3)

for (f in 1:nf) {
diff <- stats$difficulty[[f]]
if(nf < 2) {discrim <- stats$discrimination
              if(!is.null(keys)) {discrim <- discrim * abs(keys)
                              } } else {discrim <- stats$discrimination[,f] 
                                if(!is.null(keys)) {discrim <- discrim * abs(keys[,f]) 
                  }  
                    }
cat <- dim(diff)[2]


diff.vect <- as.vector(t(diff)) 
discrim.vect <- rep(discrim,each=cat)


 fit <- rep(NA,n.obs)
 theta <- rep(NA,n.obs)
 item.f <- t(items)
 item.f[abs(discrim) < cut] <- NA  #this does not change the item, just the temp version of the item
 item.f <- t(item.f)

 total <- rowMeans(t(t(item.f )* as.vector(sign(discrim))),na.rm=TRUE) #fixed 11/11/11 to be as.vector
 num.reversed <- sum(discrim < -cut)  #the number of negatively keyed items
 num.keyed <- sum(abs(discrim) > cut)

 total <- total + num.reversed * (max.item)/num.keyed + min.item
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
   
     if((!is.null(keys)) & (all(keys[,f] == -1) || (sign(cor(discrim,keys[,f])) < 0) )) {theta <- -theta
                                 total <- -total}
	scores[,f] <- theta
	scores[,nf + f] <- total
	scores[,2*nf + f] <- fit
	}
	colnames(scores) <- paste(rep(c("theta","total","fit"),each=nf),1:nf,sep="")
	return(scores)
	}

"score.irt" <- function(stats=NULL,items,keys=NULL,cut=.3,bounds=c(-5,5),mod="logistic") {
if (length(class(stats)) > 1) {
   #if(!is.null(keys) && is.vector(keys)) keys <- as.matrix(keys)
    switch(class(stats)[2],
    irt.poly = {scores <- score.irt.poly(stats$irt,items,keys,cut,bounds=bounds) },
    irt.fa =   {scores <- score.irt.2(stats$irt,items,keys,cut,bounds=bounds,mod=mod)},
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


 