#Developed January 2018 while writing an article on reliability
#modified March 2021 in response to some helpful suggestions by Jacob Lee

"testRetest" <- function(t1,t2=NULL,keys=NULL,id="id",time= "time",select=NULL,check.keys=TRUE,warnings=TRUE,lmer=TRUE,sort=TRUE) {
 cl <- match.call() 
 
 #first, some basic checks of the data
 #there are several ways we can input the data
 #x and y as time 1 and time 2
 #x includes a time variable and an id variable
 x <- t1
 
 y <- t2
 if(NCOL(x) ==1) {just.test <-TRUE } else {just.test <- FALSE}
 keys.orig <- keys
  #first check if we have a y variable, if not, then create it
 
 if(is.null(y)) {n.times <- table(x[time])    
  #sort the data by id and time
  if(id %in% colnames(x)) {
 if(sort) x <- dfOrder(x,c(time,id),) }   #don't sort if we don't have ids
 y <- x[x[,time] == names(n.times)[2],]
 x <- x[x[,time] == names(n.times)[1],]
 }  else {
 if(sort) {x <- dfOrder(x,c(time,id),) 
           y <- dfOrder(y,c(time,id),)
          }} 
 
  #only take matching cases   
if(NROW(x) != NROW(y) ) {warning("The number of subjects in x do not  match those in y")
  not.missing.x <- x[,id] %in% y[,id]
  not.missing.y <- y[,id] %in% x[,id]
  missing.idx <- x[!not.missing.x,id]
  missing.idy <- y[!not.missing.y,id]
  cat("\nThe non-matched subjects were\n",missing.idx, missing.idy,"\n I have deleted them")
  
  x <- x[not.missing.x,]
  y <- y[not.missing.y,]
  } 
 
 n.obs <- NROW(x)
 if(!just.test) {
 if(!is.null(select)) {items <- select
    x <- x[select]
    y <- y[select]
    }

 if(is.null(keys)){ items <- colnames(x) [!colnames(x) %in% c(id,time)]} else {items <- keys }
  #first check if we should reverse any items and convert location numbers (if specified) to location names
 n.items <- length(items)
  if(is.character(items)) {
   temp <- rep(1,n.items)
 
   temp [strtrim(items,1)=="-"] <- -1
   if(any( temp < 0) )  {items <- sub("-","",items) }
   } else {temp <- sign(items)
      items <- colnames(x)[abs(items)] 
    }
 #check for bad input   -- the Mollycoddle option 
if(any( !(items %in% colnames(x)) )) {
 cat("\nVariable names in keys are incorrectly specified. Offending items are ", items[which(!(items %in% colnames(x)))],"\n")
 stop("I am stopping because of improper input in the scoring keys.  See the list above for the bad item(s). ")}  
 
 x <- x[,items,drop=FALSE]
 y <- y[,items,drop=FALSE]
 
#these are the means of the unreversed items
if(NCOL(x) > 1) {mean.x <- colMeans(x,na.rm=TRUE)
mean.y <- colMeans(y,na.rm=TRUE)
}

 


     min.item <- min(x[items],na.rm=TRUE)
     max.item <- max(x[items],na.rm=TRUE)
     miny.item <- min(y[items],na.rm=TRUE)
     maxy.item <- max(y[items],na.rm=TRUE)
  if(any(temp < 0)) {
   #flip items that are negatively keyed
  	 x[items[temp <0]] <- max.item- x[items[temp < 0]] + min.item 
  	 y[items[temp <0]] <- maxy.item- y[items[temp < 0]] + miny.item 
    }
#x and y are now scored in the direction of the keys    
   select <- items
   
  if(any(colnames(x[select]) !=colnames(y[select]))) {stop("Variable names must match across tests")}
  
  
      p1 <- pca(x)
      p2 <- pca(y)
      
#Even though, if given keys, we have already flipped them, we check one more time
    #  if(is.null(keys) ) { keys <- rep(1,n.items) }  else {keys <- temp }
     keys <- rep(1,n.items)
      if((any(p1$loadings < 0)) | (any(p2$loadings < 0))) {if (check.keys) {if(warnings) message("Some items were negatively correlated with total scale and were automatically reversed.\n This is indicated by a negative sign for the variable name.") 
                    keys[p1$loadings < 0] <- -1 } else {
                       if(is.null(keys) && warnings ) {message("Some items were negatively correlated with the total scale and probably \nshould be reversed.  \nTo do this, run the function again with the 'check.keys=TRUE' option")
                       if(warnings) cat("Some items (",rownames(p1$loadings)[(p1$loadings < 0)],") were negatively correlated with the total scale and \nprobably should be reversed.  \nTo do this, run the function again with the 'check.keys=TRUE' option")
                        }} 
    }
 
if(any(keys < 0)) {   
#then find the x and y scores
newx <- t(t(x[select]) * keys + (keys < 0)*  (max.item + min.item) )  #these are now rescaled in the keyed direction  - but we have already done this if given a keys vector
newy <- t(t(y[select]) * keys + (keys < 0)* (maxy.item + miny.item))
} else {
newx <-  x[select]
newy <- y[select]
}

xscore <- rowMeans(newx,na.rm=TRUE)
yscore <- rowMeans(newy,na.rm=TRUE)
r12 <- cor(xscore,yscore,use="pairwise")
#then correlate them to get test retest r

#Now find the alpha for the x and y scales
x.alpha <- alpha.1(cov(newx,use="pairwise"))
y.alpha <- alpha.1(cov(newy,use="pairwise"))
xy.alpha <- rbind(unlist(x.alpha),unlist(y.alpha))
rownames(xy.alpha) <- c("x","y")
colnames(xy.alpha) <- c("raw G3","std G3","G6","av.r","S/N","se","lower","upper","var.r")
#then correlate each matched item

          
dxy <- dist(newx,newy) 
rqq <- dxy$rqq

#now find the item over subjects correlation
rii <- rep(NA,n.items)
for (j in (1:n.items)) {
 if(!(( is.na(sd(x[,items[j]],na.rm=TRUE))) | (is.na(sd(y[,items[j]],na.rm=TRUE)))))  {
  rii[j] <- cor(x[,items[j]],y[,items[j]],use="pairwise")}
 }  


 #ok, the data seem ok lets create a dummy variable and do the lmer on it
 n.obs <- min(NROW(newx),NROW(newy))
 xy.df <- data.frame(id = rep(1:n.obs,2), time=c(rep(1,n.obs), rep(2,n.obs)), rbind(newx,newy),row.names=1:(2*n.obs))   #this is getting it ready for mlr1 
 } else {#The case of  just two tests, no items 
 

   xy.df <- data.frame(id=rep(1:n.obs,2),time=c(time=c(rep(1,n.obs), rep(2,n.obs))),rbind(x,y)[,1],row.names=1:(2*n.obs))
  
  no.items <- TRUE
 } 
 
  
#Now, the data are ready for lmer  
#This is same as 
# ml <- mlr(xy.df,aov=FALSE,lmer=TRUE)   We now need to grab the best of multilevel reliabiity to do the next part



if(!just.test) {
  if(lmer) {ml <- mlr1(xy.df)} else {ml <- list(n.obs=n.obs,n.items=n.items)}
    if(is.null(keys.orig)) keys.orig <- rep(1,n.items)
    item.stats <- data.frame(rii=rii,p1=unclass(p1$loadings),p2=unclass(p2$loadings),mean1 = mean.x, mean2=mean.y, keys=keys,keys.orig=keys.orig)
    colnames(item.stats)[2:3] <- c("PC1", "PC2")
    key <- rownames(item.stats)
    key[item.stats$keys < 0] <- paste0("-", key[item.stats$keys < 0])
    scores <- data.frame(pca1 = p1$scores,pca2 = p2$scores,t1scores =xscore, t2scores = yscore,rqq=rqq,dxy=dxy$dxy,t1sd=dxy$sdx,t2sd=dxy$sdy)
     result <- list(r12=r12,alpha=xy.alpha,rqq=rqq,dxy=dxy,item.stats=item.stats, scores=scores,xy.df=xy.df,key=key,ml=ml,Call=cl) } else {
    if(just.test) {r12 = cor(x,y,use="pairwise")
             ml <- mlr2(xy.df)
                 result <- list(r12 =r12,ml=ml, Call=cl)} }
     class(result) <- c("psych", "testRetest")
     return(result)  
} 


########

 alpha.1 <- function(C,R=NULL) {
    n <- dim(C)[2]
    alpha.raw <- (1- tr(C)/sum(C))*(n/(n-1))
    if(is.null(R)) R <- cov2cor(C)
    sumR <- sum(R)
    alpha.std <-  (1- n/sumR)*(n/(n-1))  
    smc.R <- smc(R)
    G6 <- (1- (n-sum(smc.R))/sumR)
    av.r <- (sumR-n)/(n*(n-1))
    R.adj <- R
    diag(R.adj) <- NA
    var.r  <- var(as.vector(R.adj),na.rm=TRUE)
    mod1 <- matrix(av.r,n,n)
    Res1 <- R - mod1
    GF1 =  1- sum(Res1^2)/sum(R^2)
    Rd <- R - diag(R)
    diag(Res1) <- 0
    GF1.off <- 1 - sum(Res1^2)/sum(Rd^2)  
    sn <- n*av.r/(1-av.r)
   # Q = (2 * n^2/((n-1)^2*(sum(C)^3))) * (sum(C) * (tr(C^2) + (tr(C))^2) - 2*(tr(C) * sum(C^2))) #corrected 1/15/16 
    Q = (2 * n^2/((n - 1)^2 * (sum(C)^3))) * (sum(C) * (tr(C%*%C) +  (tr(C))^2) - 2 * (tr(C) * sum(C%*%C)))   #correction from Tamaki Hattori
    result <- list(raw=alpha.raw,std=alpha.std,G6=G6,av.r=av.r,sn=sn,Q=Q,GF1,GF1.off,var.r = var.r)
    return(result)
    }
    
    

    
print_psych.testRetest<- function(x,digits=2,short=FALSE,...) {
cat("\nTest Retest reliability ")
cat("\nCall: ")
	print(x$Call)
cat('\nNumber of subjects = ',x$ml$n.obs, " Number of items = ", x$ml$n.items)
if(x$ml$n.items  > 1)  { #The normal case
cat("\n Correlation of scale scores over time" , round(x$r12,digits))
cat("\n Alpha reliability statistics for time 1 and time 2 \n")
rownames(x$alpha) <- c("Time 1", "Time 2")
print( round(x$alpha,digits))

meanrii <- mean(x$item.stats$rii,na.rm=TRUE)
meanrqq <- mean(x$rqq,na.rm = TRUE)
sdrqq <- sd(x$rqq,na.rm = TRUE)
meandqq <- mean(x$dxy$dxy,na.rm=TRUE)

cat("\n Mean between person, across item reliability = ",round(meanrii,digits))

cat("\n Mean within person, across item reliability = ",round(meanrqq,digits))
cat(  "\nwith standard deviation of " ,round(sdrqq,digits) ,"\n")
cat("\n Mean within person, across item d2 = ",round(meandqq,digits))
temp <- x
x <- x$ml

if(!is.null(x$R1F)) cat("\nR1F  = ",round(x$R1F,digits) , "Reliability of average of all items for one  time (Random time effects)")

if(!is.null(x$RkF)) cat("\nRkF  = ",round(x$RkF,digits) , "Reliability of average of all items and both times (Fixed time effects)")
	
if(!is.null(x$R1R)) cat("\nR1R  = ",round(x$R1R,digits),"Generalizability of a single time point across all items (Random time effects)")
		
	
	if(!is.null(x$R2R)) cat("\nRkR  = ",round(x$RkR,digits),"Generalizability of average time points across all items (Fixed time effects)")
	
	if(!is.null(x$Rc))  cat("\nRc   = ",round(x$Rc,digits),"Generalizability of change (fixed time points, fixed items) ")
	  
	if(!is.null(x$RkRn) ) cat("\nRkRn = ",round(x$RkRn,digits),"Generalizability of between person differences averaged over time (time nested within people)")
	   
	if(!is.null(x$Rcn)) cat("\nRcn  = ",round(x$Rcn,digits),"Generalizability of within person variations averaged over items  (time nested within people)")
x <- temp	  
if(!is.null(x$ml$components)) {cat("\nMultilevel components of variance\n")
print(round(x$ml$components,digits))}

if(!short) {
cat("\n With Item statistics \n")
print(round(x$item.stats[-7],digits))
}  else {cat("\n To see the item.stats, print with short=FALSE. \nTo see the subject reliabilities and differences, examine the 'scores' object.") } 

} else { cat("\nTest Retest Reliability of two tests", print(round(x$r12,digits)))
  cat("\nMultilevel components of variance\n")
  print(round(x$ml$components,digits))
   }
   }
   
   
#######   
   #grab the best parts of multilevel reliability
   
mlr1 <- function(x,na.action= "na.omit" ) { 
   long <- NULL
   id <- "id"
   time <- "time"
   n.obs <- NROW(x)/2
   items <- colnames(x) [!colnames(x) %in% c("id","time")] 
   n.items <- length(items)
   n.time <- 2
 long <- data.frame(id = rep(1:n.obs,2), time=rep(1:2,each = n.obs),stack(x[items]))
 colnames(long)[4] <- "items"    #just to make it clearer
  mod.lmer <- lme4::lmer(values ~ 1 + (1 | id) + (1 | time) + (1 | items) + (1 | id:time)+ (1 | id:items)+ (1 | items :time),
  data=long,na.action=na.action)
   vc <- lme4::VarCorr(mod.lmer)
  MS_id <- vc$id[1,1]
  MS_time <- vc$time[1,1]
  MS_items <- vc$items[1,1]
 # MS_pxt <- vc[[1]][[1]]
 # MS_pxitem <- vc[[2]][[1]]
 # MS_txitem <- vc[[3]][[1]]              #changed to named locations (which were incorrect when doing by numbers)
  MS_pxt <- vc[["id:time"]][[1]]
  MS_pxitem <- vc[["id:items"]][[1]]
  MS_txitem <- vc[["items:time"]][[1]]

  error <- MS_resid <- (attributes(vc)$sc)^2
  s.lmer <- s.aov <- summary(mod.lmer)
 MS.df <- data.frame(variance= c(MS_id, MS_time ,MS_items, MS_pxt, MS_pxitem, MS_txitem, MS_resid,NA))
  rownames(MS.df) <- c("ID","Time","Items","ID x time", "ID x items", "time x items", "Residual","Total") 

 MS.df["Total",]  <- sum(MS.df[1:7,1],na.rm=TRUE)
MS.df["Percent"] <- MS.df/MS.df["Total",1]
lmer.MS <- MS.df  #save these
#fixed time, not random time,   
  R1f <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + error/( n.items)))
  #average of both times
  
  Rkf <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + error/(n.time * n.items)))
R1r <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + MS_time + MS_pxt + error/( n.items)))  #per Sean Lane
Rkr <- (MS_id + MS_pxitem/n.items)/((MS_id + MS_pxitem/n.items + MS_time/n.time + MS_pxt/n.time + error/( n.time * n.items)))
Rc <- (MS_pxt)/(MS_pxt + error/n.items)

result <- list(n.obs = n.obs, n.items=n.items, components = MS.df,R1F= R1f,RkF =Rkf,R1R = R1r,RkR = Rkr,Rc=Rc)
 return(result)
  }
  
  mlr2 <- function(x,na.action=na.omit) {
  
  #these treats the case of just two tests with no items, we want the variance components
  long <- x   #it is actually already in long format
  n.obs <- NROW(long) /2
  n.items <- 1
   mod.lmer <- lme4::lmer(values ~ 1 + (1 | id) + (1 | time) ,
  data=long,na.action=na.action)
    
    
 vc <- lme4::VarCorr(mod.lmer)
MS_id <- vc$id[1,1]


 error <- MS_resid <- (attributes(vc)$sc)^2
    MS.df <- data.frame(variance= c(MS_id, MS_resid,NA))
  rownames(MS.df) <- c("ID","Residual","Total") 

MS.df["Total",]  <- sum(MS.df[1:2,1],na.rm=TRUE)
 MS.df["Percent"] <- MS.df/MS.df["Total",1]   
 Rxx <- MS_id/MS.df["Total",1]
  result <- list(n.obs=n.obs,n.items=n.items,components = MS.df)    
  }

#find level, scatter, and pattern by subject  
dist <- function(x,y) {
x.level <- rowMeans(x,na.rm=TRUE)
y.level <- rowMeans(y,na.rm=TRUE)
n.obs <- NROW(x)
sdxi <- apply(x,1,function(xx) sd(xx,na.rm=TRUE))
sdyi <- apply(y,1,function(xx) sd(xx,na.rm=TRUE))
dxy <- rowMeans((x - y)^2,na.rm=TRUE)
rxy <- rep(NA,n.obs)
tx <- t(x)
ty <- t(y)
for(i in 1:n.obs) {
 if(!(    (is.na(sdxi[i])) | (sdxi[i]==0) | (is.na(sdyi[i]) | sdyi[i]==0) )     ) {
rxy[i]  <- cor(tx[,i],ty[,i],use="pairwise")}
}
dist.df <- data.frame(x.level=x.level,y.level=y.level,sdx=sdxi,sdy = sdyi,dxy=dxy,rqq=rxy)
return(dist.df)
}