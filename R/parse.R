 #This parses a formula like input and return the left hand side variables (y) and right hand side (x) as well as products (prod)  and partials (-)
 #  		 
 fparse <- function(expr){
      	 m <- prod <- ex <- ex2 <-  NULL
		 all.v <- all.vars(expr) 
		 te <- terms(expr)   #this will expand the expr for products
		 fac <- attributes(te)$factors
		 x <- rownames(fac)[-1] #drop the y variables
		# y <- all.v[!all.v %in% x]   
		z <- rownames(fac)[rowSums(fac) < 1]   #what does this do?
		 if(length(z) > 1)  {z <- z[-1]
		      x <- x [! x%in%z]} else {z <- NULL}
		      
		 char.exp <- as.character(expr[3])
		 #strip out exponential terms from x
		 notx <-  regmatches(char.exp, gregexpr("I\\(.*?\\)", char.exp))[[1]]
		 x <- x[!x %in%notx]
		 ex1 <- gsub("I[\\(\\)]", "", regmatches(char.exp, gregexpr("I\\(.*?\\)", char.exp))[[1]])  #look for I(x)
		 if (length(ex1)  >0) {ex <- sub("\\)","",ex1)
		 }
		 x <- x[ ! x %in% ex]
		
		 #now look for mediators
		 m <- gsub("[\\(\\)]", "", regmatches(char.exp, gregexpr("\\(.*?\\)", char.exp))[[1]])
		 if(length(m)<1) {m <- NULL} else {m <- m[! m %in% ex] }
         if(length(m) < 1) m <- NULL
         prod.terms <- sum(attributes(te)$order > 1)
         if(prod.terms > 0 ) {
          n1 <- sum(attributes(te)$order == 1)
          prod <- list()
          for(i in(1:prod.terms)) {
          prod[[i]] <- names(which(fac[,n1+i] > 0)) } 
         }
         #now, if there are ex values, get rid of the ^2
         if(!is.null(ex)) {ex <- sub("\\^2","",ex)
         }
         y <- all.v[ ! all.v %in% c(x,z,ex) ]
      return(list(y=y,x=x,m=m,prod=prod,z = z,ex=ex))
      }
      
#convert lavaan cfa instructions into matrix form
#12/11/25 to help CFA  
 #`1/15/26 added the ability to specify correlations for sims`   
lavParse <- function(model,phi=FALSE){
short <-gsub("\t","",model) #drop all tabs
short <- gsub(" ", "", short )  #drop all blanks
#first, break by lines

lines <- strsplit(short,"\n") #break into a new line for each factor
fact <- strsplit(lines[[1]],"=~") #  break into factor names and variables

#three steps
#find the factor and variable names
fact <- fact[lapply(fact,length)>0] 
covs <- fact[lapply(fact,length) ==1] 
ncovs <- length(covs)
nfact <- length(fact) - ncovs
 #drop empty lines due to bad typing
fnames <- 1:nfact
vnames<- NULL
 for (i in 1:nfact) {fnames[i]<- fact[[i]] [1]
    vect <- strsplit(fact[[i]][2],"\\+")  
    for (j in 1:length(vect[[1]])) {
    temp <- strsplit(vect[[1]][j],"\\*")
    if(length(temp[[1]]) > 1 ) {vnames <- c(vnames,temp[[1]][2])}  else {vnames <- c(vnames,vect[[1]][j])}
    }
    }
#do it for correlations as well as factor loadings
#add fnames to  the list
#for (i in 1:ncovs)  {
#cov <- strsplit(covs[[i]],"~")
#vnames <- c(vnames,cov[[1]][1])}
if(phi) vnames <- c(vnames, fnames)
    
 #we have the names, now fill the v.mat
vnames <- vnames[!duplicated(vnames)]
v.mat <- matrix(0,ncol=nfact,nrow=length(vnames))
colnames(v.mat) <- fnames
rownames(v.mat) <- vnames
#now put the 1s in
for (i in 1:nfact) {fnames[i]<- fact[[i]] [1]
    vect <- strsplit(fact[[i]][2],"\\+")
    
    for (j in 1:length(vect[[1]])) {
     temp <- strsplit(vect[[1]][j],"\\*")    #allow us to specify loading values
     if(length(temp[[1]]) > 1 ) {v.mat[temp[[1]][2],i] <- temp[[1]][1]}  else {v.mat[vect[[1]][j],i]  <- 1} }

    }
    
#by default, create the identity mattrix
if(phi){

for(i in 1:nfact) {
v.mat[fnames[i],fnames[i]] <- 1}
   if(ncovs > 0) {
 for(i in 1:ncovs){
    vect <- strsplit(covs[[i]],"~")
    fac <- vect[[1]][1]    #this is the rowvalue for the factor 
    temp <- strsplit(vect[[1]][2],"\\+")
    ntemp <- length(temp[[1]])
    coeff <- strsplit(temp[[1]],"\\*")
    for(k in 1:ntemp) { v.mat[fac,coeff[[k]][2]  ] <- coeff[[k]][1]
                        v.mat[coeff[[k]][2],fac  ] <- coeff[[k]][1]}
     }                   
    }
    }
 #
  v.mat <-  as.matrix(nchar2numeric(as.data.frame(v.mat),flag=FALSE))
return(v.mat)}



model <- 'F1 =~ .9*V1 + .8*V2 + .7*V3 
         F2 =~ .8 * V4 + .7*V5 + .6*V6
        F3 =~ .9*V7 + .7*V8 + .5*V9
        F1 ~ .3 * F2 + .6* F3
        F1 ~ .2 * F3
        F2 ~ .5*F3'
     		