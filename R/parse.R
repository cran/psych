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
     	
     	