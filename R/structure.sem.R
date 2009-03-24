#just structure.graph, without the graphics (graphics commented out)
"structure.sem" <-
function(fx,Phi=NULL,fy=NULL, out.file=NULL,labels=NULL,cut=.3,errors=TRUE,simple=TRUE,regression=FALSE) {
    	
  xmodel <- fx
  ymodel <- fy
 if(!is.null(class(xmodel)) && (length(class(xmodel))>1)) {
   if(class(xmodel)[1] =="psych" && class(xmodel)[2] =="omega") {
    Phi <- xmodel$schmid$phi
    xmodel <- xmodel$schmid$oblique} else {
   if(class(xmodel)[1] =="psych" && ((class(xmodel)[2] =="fa") | (class(xmodel)[2] =="principal"))) { if(!is.null(xmodel$Phi)) Phi <- xmodel$Phi
        xmodel <- as.matrix(xmodel$loadings)} 
         }} else {
 if(!is.matrix(xmodel) & !is.data.frame(xmodel) &!is.vector(xmodel))  {
        if(!is.null(xmodel$Phi)) Phi <- xmodel$Phi
        xmodel <- as.matrix(xmodel$loadings)
      } else {xmodel <- xmodel} 
     }
     
  digits <- 2
      
 if(!is.matrix(xmodel) ) {factors <- as.matrix(xmodel)} else {factors <- xmodel}
 
 
 
  
  #first some basic setup parameters 
  num.y <- 0   #we assume there is nothing there
  num.var  <- num.xvar <- dim(factors)[1]   #how many x variables?
  if (is.null(num.var) ){num.var <- length(factors)
                          num.factors <- 1} else {
            num.factors <- dim(factors)[2]}
   
   num.xfactors <- num.factors
  
  if(is.null(labels)) {vars <- xvars <-  rownames(xmodel)} else { xvars <-  vars <- labels}
 
 if(is.null(vars) ) {vars <-xvars <- paste("x",1:num.var,sep="")  }
  fact <- colnames(xmodel)
  if (is.null(fact)) { fact <- paste("X",1:num.factors,sep="") }
  
  
 num.yfactors <- 0
   if (!is.null(ymodel)) { 
   if(is.list(ymodel) & !is.data.frame(ymodel)  ) {ymodel <- as.matrix(ymodel$loadings)} else {ymodel <- ymodel}   
 if(!is.matrix(ymodel) ) {y.factors <- as.matrix(ymodel)} else {y.factors <- ymodel}
   
   num.y <- dim(y.factors)[1] 
      if (is.null(num.y)) {
         num.y <- length(ymodel)
         num.yfactors <- 1} else {
         num.yfactors <- dim(y.factors)[2]
        }  
      
     yvars <- rownames(ymodel)
     if(is.null(yvars)) {yvars <- paste("y",1:num.y,sep="")  }
      if(is.null(labels)) {vars <- c(xvars,yvars)} else {yvars <- labels[(num.xvar+1):(num.xvar+num.y)]} 
    
     vars <- c(vars,yvars)
     
      yfact <- colnames(ymodel)
      if(is.null(yfact)) {yfact <- paste("Y",1:num.yfactors,sep="") }
      fact <- c(fact,yfact)
     num.var <- num.xvar + num.y
     num.factors <- num.xfactors + num.yfactors
     }
    sem <- matrix(rep(NA),6*(num.var*num.factors + num.factors),ncol=3)
    colnames(sem) <- c("Path","Parameter","Value")
   
   if(!regression) {    #the normal condition is to draw a latent model
           k <- num.factors  
      
  	if (num.xfactors ==1) {    
    	for (i in 1:num.xvar) { 
                          sem[i,1] <- paste(fact[1],"->",vars[i],sep="")
                          if(is.numeric(factors[i])) {sem[i,2] <- vars[i]} else {sem[i,2] <- factors[i] }
                                }  
                        k <- num.xvar+1 
     	} else {         #end of if num.xfactors ==1 
        #all loadings > cut in absolute value
                   k <- 1
                   for (i in 1:num.xvar) {
                   	for (f in 1:num.xfactors) { 
                   		if((!is.numeric(factors[i,f] ) && (factors[i,f] !="0"))||  ((is.numeric(factors[i,f]) && abs(factors[i,f]) > cut ))) {
               
                             sem[k,1] <- paste(fact[f],"->",vars[i],sep="")
                             if(is.numeric(factors[i,f])) {sem[k,2] <- paste("F",f,vars[i],sep="")} else {sem[k,2] <- factors[i,f]}
                              k <- k+1 }   #end of if 
                         }  
                        }      
       } 
        if(errors) {  for (i in 1:num.xvar) { 
                                          sem[k,1] <- paste(vars[i],"<->",vars[i],sep="")
                                          sem[k,2] <- paste("x",i,"e",sep="")
                                             k <- k+1 }
                    }
    } else  {   #the regression case
       if (title=="Structural model") title <- "Regression model"
       k <- num.var+1
       yvars <- "Y1"
              }   
       
  #now, if there is a ymodel, do it for y model 
  
  if(!is.null(ymodel)) { 
  if (num.yfactors ==1) {     
    for (i in 1:num.y) {
                          sem[k,1] <- paste(fact[1+num.xfactors],"->",yvars[i],sep="")
                           if(is.numeric(y.factors[i] ) ) {sem[k,2] <- paste("Fy",yvars[i],sep="")} else {sem[k,2] <- y.factors[i]}
                         k <- k +1
                        }                      
      } else {   #end of if num.yfactors ==1 
        #all loadings > cut in absolute value
                   for (i in 1:num.y) {
                   for (f in 1:num.yfactors) { #if (!is.numeric(y.factors[i,f]) || (abs(y.factors[i,f]) > cut)) 
                    if((!is.numeric(y.factors[i,f] ) && (y.factors[i,f] !="0"))||  ((is.numeric(y.factors[i,f]) && abs(y.factors[i,f]) > cut ))) {
                    
                        
                           sem[k,1] <- paste(fact[f+num.xfactors],"->",vars[i+num.xvar],sep="")
                           if(is.numeric(y.factors[i,f])) { sem[k,2] <- paste("Fy",f,vars[i+num.xvar],sep="")} else {sem[k,2] <- y.factors[i,f]}
                         k <- k+1 }   #end of if 
                         }  #end of factor
                        } # end of variable loop
           
       }  
       if(errors) {  for (i in 1:num.y) {
                    sem[k,1] <- paste(vars[i+num.xvar],"<->",vars[i+num.xvar],sep="")
                    sem[k,2] <- paste("y",i,"e",sep="")
                    k <- k+1 }}
  
  }   #end of !is.null(ymodel)
 
 if (!is.null(labels)) {var.labels <- c(labels,fact)
  } 

if(!regression) {
  if(!is.null(Phi)) {if (!is.matrix(Phi))  Phi <- matrix(c(1,Phi,0,1),ncol=2) 
  if(num.xfactors>1) {for (i in 2:num.xfactors) {
                         for (j in 1:(i-1)) {if((!is.numeric(Phi[i,j] ) && (Phi[i,j] !="0"))||  ((is.numeric(Phi[i,j]) && abs(Phi[i,j]) > cut ))) {
                          
                            if(Phi[i,j] != Phi[j,i]){
                            sem[k,1]  <- paste(fact[j],"->",fact[i],sep="")
                            if (is.numeric(Phi[i,j])) {sem[k,2] <-  paste("rF",j,"F",i,sep="")} else {sem[k,2] <- Phi[i,j] }
                            } else {
                            sem[k,1]  <- paste(fact[i],"<->",fact[j],sep="")
                             if (is.numeric(Phi[i,j])) {sem[k,2] <-  paste("rF",i,"F",j,sep="")} else {sem[k,2] <- Phi[i,j] }
                           }
                          k <- k + 1
                              }              }
                            }
                       }
                       } #end of correlations within x set
               if(!is.null(ymodel)) {
                  for (i in 1:num.xfactors) { 
                      for (j in 1:num.yfactors) {
                      if((!is.numeric(Phi[i,j+num.xfactors] ) && (Phi[i,j+num.xfactors] !="0"))||  ((is.numeric(Phi[i,j+num.xfactors]) && abs(Phi[i,j+num.xfactors]) > cut ))) {
                      # clust.graph <- addEdge( fact[j+num.xfactors],fact[i],clust.graph,1)
                      # if (is.numeric(Phi[i,j+num.xfactors])) { edge.label[k] <- round(Phi[i,j+num.xfactors],digits)} else {edge.label[k] <- Phi[i,j+num.xfactors]}
                      # edge.dir[k] <- "back"
                      # edge.name[k] <- paste(fact[j+num.xfactors],"~",fact[i],sep="")
                        sem[k,1]  <- paste(fact[i],"->",fact[j+num.xfactors],sep="")
                        if (is.numeric(Phi[i,j+num.xfactors])) {sem[k,2] <-  paste("rX",i,"Y",j,sep="")} else {sem[k,2] <- Phi[i,j+num.xfactors] } 
                       k <- k + 1 }
                        }
                      
                  
                  }
                          
      }  
     } else {if(!is.null(Phi)) {if (!is.matrix(Phi))  Phi <- matrix(c(1,Phi,0,1),ncol=2) 
        for (i in 2:num.xvar) {
             for (j in 1:(i-1))  {

                if(Phi[i,j] != Phi[j,i]){edge.dir[k] <- "back"} else {edge.dir[k] <- "both"}
               k <- k + 1            }}
                              }          
     }
     
     for(f in 1:num.factors) {
        sem[k,1]  <- paste(fact[f],"<->",fact[f],sep="")
        sem[k,3] <-  "1"
         k <- k+1
     }
  
model=sem[1:(k-1),]
class(model) <- "mod"  #to make for pretty output when using sem package -- suggested by John Fox
return(model)
   }
