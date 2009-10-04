#Created September 25, 2009
#based upon structure graph but not using Rgraphviz
#creates a structural equation path diagram, draws it, and saves sem commands

"structure.diagram" <-
function(fx,Phi=NULL,fy=NULL,labels=NULL,cut=.3,errors=FALSE,simple=TRUE,regression=FALSE,
   digits=1,e.size=.1,main="Structural model", ...){
    	
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
 
  #first some basic setup parameters    -- these just convert the various types of input 
 if(!is.matrix(xmodel) ) {factors <- as.matrix(xmodel)} else {factors <- xmodel}
 
  num.y <- 0   #we assume there is nothing there
  num.var <- num.xvar <- dim(factors)[1]   #how many x variables?
  
  if (is.null(num.xvar) ){num.xvar <- length(factors)
                          num.xfactors <- 1} else {
  						 num.factors <-  num.xfactors <- dim(factors)[2]}
  	
  if(is.null(labels)) {vars <- xvars <-  rownames(xmodel)} else { xvars <-  vars <- labels}
  
  if(is.null(vars) ) {vars <- xvars <- paste("x",1:num.xvar,sep="")  }
  fact <- colnames(xmodel)
  if (is.null(fact)) { fact <- paste("X",1:num.xfactors,sep="") }
   if(is.numeric(factors)) {factors <- round(factors,digits)
                            }
   
   num.yfactors <- 0
   num.yvar <- 0   
  
   if (!is.null(ymodel)) { 
   		if(is.list(ymodel) & !is.data.frame(ymodel)  ) {ymodel <- as.matrix(ymodel$loadings)} else {ymodel <- ymodel}   
	 	if(!is.matrix(ymodel) ) {y.factors <- as.matrix(ymodel)} else {y.factors <- ymodel}
      num.y <- dim(y.factors)[1] 
      if (is.null(num.y)) {
         num.y <- length(ymodel)
         num.yfactors <- 1} else {
         num.yfactors <- dim(y.factors)[2]
                          }  
     num.yvar <- num.y
     yvars <- rownames(ymodel)
     if(is.null(yvars)) {yvars <- paste("y",1:num.y,sep="")  }
     if(is.null(labels)) {vars <- c(xvars,yvars)} else {yvars <- labels[(num.xvar+1):(num.xvar+num.y)]} 
     
      yfact <- colnames(ymodel)
      if(is.null(yfact)) {yfact <- paste("Y",1:num.yfactors,sep="") }
      fact <- c(fact,yfact)
      if(is.numeric(y.factors)) {y.factors <- round(y.factors,digits)
                                 }
   }#end of if(null(y.model))
   
    num.var <- num.xvar + num.y
     num.factors <- num.xfactors + num.yfactors
    
    sem <- matrix(rep(NA),6*(num.var*num.factors + num.factors),ncol=3)    #this creates an output model for sem analysis
    colnames(sem) <- c("Path","Parameter","Value")
    var.rect <- list()
    fact.rect <- list()
    
    if(is.numeric(Phi) ) { Phi <- round(Phi,digits)}
                           

 
###create the basic figure

length.labels <-  0    # a filler for now
limx <- c(-length.labels,max(num.xvar,num.yvar)+2)
limy <-  c(0,max(num.xvar,num.yvar)+1)


plot(0,type="n",xlim=limx,ylim=limy,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
    
   #now draw the x part
   
   k <- num.factors  
   
    for (v in 1:num.xvar) {   
 	 var.rect[[v]] <- dia.rect(0,num.xvar-v+1,xvars[v],xlim=limx,ylim=limy,...)
     }
  nvar <- num.xvar
  f.scale <- (num.xvar+ 1)/(num.xfactors+1)
 for (f in 1:num.xfactors) {
   		fact.rect[[f]] <- dia.ellipse(limx[2]/3,(num.xfactors+1-f)*f.scale,fact[f],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
     		for (v in 1:num.xvar)  {
     		    if(is.numeric(factors[v,f])) {
     		    if(simple && (abs(factors[v,f]) == max(abs(factors[v,])) )  && (abs(factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) { if (!regression) {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =factors[v,f],col=((sign(factors[v,f])<0) +1)) 
    			
    		} else {dia.arrow(to=fact.rect[[f]],from=var.rect[[v]]$right,labels =factors[v,f],col=((sign(factors[v,f])<0) +1))} }
                                    }  else {
                if (factors[v,f] !="0") {
               if (!regression) { dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =factors[v,f]) 
    			} else {dia.arrow(to=fact.rect[[f]],from=var.rect[[v]]$right,labels =factors[v,f])} }
                                    } }
                              } 
                              
                              
            if (num.xfactors ==1) { 
                      for(i in 1:num.xvar) {
                       sem[i,1] <- paste(fact[1],"->",vars[i],sep="")
                          if(is.numeric(factors[i])) {sem[i,2] <- vars[i]} else {sem[i,2] <- factors[i] }
                        }}  #end of if num.xfactors ==1 
       k <- num.xvar+1 
  
                   k <- 1
                   for (i in 1:num.xvar) {
                   for (f in 1:num.xfactors) { #if (!is.numeric(factors[i,f]) ||  (abs(factors[i,f]) > cut))
                   if((!is.numeric(factors[i,f] ) && (factors[i,f] !="0"))||  ((is.numeric(factors[i,f]) && abs(factors[i,f]) > cut ))) {
                              
                              sem[k,1] <- paste(fact[f],"->",vars[i],sep="")
                             if(is.numeric(factors[i,f])) {sem[k,2] <- paste("F",f,vars[i],sep="")} else {sem[k,2] <- factors[i,f]}
                              k <- k+1 }   #end of if 
                         }  
                        }      
        
  if(errors) {  for (i in 1:num.xvar) { dia.self(var.rect[[i]],side=3) 
                                          sem[k,1] <- paste(vars[i],"<->",vars[i],sep="")
                                          sem[k,2] <- paste("x",i,"e",sep="")
                                    k <- k+1 }
                    }  
     
    
       
 ## #now, if there is a ymodel, do it for y model 
  
  if(!is.null(ymodel)) { 
  	 y.adj <-  num.yvar - num.xvar +1
   	 f.yscale <- limy[2]/(num.yfactors+1)
   	 y.fadj <- 0
   for (v in 1:num.yvar) {    var.rect[[v+num.xvar]] <- dia.rect(limx[2],limy[2]-v + y.adj,yvars[v],xlim=limx,ylim=limy,...) }
  
   for (f in 1:num.yfactors) {
   		fact.rect[[f+num.xfactors]] <- dia.ellipse(2*limx[2]/3,(num.yfactors+1-f)*f.yscale +y.fadj,yfact[f],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
   		  
     		for (v in 1:num.yvar) {if(is.numeric(y.factors[v,f])) {
     		{if(simple && (abs(y.factors[v,f]) == max(abs(y.factors[v,])) )  && (abs(y.factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) {dia.arrow(from=fact.rect[[f+num.xfactors]],to=var.rect[[v+num.xvar]]$left,labels =y.factors[v,f],col=((sign(y.factors[v,f])<0) +1)) }
                                    }
                   } else {if(factors[v,f] !="0")  {dia.arrow(from=fact.rect[[f+num.xfactors]],to=var.rect[[v+num.xvar]]$left,labels =y.factors[v,f]) }
             }}
                             } 
                       
  
  	if (num.yfactors ==1) {     
    for (i in 1:num.y) { sem[k,1] <- paste(fact[1+num.xfactors],"->",yvars[i],sep="")
                           if(is.numeric(y.factors[i] ) ) {sem[k,2] <- paste("Fy",yvars[i],sep="")} else {sem[k,2] <- y.factors[i]}
                         k <- k +1
                        }                      
                         } else {   #end of if num.yfactors ==1 
        
                   for (i in 1:num.y) {
                   for (f in 1:num.yfactors) { 
                    if( (y.factors[i,f] !="0")  && (abs(y.factors[i,f]) > cut )) {
                           sem[k,1] <- paste(fact[f+num.xfactors],"->",vars[i+num.xvar],sep="")
                           if(is.numeric(y.factors[i,f])) { sem[k,2] <- paste("Fy",f,vars[i+num.xvar],sep="")} else {sem[k,2] <- y.factors[i,f]}
                         k <- k+1 }   #end of if 
                         }  #end of factor
                        } # end of variable loop      
                               } #end of else if 
                  # }
       
     if(errors) {  for (i in 1:num.y) { 
                      dia.self(var.rect[[i+num.xvar]],side=3) 
                    sem[k,1] <- paste(vars[i+num.xvar],"<->",vars[i+num.xvar],sep="")
                    sem[k,2] <- paste("y",i,"e",sep="")
                    k <- k+1 }}
  
    } #end of if.null(ymodel)

if(!regression) {
  if(!is.null(Phi)) {if (!is.matrix(Phi)) { if(!is.null(fy)) {Phi <- matrix(c(1,0,Phi,1),ncol=2)} else {Phi <- matrix(c(1,Phi,Phi,1),ncol=2)}} 
 
 if(num.xfactors>1) {for (i in 2:num.xfactors) { #first do the correlations within the f set
      for (j in 1:(i-1)) {
                 {if((!is.numeric(Phi[i,j] ) && ((Phi[i,j] !="0")||(Phi[j,i] !="0")))||  ((is.numeric(Phi[i,j]) && abs(Phi[i,j]) > cut ))) {
                           
                          if(Phi[i,j] == Phi[j,i] ) {
                                                      	dia.curve(from=fact.rect[[i]]$right,to=fact.rect[[j]]$right, labels = Phi[i,j],scale=2*(i-j)/num.xvar)
                                                     	sem[k,1]  <- paste(fact[i],"<->",fact[j],sep="")
                                                     	sem[k,2] <-  paste("rF",i,"F",j,sep="")} else {
                                                     	
                                                     	    if(Phi[i,j] !="0") { 	
                                                      		sem[k,1]  <- paste(fact[i]," ->",fact[j],sep="")
                                                      		sem[k,2] <-  paste("rF",i,"F",j,sep="")} else {
                                                      
                                                      sem[k,1]  <- paste(fact[i],"<-",fact[j],sep="")
                                                      sem[k,2] <-  paste("rF",i,"F",j,sep="")} 
                                                     	} } else { 
                           
                            sem[k,1]  <- paste(fact[i],"<->",fact[j],sep="")
                            
                             if (is.numeric(Phi[i,j])) {sem[k,2] <-  paste("rF",i,"F",j,sep="")} else {sem[k,2] <- Phi[i,j] } }
                            
                                   
                            k <- k + 1} }
                            
                            }
                            }  #end of correlations within the fx set
      if(!is.null(ymodel)) {
                  for (i in 1:num.xfactors) { 
                      for (j in 1:num.yfactors) {
                      		if((!is.numeric(Phi[j+num.xfactors,i] ) && (Phi[j+num.xfactors,i] !="0"))||  ((is.numeric(Phi[j+num.xfactors,i]) && abs(Phi[j+num.xfactors,i]) > cut ))) {
                      
                       dia.arrow(from=fact.rect[[i]],to=fact.rect[[j+num.xfactors]],Phi[j+num.xfactors,i])
                       sem[k,1]  <- paste(fact[i],"->",fact[j+num.xfactors],sep="") } else {
                      
                       sem[k,1]  <- paste(fact[i],"<->",fact[j+num.xfactors],sep="")}
                       
                        if (is.numeric(Phi[j+num.xfactors,i])) {sem[k,2] <-  paste("rX",i,"Y",j,sep="")} else {sem[k,2] <- Phi[j+num.xfactors,i] } 
                       k <- k + 1 }
                        } }
               } }          
     for(f in 1:num.factors) {
       	 	sem[k,1]  <- paste(fact[f],"<->",fact[f],sep="")
        	sem[k,3] <-  "1"
        	 k <- k+1
     		} 
model=sem[1:(k-1),]
class(model) <- "mod"   #suggested by John Fox to make the output cleaner
return(model)
   }
   
 
