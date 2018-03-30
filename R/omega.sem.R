#January 22, 2008
#omega.graph without the graph

"omega.sem" <-
function(om.results,out.file=NULL,sl=TRUE,labels=NULL,nf=3){
if(is.list(om.results)) {
   if (sl) {factors <- as.matrix(om.results$schmid$sl)   } else {factors <- as.matrix(om.results$schmid$oblique)}
   
  #first some basic setup parameters 
  
   num.var <- dim(factors)[1]   #how many variables?
  if (sl) {num.factors <- dim(factors)[2] -4 } else {num.factors <- dim(factors)[2]}  #g, h2,u2,p2
 # if(num.factors ==1) {warning("giving lavaan code for a 1 factor omega doesn't really make sense.")}
  # return(list(sem=NULL,lavaan=NULL))}
   gloading <- om.results$schmid$gloading
   
    } else {factors <- om.results
    num.var <- nrow(factors)
    gloading <- factors[,1]
    num.factors <-nf+1}
    
    if(sl) {fact <- c("g",paste("F",1:num.factors,"*",sep="")) } else {fact <- c("g",paste("F",1:num.factors,sep="")) }   # e.g.  "g"  "F'1" "F2" "F3"
      vars <- paste("V",1:num.var,sep="")  
   if (!is.null(labels)) {vars <- paste(labels)} else{vars <- rownames(factors) }
   
   lavaan <- vector("list",nf+1)
   if (sl) {
            sem <- matrix(rep(NA,6*(2*num.var + num.factors)),ncol=3)  #used for sem diagram
            
            } else {
            
            sem <- matrix(rep(NA,6*(num.var + num.factors)+3),ncol=3)  #used for sem diagram
            
                       }
            
  #show the cluster structure with ellipses
   if (sl) {
   l <- matrix(factors[,2:(num.factors+1)],ncol=num.factors) } else { l <- factors }
   
   m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max),  ncol = 1)
     lavaan[[1]] <- 'g =~ '    
  if (sl) { 
          k <- num.var
         
         for (i in 1:num.var) {  
                                 sem[i,1] <- paste(fact[1],"->",vars[i],sep="")
                                 sem[i,2] <- vars[i]
                                 lavaan[[1]] <- paste0(lavaan[[1]], '+', vars[i])               
            } 
            
         } else { 
          
          k <- num.factors
          for (j in 1:num.factors) {
           sem[j,1] <- paste(fact[1],"->",fact[1+j],sep="")
           sem[j,2] <- paste("g",fact[1+j],sep="")
        lavaan[[1]] <- paste0(lavaan[[1]], ' + ', fact[1+j]) 
        
          
                 } 
                 }
  if(num.factors > 1)  for(j in 1:num.factors) { lavaan[[1+j]] <- paste0('F',j ,"=~ ")}
   

   for (i in 1:num.var) { 
                         sem[i+k,1] <- paste(fact[1+m1[i]],"->",vars[i],sep="")
                        sem[i+k,2] <- paste(fact[1+m1[i]],vars[i],sep="")
                        if(num.factors > 1) lavaan[[1+ m1[i ]]] <- paste0(lavaan[[1 + m1[i] ]] ,' + ', vars[i])
                        } 
                    

 if(sl) {
       k <- num.var*2
      for (i in 1:num.var) {
            sem[i+k,1] <- paste(vars[i],"<->",vars[i],sep="")
            sem[i+k,2] <- paste("e",i,sep="")
                            }
        k <- k + num.var
       for (f in 1:num.factors) {
             sem[f+k,1] <- paste(fact[1+f],"<->",fact[1+f],sep="")
             sem[f+k,3] <- "1"
                                 }
        k <- k+ num.factors
           sem[k+1,1] <- paste("g <->g")
           sem[k+1,3] <- "1"
           k<- k+1
           
          } else { 
          
          k <- num.var + num.factors
          for (i in 1:num.var) {
            sem[i+k,1] <- paste(vars[i],"<->",vars[i],sep="")
            sem[i+k,2] <- paste("e",i,sep="")
                            }
            k <- 2*num.var + num.factors
          for (f in 1:num.factors) {
            sem[f+k,1] <- paste(fact[f+1],"<->",fact[f+1],sep="")
            sem[f+k,3] <- "1"
                            }
             k <- 2*num.var + 2*num.factors
             sem[k+1,1] <- paste("g<->g")
             sem[k+1,3] <- "1"
             k <- k+1
          }
          
colnames(sem) <- c("Path","Parameter","Initial Value")
sem[,1] <- gsub("-","",sem[,1]) #get rid of the negative signs for variables
sem[,1] <- gsub(">","->",sem[,1])  #but put back in the arrows
lavaan <- gsub("-","",unlist(lavaan))
lavaan <- noquote(lavaan)
return(list(sem=sem[1:k,],lavaan=lavaan))
#return(list(sem=sem[1:k,],)
   }
 
